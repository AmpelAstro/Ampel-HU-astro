#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/TransientTablePublisher.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                06.05.2021
# Last Modified Date:  16.01.2024
# Last Modified By:    ernstand@physik.hu-berlin.de

import io
import os
import re
import numpy as np
from collections.abc import Generator
from functools import reduce
from typing import Any, Optional, Union

import json

import backoff
import pandas as pd
import requests
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.util.mappings import get_by_path
from ampel.view.SnapView import SnapView
from ampel.view.TransientView import TransientView


class TransientTablePublisher(AbsPhotoT3Unit):
    """

    Construct a table based on selected T2 output values.
    Current output format can be csv or latex.
    Table can optionally saved to a local file or submitted to slack (future: posted to desy web).

    Config parameters:
    include_stock (bool)
    include_channels (bool)

    If one wants to convert the AMPEL stock ID to external IDs, define
    convert_stock_to (str|None)
    For ZTF-IDs, pass 'convert_stock_to: ztf'

    How to deal with names. Will search each transients names for entries containing "value",
    and return any output under "key"
    name_filter = { 'ZTF name' : 'ZTF', 'TNS ID' : 'TNS' }

    Selection of fields to save. Matches structure of t2document result dict, e.g.:
    table_schema = { { 't2_unit'  : {
                'table_label_1' : ['path','to','val'],
            'table_label_2' : ['other','path']
            },
        } }
    transient_table_schema = { { 'point_t2_unit'  : {
                'table_label_1' : ['path','to','val'],
            'table_label_2' : ['other','path']
            },
        } }

    Output format (converted through pandas)
    fmt = 'csv'     # Current options 'csv', 'latex'.

    Destination attempted if the appropriate  parameters are set for
    file_name
    slack:
      slack_channel
      slack_token
    local save:
      local_path

    Todo:
    - save to desy webb?
    - include format option for printing

    """

    # Two tables describing what information to save into the table.
    # Schema for state dependent T2s (one row for each)
    table_schema: dict[str, Any] = {}
    # Schema for transient dependent T2s (added to each row together with base info)
    transient_table_schema: dict[str, Any]

    name_filter: dict[str, str] = {"ZTF name": "ZTF", "TNS ID": "TNS"}
    include_stock: bool = False
    convert_stock_to: str | None = None

    sort_by_key: str | None = "kilonovaness"
    sort_ascending: bool = False

    include_pos: bool = True
    include_channels: bool = True
    # Add also transients lacking any T2 info
    save_base_info: bool = False

    fmt: str = "csv"
    write_mode: str = "a"
    rename_files: bool = False

    dir_name: str = "TransientTable"
    file_name: str = dir_name
    slack_channel: None | str = None
    slack_token: None | NamedSecret[str]
    local_path: None | str = None

    move_files: bool = True

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        #    def process(self, gen: Generator[SnapView, T3Send, None], t3s: T3Store) -> None:
        """
        Loop through provided TransientViews and extract data according to the
        configured schema.
        """

        table_rows: list[dict[str, Any]] = []
        for k, tran_view in enumerate(gen, 1):
            basetdict: dict[str, Any] = {}
            # Assemble t2 information bound to the transient (e.g. Point T2s)
            for t2unit, table_entries in self.transient_table_schema.items():
                # New SnapView has method for directly retrieve result.
                # Possibly use this.
                if isinstance(t2res := tran_view.get_latest_t2_body(unit=t2unit), dict):
                    for label, path in table_entries.items():
                        basetdict[label] = get_by_path(t2res, path)

            # Assemble info which could vary from state to state
            # Should add config to labels if multiple exports
            # from same unit is requested.
            stateinfo = []
            for t1_document in tran_view.t1 or []:
                t1_link = t1_document["link"]
                tdict = {}
                for t2unit, table_entries in self.table_schema.items():
                    if isinstance(
                        t2res := tran_view.get_latest_t2_body(
                            unit=t2unit, link=t1_link
                        ),
                        dict,
                    ):
                        for label, path in table_entries.items():
                            tdict[label] = get_by_path(t2res, path)
                if len(tdict) > 0:
                    stateinfo.append(tdict)

            if (
                len(stateinfo) == 0
                and len(basetdict.keys()) == 0
                and not self.save_base_info
            ):
                continue

            # Collect base information applying to all states
            # If here, add stock info (name, channel etcs)
            if names := (tran_view.stock or {}).get("name", []):
                for label, name_str in self.name_filter.items():
                    r = re.compile(name_str)
                    # While names will mostly be unique, it might not always be the case.
                    basetdict[label] = list(filter(r.match, names))  # type: ignore[arg-type]
                    # Avoid list when possible
                    if (
                        isinstance((item := basetdict[label]), (list, tuple))
                        and len(item) == 1
                    ):
                        basetdict[label] = item[0]

            if self.include_stock:
                basetdict["stock"] = tran_view.id

            if self.convert_stock_to is not None:
                assert self.convert_stock_to in ["ztf"]

                if self.convert_stock_to == "ztf":
                    from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper

                    stock_id = tran_view.id
                    ztf_id = ZTFIdMapper.to_ext_id(stock_id)
                    basetdict["ztf_id"] = ztf_id

            if self.include_pos:
                lcurve = tran_view.get_lightcurves()
                if lcurve is not None:
                    pos = lcurve[0].get_pos(ret="brightest")
                    if pos is not None:
                        basetdict["ra"] = pos[0]
                        basetdict["dec"] = pos[1]
            if self.include_channels and tran_view.stock:
                channels = tran_view.stock.get("channel")
                # Allow for both single (most common) and duplacte channels.
                basetdict["channels"] = (
                    channels[0]
                    if isinstance(channels, (list, tuple)) and len(channels) == 1
                    else channels
                )

            # Collect and add to table
            if len(stateinfo) > 0:
                for tdict in stateinfo:
                    tdict.update(basetdict)
                    table_rows.append(tdict)
            else:
                # Only transient info
                table_rows.append(basetdict)

        self.logger.info("", extra={"table_count": len(table_rows)})
        if len(table_rows) == 0:
            return None

        # Export assembled information
        # Convert
        df = pd.DataFrame.from_dict(table_rows)

        # if "map_name" in df.columns and "map_seed" in df.columns:
        # df["map_name"] = np.char.replace(np.array(df["map_name"], dtype=str), "random", "random"+df["map_seed"])

        # print(df["map_name"].iloc[0])
        if "random" in df["map_name"].iloc[0] or self.rename_files:
            # print("transienttablepublisher:: ", df["map_seed"].iloc[0])
            tmp_seed_name = df["map_seed"].iloc[0]
            if type(tmp_seed_name) == str:
                self.file_name += "_" + tmp_seed_name
            else:
                self.file_name += "_" + str(int(tmp_seed_name))

        # sort dataframe by key
        if self.sort_by_key in df.keys():
            df = df.sort_values(by=self.sort_by_key, ascending=self.sort_ascending)
        else:
            self.logger.warn(
                f"Cannot sort table by {self.sort_by_key} - legal keys: {df.keys()}"
            )

        # Local save
        if self.local_path is not None:
            path_name = os.path.join(self.local_path, self.dir_name)
            # print("PATHNAME::", path_name)
            if not os.path.exists(path_name):
                os.makedirs(path_name, exist_ok=True)
            full_path = os.path.join(path_name, self.file_name)
            # print("FILE PATH::", full_path)

            with open(full_path + "." + self.fmt, "w") as tmp_file:
                tmp_file.close()
            if self.fmt == "csv":
                # print(self.write_mode)
                df.to_csv(full_path + ".csv", sep=";", mode=self.write_mode)
            elif self.fmt == "latex":
                df.to_latex(full_path + ".tex")
            elif self.fmt == "json":
                json_dumps = df.to_json(indent=2)
                with open(full_path + ".json", self.write_mode) as json_file:
                    json_file.write(json_dumps)
                    json_file.close()
            self.logger.info("Exported", extra={"path": full_path})

        # Export to slack if requested
        self._slack_export(df)

        # Could potentially return a document to T3 collection detailing
        # what was done, as well as the table itself.

        # take everything local_path and put it into new folder named after skymap
        # print(df.keys)
        map_name_key = "map_name"
        if map_name_key in list(df.keys()) and self.move_files:
            files_local_path = os.listdir(self.local_path)
            skymap_name = df[map_name_key][
                0
            ]  # need to change if for some reason several maps get saved in same file
            skymap_dir_name = skymap_name  # [: skymap_name.find(".")]  # bare name
            if skymap_name[-1] != "z":  # if non trivial rev version (hacky)
                skymap_dir_name += (
                    "_rev_" + skymap_name[skymap_name.find(",") + 1 :]
                )  # find "," and add rev version after that

            print("TransientTablePublisher: TMP FILES MOVED TO " + skymap_dir_name)

            if self.local_path is not None:
                skymap_directory = os.path.join(
                    self.local_path + "/../" + skymap_dir_name
                )
                # print(skymap_directory)
                os.makedirs(skymap_directory, exist_ok=True)
                for file in files_local_path:
                    if file.find(".fits.gz") == -1:
                        tmp_file_path = os.path.join(self.local_path, file)
                        if not (os.path.isfile(tmp_file_path)):
                            continue
                        os.replace(tmp_file_path, os.path.join(skymap_directory, file))

        return None

    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: not isinstance(e, requests.HTTPError)
        or e.response.status_code not in {503, 429},
        max_time=60,
    )
    def _slack_export(self, df):
        """
        Export content of Pandas dataframe to slack.
        """
        if self.slack_channel is None or self.slack_token is None:
            return

        # Slack summary
        buffer = io.StringIO(self.file_name)
        if self.fmt == "csv":
            df.to_csv(buffer, sep=";")
        elif self.fmt == "latex":
            df.to_latex(buffer)

        param = {
            "token": self.slack_token.get(),
            "channels": self.slack_channel,
            "title": "From the Table Publisher",
            "username": "AMPEL-live",
            "as_user": "false",
            "filename": self.file_name,
        }

        ret = requests.post(
            "https://slack.com/api/files.upload",
            params=param,
            files={"file": buffer.getvalue()},
        )
        ret.raise_for_status()

        self.logger.info(ret.text)

        return
