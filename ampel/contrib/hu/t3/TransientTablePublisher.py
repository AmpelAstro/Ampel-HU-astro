#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/TransientTablePublisher.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                06.05.2021
# Last Modified Date:  05.12.2021
# Last Modified By:    jnordin@physik.hu-berlin.de

import re, os, requests, io, backoff
import pandas as pd

from functools import reduce
from typing import Any, Optional, Generator

from ampel.types import T3Send
from ampel.secret.NamedSecret import NamedSecret
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit
from ampel.view.T3Store import T3Store
from ampel.view.SnapView import SnapView


class TransientTablePublisher(AbsT3ReviewUnit):
    """

    Construct a table based on selected T2 output values.
    Current output format can be csv or latex.
    Table can optionally saved to a local file or submitted to slack (future: posted to desy web).

    Config parameters:
    include_stock (bool)
    include_channels (bool)

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
    - include format option for prointing

    """

    # Two tables describing what information to save into the table.
    # Schema for state dependent T2s (one row for each)
    table_schema: dict['str', Any]
    # Schema for transient dependent T2s (added to each row together with base info)
    transient_table_schema: dict['str', Any]

    name_filter: dict[str, str] = {'ZTF name': 'ZTF', 'TNS ID': 'TNS'}
    include_stock: bool = False
    include_channels: bool = True
    # Add also transients lacking any T2 info
    save_base_info: bool = False

    fmt: str = 'csv'

    file_name: str = 'TransientTable.csv'
    slack_channel: Optional[str] = None
    slack_token: Optional[NamedSecret[str]]
    local_path: Optional[str] = None


    def process(self, gen: Generator[SnapView, T3Send, None], t3s: T3Store) -> None:
        """
        Loop through provided TransientViews and extract data according to the
        configured schema.
        """


        table_rows: list[dict[str, Any]] = []
        for k, tran_view in enumerate(gen, 1):
            basetdict = {}
            # Assemble t2 information bound to the transient (e.g. Point T2s)
            for t2unit, table_entries in self.transient_table_schema.items():
                # New SnapView has method for directly retrieve result.
                # Possibly use this.
                if t2res := tran_view.get_latest_t2_body(unit_id=t2unit):
                    for label, path in table_entries.items():
                        try:
                            basetdict[label] = reduce(dict.get, path, t2res)
                        except (KeyError, TypeError):
                            pass # We leave missing entries empty

            # Assemble info which could vary from state to state
            # Should add config to labels if multiple exports
            # from same unit is requested.
            stateinfo = []
            for t1_document in tran_view.t1 or []:
                t1_link = t1_document['link']
                tdict = {}
                for t2unit, table_entries in self.table_schema.items():
                    if t2res := tran_view.get_latest_t2_body(unit=t2unit, link=t1_link):
                        for label, path in table_entries.items():
                            try:
                                tdict[label] = reduce(dict.get, path, t2res)
                            except (KeyError, TypeError):
                                pass # We leave missing entries empty
                if len(tdict) > 0:
                    stateinfo.append(tdict)

            if len(stateinfo) == 0 and len(basetdict.keys()) == 0 and not self.save_base_info:
                continue

            # Collect base information applying to all states
            # If here, add stock info (name, channel etcs)
            if names := (tran_view.stock or {}).get("name"):
                for label, name_str in self.name_filter.items():
                    r = re.compile(name_str)
                    # While names will mostly be unique, it might not always be the case.
                    basetdict[label] = list(filter(r.match, names))
                    # Avoid list when possible
                    if len(basetdict[label]) == 1:
                        basetdict[label] = basetdict[label][0]

            if self.include_stock:
                basetdict['stock'] = tran_view.id
            if self.include_channels:
                basetdict['channels'] = tran_view.stock.get("channel")
            # Allow for both single (most common) and duplacte channels.
                if len(basetdict['channels']) == 0:
                    basetdict['channels'] = basetdict['channels'][0]

            # Collect and add to table
            if len(stateinfo) > 0:
                for tdict in stateinfo:
                    tdict.update(basetdict)
                    table_rows.append(tdict)
            else:
                # Only transient info
                table_rows.append(basetdict)

        self.logger.info("", extra={'table_count': len(table_rows)})
        if len(table_rows) == 0:
            return None

        # Export assembled information
        # Convert
        df = pd.DataFrame.from_dict(table_rows)

        # Local save
        if self.local_path is not None:
            full_path = os.path.join(self.local_path, self.file_name)
            if self.fmt == 'csv':
                df.to_csv(full_path)
            elif self.fmt == 'latex':
                df.to_latex(full_path)
            self.logger.info('Exported', extra={'path': full_path})

        # Export to slack if requested
        if self.slack_channel is not None and self.slack_token is not None:
            self._slack_export(df)

        # Could potentially return a document to T3 collection detailing
        # what was done, as well as the table itself.
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
        giveup=lambda e: e.response.status_code not in {503, 429},
        max_time=60,
    )
    def _slack_export(self, df):
        """
        Export content of Pandas dataframe to slack.
        """

        # Slack summary
        buffer = io.StringIO(self.file_name)
        if self.fmt == 'csv':
            df.to_csv(buffer)
        elif self.fmt == 'latex':
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
