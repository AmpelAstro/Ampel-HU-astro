#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/SlackSummaryPublisher.py
# License           : BSD-3-Clause
# Author            : robert stein
# Date              : 11.03.2018
# Last Modified Date: 14.11.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import collections
import datetime
import io
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import requests
from slack import WebClient
from slack.errors import SlackClientError

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.log.utils import log_exception
from ampel.model.Secret import Secret
from ampel.util.collections import ampel_iter
from ampel.view.TransientView import TransientView
from ampel.ztf.utils import to_ztf_id


class SlackSummaryPublisher(AbsT3Unit):
    """"""

    dry_run: bool = False
    quiet: bool = False
    slack_channel: str
    slack_token: Secret
    excitement: Dict[str, int] = {"Low": 50, "Mid": 200, "High": 400}
    full_photometry: bool = False
    #: Fields to extract from each detection (if full_photometry enabled)
    cols: List[str] = [
        "ra",
        "dec",
        "magpsf",
        "sgscore1",
        "rb",
        "distnr",
        "distpsnr1",
        "isdiffpos",
    ]

    def post_init(self) -> None:
        """ """
        self.frames: List[pd.DataFrame] = []
        self.photometry: List[pd.DataFrame] = []
        self.channels: Set[str] = set()

    def add(self, transients: Tuple[TransientView, ...]) -> None:
        """"""
        summary, full = self.combine_transients(transients)
        self.frames += summary
        self.photometry += full

    def done(self):
        """"""
        if len(self.frames) == 0 and self.quiet:
            return

        date = str(datetime.date.today())

        sc = WebClient(self.slack_token.get())

        m = calculate_excitement(
            len(self.frames), date=date, thresholds=self.excitement
        )

        if self.dry_run:
            self.logger.info(m)
        else:
            api = sc.api_call(
                "chat.postMessage",
                json={
                    "channel": self.slack_channel,
                    "text": m,
                    "username": "AMPEL-live",
                    "as_user": False,
                },
            )
            if not api["ok"]:
                raise SlackClientError(api["error"])

        if len(self.frames) > 0:

            df = pd.DataFrame.from_records(self.frames)
            # Set fill value for channel columns to False
            for channel in self.channels:
                df[channel].fillna(False, inplace=True)
            # Set fill value for all other columns to MISSING
            for field in set(df.columns.values).difference(self.channels):
                df[field].fillna("MISSING", inplace=True)
            # Move channel info at end
            df = df.reindex(
                copy=False,
                columns=[c for c in df.columns if not c in self.channels]
                + list(self.channels),
            )

            filename = "Summary_%s.csv" % date

            buffer = io.StringIO(filename)
            df.to_csv(buffer)

            param = {
                "token": self.slack_token.get(),
                "channels": self.slack_channel,
                "title": "Summary: " + date,
                "username": "AMPEL-live",
                "as_user": "false",
                "filename": filename,
            }

            if self.dry_run:
                # log only first two lines
                csv = buffer.getvalue()
                idx = 0
                for _ in range(2):
                    idx = csv.find("\n", idx) + 1
                self.logger.info(
                    {
                        "files": {"file": csv[:idx] + "..."},
                        "token": param["token"][:8] + "...",
                        **{k: v for (k, v) in param.items() if k != "token"},
                    }
                )
            else:
                r = requests.post(
                    "https://slack.com/api/files.upload",
                    params=param,
                    files={"file": buffer.getvalue()},
                )
                self.logger.info(r.text)

            if self.full_photometry:
                photometry = pd.concat(self.photometry, sort=False)
                # Set fill value for channel columns to False
                for channel in self.channels:
                    photometry[channel].fillna(False, inplace=True)
                # Set fill value for all other columns to MISSING
                for field in set(photometry.columns.values).difference(self.channels):
                    photometry[field].fillna("MISSING", inplace=True)
                # Move channel info at end
                photometry = photometry.reindex(
                    copy=False,
                    columns=[c for c in df.columns if not c in self.channels]
                    + list(self.channels),
                )

                filename = "Photometry_%s.csv" % date

                buffer = io.StringIO(filename)
                photometry.to_csv(buffer)

                param = {
                    "token": self.slack_token.get(),
                    "channels": self.slack_channel,
                    "title": "Full Photometry: " + date,
                    "username": "AMPEL-live",
                    "as_user": "false",
                    "filename": filename,
                }

                if self.dry_run:
                    csv = buffer.getvalue()
                    idx = 0
                    for _ in range(2):
                        idx = csv.find("\n", idx) + 1
                    self.logger.info({"files": {"file": csv[:idx] + "..."}, **param})
                else:
                    r = requests.post(
                        "https://slack.com/api/files.upload",
                        params=param,
                        files={"file": buffer.getvalue()},
                    )
                    self.logger.info(r.text)

    def combine_transients(
        self, transients: Tuple[TransientView, ...]
    ) -> Tuple[List[pd.DataFrame], List[pd.DataFrame]]:
        """"""

        frames = []
        photometry = []

        for transient in transients:

            mycols = set(self.cols)

            frame = {
                "tranId": transient.id,
                "ztf_name": to_ztf_id(transient.id),
            }

            if (summary := transient.get_t2_result(unit_id="T2LightCurveSummary")) :
                frame.update(summary)

            # include other T2 results, flattened
            for t2record in ampel_iter(transient.t2):
                if (
                    t2record["unit"] == "T2LightCurveSummary"
                    or not t2record.get("body")
                    or not (output := t2record["body"][-1].get("result"))
                ):
                    continue

                # Flatten T2 output dictionary
                # If this is not a dictionry it will be through
                res_flat = flat_dict(output, prefix="T2-")

                # Add these to the dataframe (could we just join the dictionaries?)
                for key, value in res_flat.items():
                    try:
                        frame[key] = str(value)
                    except ValueError as ve:
                        log_exception(self.logger, ve)

            assert transient.stock
            frame.update(
                {channel: True for channel in ampel_iter(transient.stock["channel"])}
            )
            self.channels.update((str(c) for c in transient.stock["channel"] or []))

            frames.append(frame)

            if self.full_photometry:
                if transient.t0 is None:
                    raise ValueError(
                        "Full photometry requested, but no T0 records loaded"
                    )
                photometry.append(
                    pd.DataFrame(
                        [
                            {
                                "ztf_name": frame["ztf_name"],
                                **{
                                    k: v
                                    for pp in transient.t0
                                    if (
                                        not pp.get("tag")
                                        or not "SUPERSEDED" in pp["tag"]
                                    )
                                    for k, v in pp["body"].items()
                                    if k in mycols
                                },
                            }
                        ]
                    )
                )

        return frames, photometry


def calculate_excitement(n_transients, date, thresholds, n_alerts=np.nan):
    """"""

    message = "UPDATE! Alert summary for " + date + ". "

    if n_alerts == 0:
        message += "\n DISASTER! No alerts were ingested last night. :sob: "
        message += "Perhaps ZTF had technical problems? :telescope: \n"
        message += (
            "It might be worth checking the ZTF-General Slack channel "
            "for detector issues."
        )

        return message

    else:
        if not np.isnan(n_alerts):
            message += " In total, " + str(n_alerts) + " alerts were ingested. \n"
        else:
            message += "\n"

        if n_transients == 0:
            message += (
                "TRAGEDY! Sadly, it seems that no transients passed the "
                "filters last night"
            )
            return message

        elif n_transients == 1:
            message += (
                "LONE WOLF! It seems that there was just one transient"
                " which passed the filters. Guess it's better than "
                "nothing... :shrug: \n ~The results are summarised "
                "below.~ That one result is shown below."
            )
            return message

        elif n_transients < thresholds["Low"]:
            message += (
                "MEH! We only found " + str(n_transients) + " "
                "transients last night. :unamused: That's "
                "disappointing. Hopefully we get more tomorrow!"
            )

        elif n_transients < thresholds["Mid"]:
            message += (
                "NOT BAD! We found " + str(n_transients) + " "
                "transients last night. :slightly_smiling_face: "
                "Let's keep up the good work!"
            )

        elif n_transients < thresholds["High"]:
            message += (
                "IMPRESSIVE! We found "
                + str(n_transients)
                + " transients last night. :grin: That's exciting! Good "
                "luck to anyone trying to check all of those by hand..."
            )

        else:
            message += (
                "WOW!!! We found "
                + str(n_transients)
                + " transients last night. :tada: That's loads! Now we "
                "just need to figure out what to do with all of them..."
            )

        message += "\n The results are summarised below. "

    return message


def flat_dict(d, prefix=""):
    """
    Loop through dictionary d
    Append any key, val pairs to the return list ret
    Add the prefix to any key param
    Recurse if encountered value is a nested dictionary.
    """

    if not isinstance(d, collections.Mapping):
        return d

    ret = {}

    for key, val in d.items():
        if isinstance(val, collections.Mapping):
            ret = {**ret, **flat_dict(val, prefix=prefix + str(key) + "_")}
        else:
            ret[prefix + str(key)] = val

    return ret
