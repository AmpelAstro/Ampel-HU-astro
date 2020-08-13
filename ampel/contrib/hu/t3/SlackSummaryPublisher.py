#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/SlackSummaryPublisher.py
# License           : BSD-3-Clause
# Author            : robert stein
# Date              : 11.03.2018
# Last Modified Date: 14.11.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import pandas as pd
import numpy as np
import collections
import io, datetime, requests
from slack import WebClient
from slack.errors import SlackClientError
from typing import Dict, List, Union, Optional, Any
from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.ztf.utils import to_ampel_id, to_ztf_id
from ampel.util.collections import ampel_iter
from ampel.model.EncryptedDataModel import EncryptedDataModel


class SlackSummaryPublisher(AbsT3Unit):
    """
    """

    dry_run: bool = False
    quiet: bool = False
    slack_channel: str
    slack_token: Union[str, EncryptedDataModel]
    excitement: Dict[str, int] = {"Low": 50, "Mid": 200, "High": 400}
    full_photometry: bool = False
    cols: List[str] = [
        "ztf_name", "ra", "dec", "magpsf", "sgscore1", "rb",
        "last_significant_nondet", "first_detection",
        "most_recent_detection", "n_detections",
        "distnr", "distpsnr1", "isdiffpos", "_id"
    ]
    requireNoAGN: bool = False
    requireNoSDSStar: bool = False
    requireNEDz: bool = False


    def post_init(self) -> None:
        """ """
        self.frames = []
        self.photometry = []
        self.channels = set()
        assert self.dry_run


    def add(self, transients):
        """
        """
        summary, full = self.combine_transients(transients)
        self.frames += summary
        self.photometry += full


    def done(self):
        """
        """
        if len(self.frames) == 0 and self.quiet:
            return

        date = str(datetime.date.today())

        sc = WebClient(self.slack_token)

        m = calculate_excitement(len(self.frames), date=date,
            thresholds=self.excitement
        )

        if self.dry_run:
            print(m, self.logger)
            self.logger.info(m)
        else:
            api = sc.api_call(
                "chat.postMessage",
                channel=self.slack_channel,
                text=m,
                username="AMPEL-live",
                as_user=False
            )
            if not api['ok']:
                raise SlackClientError(api['error'])

        if len(self.frames) > 0:

            df = pd.concat(self.frames, sort=False)
            # Set fill value for channel columns to False
            for channel in self.channels:
                df[channel].fillna(False, inplace=True)
            # Set fill value for all other columns to MISSING
            for field in set(df.columns.values).difference(self.channels):
                df[field].fillna('MISSING', inplace=True)

            filename = "Summary_%s.csv" % date

            buffer = io.StringIO(filename)
            df.to_csv(buffer)

            param = {
                'token': self.slack_token,
                'channels': self.slack_channel,
                'title': 'Summary: ' + date,
                "username": "AMPEL-live",
                "as_user": "false",
                "filename": filename

            }

            if self.dry_run:
                # log only first two lines
                csv = buffer.getvalue()
                idx = 0
                for _ in range(2):
                    idx = csv.find('\n', idx) + 1
                self.logger.info({"files": {"file": csv[:idx] + '...'}, **param})
            else:
                r = requests.post(
                    "https://slack.com/api/files.upload",
                    params=param,
                    files={"file": buffer.getvalue()}
                )
                self.logger.info(r.text)

            if self.full_photometry:
                photometry = pd.concat(self.photometry, sort=False)
                # Set fill value for channel columns to False
                for channel in self.channels:
                    photometry[channel].fillna(False, inplace=True)
                # Set fill value for all other columns to MISSING
                for field in set(photometry.columns.values).difference(self.channels):
                    photometry[field].fillna('MISSING', inplace=True)

                filename = "Photometry_%s.csv" % date

                buffer = io.StringIO(filename)
                photometry.to_csv(buffer)

                param = {
                    'token': self.slack_token,
                    'channels': self.slack_channel,
                    'title': 'Full Photometry: ' + date,
                    "username": "AMPEL-live",
                    "as_user": "false",
                    "filename": filename
                }

                if self.dry_run:
                    csv = buffer.getvalue()
                    idx = 0
                    for _ in range(2):
                        idx = csv.find('\n', idx) + 1
                    self.logger.info({"files": {"file": csv[:idx] + '...'}, **param})
                else:
                    r = requests.post(
                        "https://slack.com/api/files.upload",
                        params=param,
                        files={"file": buffer.getvalue()}
                    )
                    self.logger.info(r.text)


    def combine_transients(self, transients):
        """
        """

        frames = []
        photometry = []

        for transient in transients:

            mycols = list(self.cols)

            if not transient.t0:
                continue

            tdf = pd.DataFrame(
                [x["body"] for x in transient.t0]
            )

            # NB: photopoint fields like jd, diffmaglim, fid are ZTF (IPAC) specific
            # use tranId from parent view to compute ZTF name
            tdf['tranId'] = transient.id
            tdf['ztf_name'] = tdf['tranId'].apply(to_ztf_id)
            tdf["most_recent_detection"] = max(tdf["jd"])
            tdf["first_detection"] = min(tdf["jd"])
            tdf["n_detections"] = len(tdf["jd"])

            # Parse upper limits if present for the last upper limit prior to detection
            # As is, an upper limit between two detections (so after the first) will not reset this clock
            # It can be discussed whether this is the requested behaviour (see jd>min(tdf["first_detection"]) ) below
            # For example case, look at ZTF19abejaiy
            # Only include "significant" (limit deeper than 19.5)
            upper_limits = [pp for pp in transient.t0 if pp['_id'] < 0]
            if upper_limits is not None:
                jd_last_nondet = 0
                mag_last_nondet = 99
                filter_last_nondet = 99
                for ulim in upper_limits:
                    jd = ulim["body"]["jd"]
                    if jd < jd_last_nondet or jd > min(tdf["first_detection"]):
                        continue
                    ul = ulim["body"]["diffmaglim"]
                    if ul < 19.5:
                        continue
                    jd_last_nondet = jd
                    mag_last_nondet = ul
                    filter_last_nondet = ulim["body"]["fid"]
                if jd_last_nondet > 0:
                    tdf["last_significant_nondet_jd"] = jd_last_nondet
                    tdf["last_significant_nondet_mag"] = mag_last_nondet
                    tdf["last_significant_nondet_fid"] = filter_last_nondet

            if transient.t2:
                for j, t2record in enumerate(transient.t2):
                    if not (output := t2record.get("body", [{}])[-1].get("result")):
                        continue

                    # Flatten T2 output dictionary
                    # If this is not a dictionry it will be through
                    res_flat = flat_dict(output, prefix='T2-')

                    # Add these to the dataframe (could we just join the dictionaries?)
                    for key, value in res_flat.items():
                        try:
                            tdf[key] = str(value)
                            mycols.append(key)
                        except ValueError as ve:
                            self.logger.error(ve)

            if self.requireNEDz:
                if "T2-NEDz_z" not in mycols:
                    continue

            for channel in ampel_iter(transient.stock["channel"]):
                tdf[channel] = True
                self.channels.add(channel)

            mycols += list(self.channels)
            missing = set(mycols).difference(tdf.columns.values)

            # deduplicate mycols, preserving order
            mycols = list(dict.fromkeys([x for x in mycols if x not in missing]))
            # remove stupid columns and save to table
            frames.append(tdf[mycols][:1])
            photometry.append(tdf)
            # frames.append(tdf[:1])

        return frames, photometry


def calculate_excitement(n_transients, date, thresholds, n_alerts=np.nan):
    """
    """

    message = "UPDATE! Alert summary for " + date + ". "

    if n_alerts == 0:
        message += "\n DISASTER! No alerts were ingested last night. :sob: "
        message += "Perhaps ZTF had technical problems? :telescope: \n"
        message += "It might be worth checking the ZTF-General Slack channel " \
                   "for detector issues."

        return message

    else:
        if not np.isnan(n_alerts):
            message += " In total, " + str(n_alerts) + \
                       " alerts were ingested. \n"
        else:
            message += "\n"

        if n_transients == 0:
            message += "TRAGEDY! Sadly, it seems that no transients passed the " \
                       "filters last night"
            return message

        elif n_transients == 1:
            message += "LONE WOLF! It seems that there was just one transient" \
                       " which passed the filters. Guess it's better than " \
                       "nothing... :shrug: \n ~The results are summarised " \
                       "below.~ That one result is shown below."
            return message

        elif n_transients < thresholds["Low"]:
            message += "MEH! We only found " + str(n_transients) + " " \
                           "transients last night. :unamused: That's " \
                           "disappointing. Hopefully we get more tomorrow!"

        elif n_transients < thresholds["Mid"]:
            message += "NOT BAD! We found " + str(n_transients) + " " \
                      "transients last night. :slightly_smiling_face: " \
                      "Let's keep up the good work!"

        elif n_transients < thresholds["High"]:
            message += "IMPRESSIVE! We found " + str(n_transients) + \
                       " transients last night. :grin: That's exciting! Good " \
                       "luck to anyone trying to check all of those by hand..."

        else:
            message += "WOW!!! We found " + str(n_transients) + \
                       " transients last night. :tada: That's loads! Now we " \
                       "just need to figure out what to do with all of them..."

        message += "\n The results are summarised below. "

    return message


def flat_dict(d, prefix = ''):
    '''
    Loop through dictionary d
    Append any key, val pairs to the return list ret
    Add the prefix to any key param
    Recurse if encountered value is a nested dictionary.
    '''

    if not isinstance(d, collections.Mapping):
        return d

    ret = {}

    for key, val in d.items():
        if isinstance(val, collections.Mapping):
            ret = {**ret, **flat_dict(val, prefix = prefix + str(key) + '_')}
        else:
            ret[prefix + str(key)] = val

    return ret
