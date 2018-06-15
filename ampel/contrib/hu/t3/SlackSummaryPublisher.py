#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/SlackPublisher.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 11.03.2018
# Last Modified Date: 17.03.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>
from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.pipeline.logging.LoggingUtils import LoggingUtils
from slackclient import SlackClient
from datetime import datetime
import requests
import datetime
from ampel.pipeline.common.AmpelUtils import AmpelUtils
import pandas as pd
import numpy as np
import io


class SlackSummaryPublisher(AbsT3Unit):
    """
    """

    version = 2.2

    def __init__(self, logger, base_config=None, run_config=None, global_info=None):
        """
        """
        self.logger = LoggingUtils.get_logger() if logger is None else logger
        self.base_config = base_config
        self.run_config = run_config
        self.frames = []

    def add(self, transients):
        self.frames += self.combine_transients(transients)

    def run(self):

        df = pd.concat(self.frames, sort=False)

        try:
            date = self.run_config["date"]
        except KeyError:
            date = str(datetime.date.today())

        sc = SlackClient(self.run_config["Slack_token"])

        m = calculate_excitement(len(df), date=date,
                                 thresholds=self.run_config["excitement_levels"]
                                 )

        api = sc.api_call(
                    "chat.postMessage",
                    channel=self.run_config["Slack_channel"],
                    text=m,
                    username="AMPEL-live",
                    as_user=False
                )

        self.logger.info(api)

        filename = "Summary_" + date + ".csv"

        buffer = io.StringIO(filename)
        df.to_csv(buffer)

        param = {
            'token': self.run_config["Slack_token"],
            'channels': self.run_config["Slack_channel"],
            'title': 'Summary: ' + date,
            "username": "AMPEL-live",
            "as_user": "false",
            "filename": filename

        }

        r = requests.post(
            "https://slack.com/api/files.upload",
            params=param,
            files={"file": buffer.getvalue()}
        )
        self.logger.info(r.text)

    def combine_transients(self, transients):

        mycols = list(self.run_config["mycols"]) + self.run_config["channel(s)"]

        frames = []

        for transient in transients:
            if transient.photopoints is None:
                continue

            tdf = pd.DataFrame(
                [x.content for x in transient.photopoints])

            # compute ZTF name
            tdf['ztf_name'] = tdf['tranId'].apply(AmpelUtils.get_ztf_name)
            tdf["most_recent_detection"] = max(tdf["jd"])
            tdf["first_detection"] = min(tdf["jd"])
            tdf["n_detections"] = len(tdf["jd"])

            for channel in self.run_config["channel(s)"]:
                if channel in transient.channel:
                    tdf[channel] = True
                else:
                    tdf[channel] = False

            # remove stupid columns and save to table
            existing = set(tdf.keys())
            frames.append(tdf[[k for k in mycols if k in existing]][:1])

        return frames


def calculate_excitement(n_transients, date, thresholds, n_alerts=np.nan):
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
                       "filters last night, even though alerts were ingested.  " \
                       ":cry: \n Maybe we should check if the other SWGs " \
                       "found things?"
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
