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
from datetime import datetime, timedelta
from astropy import time
import os
import requests
import datetime
import time as t
from ampel.pipeline.db.DBUtils import DBUtils
import pandas as pd
import numpy as np


class SlackSummaryPublisher(AbsT3Unit):
    """
    """

    version = 2.0

    def __init__(self, logger, base_config=None):
        """
        """
        self.logger = LoggingUtils.get_logger() if logger is None else logger

    def run(self, run_config, raw_transients=None):

        if raw_transients is not None:
            transients = combine_transients(raw_transients, run_config)
        else:
            return

        try:
            date = run_config["date"]
        except KeyError:
            date = str(datetime.date.today())

        sc = SlackClient(run_config["Slack_token"])

        m = calculate_excitement(len(transients), date=date,
                                 thresholds=run_config["excitement_levels"],
                                 )

        api = sc.api_call(
                    "chat.postMessage",
                    channel=run_config["Slack_channel"],
                    text=m,
                    username="AMPEL-live",
                    as_user=False
                )

        self.logger.info(api)

        tmp_dir = "/Users/avocado/Ampel_output/"

        filename = run_config["name"] + "_" + date + ".csv"

        path = tmp_dir + filename

        transients.to_csv(path, index=False)

        with open(path, 'rb') as f:
            param = {
                'token': run_config["Slack_token"],
                'channels': run_config["Slack_channel"],
                'title': 'Summary: ' + date,
                "username": "AMPEL-live",
                "as_user": "false"

            }
            r = requests.post(
                "https://slack.com/api/files.upload",
                params=param,
                files={'file': f}
            )
            self.logger.info(r.text)

        os.remove(path)


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


def combine_transients(transients, t3_run_config):
    mycols = list(t3_run_config["mycols"])

    frames = []
    all_transients = []

    for transient in transients:

        if transient not in all_transients:
            all_transients.append(transient)
            tdf = pd.DataFrame([x.content for x in transient.get_photopoints()])

            # compute ZTF name
            tdf['ztf_name'] = tdf['tranId'].apply(DBUtils.get_ztf_name)
            tdf["most_recent_detection"] = max(tdf["jd"])
            tdf["first_detection"] = min(tdf["jd"])
            tdf["n_detections"] = len(tdf["jd"])

            # remove stupid columns and save to table
            frames.append(tdf[mycols][:1])

    #         else:
    #             index = all_transients.index(transient)
    #             frames[index][channel] = True

    df = pd.concat(frames)
    return df