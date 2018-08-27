#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/SlackSummaryPublisher.py
# License           : BSD-3-Clause
# Author            : robert stein
# Date              : 11.03.2018
# Last Modified Date: 04.08.2018
# Last Modified By  : vb

import pandas as pd
import numpy as np
import collections
import io, pickle, datetime, requests
from slackclient import SlackClient

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.pipeline.common.AmpelUtils import AmpelUtils
from ampel.pipeline.logging.LoggingUtils import LoggingUtils

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
        self.photometry = []


    def add(self, transients):
        """
        """
        summary, full = self.combine_transients(transients)
        self.frames += summary
        self.photometry += full


    def done(self):
        """
        """
        if len(self.frames) == 0 and self.run_config.get('quiet', False):
            return
   
        try:
            date = self.run_config["date"]
        except KeyError:
            date = str(datetime.date.today())

        sc = SlackClient(self.run_config["Slack_token"])

        m = calculate_excitement(len(self.frames), date=date,
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
        
        if len(self.frames) > 0:

            df = pd.concat(self.frames, sort=False)
            photometry = pd.concat(self.photometry, sort=False)

            filename = "Summary_%s.csv" % date

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

            if self.run_config["full_photometry"]:

                filename = "Photometry_%s.csv" % date

                buffer = io.StringIO(filename)
                photometry.to_csv(buffer)

                param = {
                    'token': self.run_config["Slack_token"],
                    'channels': self.run_config["Slack_channel"],
                    'title': 'Full Photometry: ' + date,
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
        """
        """

        frames = []
        photometry = []

        for transient in transients:

            mycols = list(self.run_config["mycols"]) + list(self.run_config["channel(s)"])

            if not transient.photopoints:
                continue

            tdf = pd.DataFrame(
                [x.content for x in transient.photopoints]
            )

            # compute ZTF name
            tdf['ztf_name'] = tdf['tranId'].apply(AmpelUtils.get_ztf_name)
            tdf["most_recent_detection"] = max(tdf["jd"])
            tdf["first_detection"] = min(tdf["jd"])
            tdf["n_detections"] = len(tdf["jd"])
            
            try:
                if transient.t2records is not None:
                    for j, t2record in enumerate(transient.t2records):
                        if not t2record.results:
                            continue
                        res = (t2record.results[-1])
                        if not "output" in res:
                            continue

                        # Flatten T2 output dictionary
                        # If this is not a dictionry it will be through
                        res_flat = flat_dict(res['output'],prefix='T2-')

                        # Add these to the dataframe (could we just join the dictionaries?)                        
                        for key, value in res_flat.items():
                            try:
                                tdf[key] = value
                                mycols.append(key)
                            except ValueError as ve:
                                self.logger.error(ve)
            except:
                pass

            for channel in self.run_config["channel(s)"]:
                if channel in transient.channel:
                    tdf[channel] = True
                else:
                    tdf[channel] = False

            dfcols = list(tdf.columns.values)
            missing = [x for x in mycols if x not in dfcols]

            for col in missing:
                tdf[col] = "MISSING"

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


def flat_dict(d, prefix = ''):
    '''
    Loop through dictionary d
    Append any key, val pairs to the return list ret
    Add the prefix to any key param
    Recurse if encountered value is a nested dictionary.
    '''

    
    if not isinstance(d,collections.Mapping):
        return d
    
    ret = {}

    for key, val in d.items():
        if isinstance(val, collections.Mapping):
            ret = {**ret, **flat_dict(val, prefix = prefix+str(key)+'_') }
        else:
            ret[ prefix+str(key) ] = val

    return ret
