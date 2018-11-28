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
import io, pickle, datetime, requests
from slackclient import SlackClient
from slackclient.exceptions import SlackClientError
from pydantic import BaseModel
from typing import Dict, List, Union
from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils
from ampel.pipeline.common.AmpelUtils import AmpelUtils
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.pipeline.config.EncryptedConfig import EncryptedConfig

class SlackSummaryPublisher(AbsT3Unit):
    """
    """

    class RunConfig(BaseModel):
        class Config:
            ignore_extra = False
        dryRun: bool = False
        quiet: bool = False
        slackChannel: str
        slackToken: Union[str, EncryptedConfig]
        excitement: Dict[str, int] = {"Low": 50,"Mid": 200,"High": 400}
        fullPhotometry: bool = False
        cols: List[str] = [
            "ztf_name","ra","dec","magpsf","sgscore1","rb",
            "most_recent_detection","first_detection","n_detections",
            "distnr","distpsnr1","isdiffpos","_id"
            ]
        requireNoAGN: bool = False
        requireNoSDSStar: bool = False
        requireNEDz: bool = False


    def __init__(self, logger, base_config=None, run_config=None, global_info=None):
        """
        """
        self.logger = AmpelLogger.get_logger() if logger is None else logger
        self.run_config = run_config
        self.frames = []
        self.photometry = []
        self.channels = set()

    def add(self, transients):
        """
        """
        summary, full = self.combine_transients(transients)
        self.frames += summary
        self.photometry += full


    def done(self):
        """
        """
        if len(self.frames) == 0 and self.run_config.quiet:
            return
   
        date = str(datetime.date.today())

        sc = SlackClient(self.run_config.slackToken)

        m = calculate_excitement(len(self.frames), date=date,
            thresholds=self.run_config.excitement
        )
        
        if self.run_config.dryRun:
            self.logger.info(m)
        else:
            api = sc.api_call(
                "chat.postMessage",
                channel=self.run_config.slackChannel,
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
                'token': self.run_config.slackToken,
                'channels': self.run_config.slackChannel,
                'title': 'Summary: ' + date,
                "username": "AMPEL-live",
                "as_user": "false",
                "filename": filename

            }

            if self.run_config.dryRun:
                # log only first two lines
                csv = buffer.getvalue()
                idx = 0
                for _ in range(2):
                    idx = csv.find('\n',idx)+1
                self.logger.info({"files": {"file": csv[:idx]+'...'}, **param})
            else:
                r = requests.post(
                    "https://slack.com/api/files.upload",
                    params=param,
                    files={"file": buffer.getvalue()}
                )
                self.logger.info(r.text)

            if self.run_config.fullPhotometry:
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
                    'token': self.run_config.slackToken,
                    'channels': self.run_config.slackChannel,
                    'title': 'Full Photometry: ' + date,
                    "username": "AMPEL-live",
                    "as_user": "false",
                    "filename": filename
                }

                if self.run_config.dryRun:
                    csv = buffer.getvalue()
                    idx = 0
                    for _ in range(2):
                        idx = csv.find('\n',idx)+1
                    self.logger.info({"files": {"file": csv[:idx]+'...'}, **param})
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

            mycols = list(self.run_config.cols)

            if not transient.photopoints:
                continue

            tdf = pd.DataFrame(
                [x.content for x in transient.photopoints]
            )

            # use tranId from parent view to compute ZTF name
            tdf['tranId'] = transient.tran_id
            tdf['ztf_name'] = tdf['tranId'].apply(ZTFUtils.to_ztf_id)
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
                                tdf[key] = str(value)
                                mycols.append(key)
                            except ValueError as ve:
                                self.logger.error(ve)
            except:
                pass

            if self.run_config.requireNEDz:
                if not "T2-NEDz_z" in mycols:
                    continue

            for channel in AmpelUtils.iter(transient.channel):
                tdf[channel] = True
                self.channels.add(channel)

            mycols += list(self.channels)
            missing = set(mycols).difference(tdf.columns.values)

            # deduplicate mycols, preserving order
            mycols = list(dict.fromkeys([x for x in mycols if not x in missing]))
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
