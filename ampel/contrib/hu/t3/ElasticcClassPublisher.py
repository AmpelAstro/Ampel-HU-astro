#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ElasticcClassPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                11.04.2022
# Last Modified Date:  11.04.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from typing import Any, Sequence, Dict
from collections.abc import Generator

from ampel.struct.StockAttributes import StockAttributes
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit, T3Send
from ampel.secret.NamedSecret import NamedSecret
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store

from ampel.contrib.hu.t3.ElasticcTomClient import ElasticcTomClient


# Tying classification output classes with ELAsTICC taxonomy classes
parsnip_taxonomy = {
    # Correct Elasticc outputs but they changed?
    # https://github.com/plasticc/taxonomy/blob/main/taxonomy.ipynb ?
    'SLSN':  131,
    'SNII':  113,
    'SNIa':  111,
    'SNIbc': 112,
    'TDE':  132,
    # Added for compatibility with old parsnip models
    'KN': 121,
    'SLSN-I': 131,
    'SNIa-91bg': 115,
    'SNIax': 114,
    }

class ElasticcClassPublisher(AbsT3ReviewUnit):
    """

    This unit is intended to submit classifications to the DESC TOM db during
    the ELAsTICC LSST alert simulation.

    This will have to proceed through the following stages for each stock:
    - Retrieve all of the states and associate these to elasticc alertId.
    - Check logs for which of these states a report was already (successfully) sent.
    - For still unsubmitted states, check whether T2 results exists for all
      classification units listed in the config.
    - If so, parse the T2Documents for the classifications and shape as ClassificationDict.
    - Try to submit using the ElasticcTomClient.
    - If successful, save info so state is not resubmitted next time.
    - If not successful, save not in log and wait for more T2s to complete.

    If T2 units have different run-times, we could have different instances
    running which looks for different units (fast/slow).
    Do we need to have some failsafe if some classifier does not return anything?
    (like submitting the rest after some time?)

    Note: still not sure whether we need to track also the diaSourceId and
    issue timestamp.

    """

    broker_name: str = 'AMPEL'
    broker_version: str = 'v0.1'
    desc_user: NamedSecret[str]
    desc_password: NamedSecret[str]
    tomclient = None


    t2classifiers: Sequence[str]
    # If we have classifiers running in multiple configurations we will
    # need to also sort based on config. Solve this if needed.

    def post_init(self) -> None:
        self.tomclient = ElasticcTomClient(self.desc_user, self.desc_password, self.logger)

    def search_journal_elasticc(
        self, tran_view: TransientView
    ) -> Dict[int,Any]:
        """
        Look through the journal for mapping between alert ID, timestampe and state id.

        Assuming journal entries from this unit has the following layout
        extra = {
            "t1State": t1_link,
            "descPutResponse": response,
            "descPutComplete": True,
            "descPutUnits": self.t2classifiers,
            }

        Returns dict:
        {state:{
              'alertId':111,         # elasticc alertId
              'diaSourceId':         # diaSourceId - where is this found?
              'elasticcPublishTimestamp':111,   # timestamp of kafka server
              'brokerIngestTimestamp': 111,    # timestamp of t0 ingestion
              'submitted': bool      # was this state submitted for these classifiers
              },
         ...}

        """

        # Create a list of states for which the list of units
        def select(entry: "JournalRecord") -> bool:
            return bool(
                (entry["extra"] is not None and ("descPutComplete" in entry["extra"])
                and (entry["extra"]["descPutComplete"]) )
                and (entry["extra"]["descPutUnits"]==tuple(self.t2classifiers))
                and entry["unit"] and entry["unit"] == self.__class__.__name__
            )

        done_t1states = []

        # All entries which we found should correspond to correctly sumitted classifications
        if jentries := list(tran_view.get_journal_entries(tier=3, filter_func=select)):
            for entry in jentries:
                done_t1states.append( entry['extra']['t1State'])

        # Next section would look through the journal and find the elasticc alert
        # data needed. Here we are doing some short version of it
        # Perform another journal search and check whether this unit was run
        # for this state
        state_map = {}
        def select(entry: "JournalRecord") -> bool:
            return bool(
                entry["alert"] is not None and entry["link"] is not None
            )
        if jentries := list(tran_view.get_journal_entries(tier=0, filter_func=select)):
            for entry in jentries:
                state_map[entry['link']] = {
                    'alertId': entry['alert'],
                    'diaSourceId': 666,
                    'brokerIngestTimestamp': entry['ts'],
                    'elasticcPublishTimestamp': 666 }
                if entry['link'] in done_t1states:
                    state_map[entry['link']]['submitted'] = True
                else:
                    state_map[entry['link']]['submitted'] = False
        return state_map



    def _get_parsnip_classification(self, t2view: TransientView) -> Dict:
        """
        Parse a T2Document from parsnip and produce a ClassificationDict
        as expected by ELAsTICC.

        Returns a list of dicts with the following structure:
        {
            'classifierName': 'Parsnip',
            'classifierParams':str(config),
            'classId': parsnip_taxonomy[parsnip_type],
            'probability': probability,
        }

        # Notes:
        # 1. Assumes that all redshift weighting and selection was done in T2RunParsnip
        # 2.

        """

        # Have already verified this to have completed.
        t2body = t2view.body[-1]

        # Check if Parsnip ran but not produce classification
        if not 'classification' in t2body:
            return []

        return [ {
                'classifierName': 'Parsnip',
                'classifierParams':str(t2view.config),
                'classId': parsnip_taxonomy[parsnip_type],
                'probability': prob,
                } for parsnip_type, prob in t2body['classification'].items() ]


    def process(self, gen: Generator[TransientView, T3Send, None], t3s: T3Store) -> None:
        """

        """


        for k, tran_view in enumerate(gen, 1):

            # Check journal for state/alertId combos and whether already
            # submitted (for this t2classifiers list).
            state_alert = self.search_journal_elasticc(tran_view)

            for t1_link, alertinfo in state_alert.items():
                if alertinfo['submitted']:
                    self.logger.info('submitted', extra={'t1':t1_link})
                    continue
                classifications = {}
                # Check whether we can get classification for all requested
                # units for this state
                for t2unit in self.t2classifiers:
                    if t2views := tran_view.get_t2_views(unit=t2unit, link=t1_link):
                        t2view = next(t2views, None)
                        if t2view is None:
                            self.logger.info('No T2Doc found', extra={'unit':t2unit})
                            continue   # No T2 ticket found
                        # Only reason there could be multiple views here is if we
                        # are running different configs... if so this unit wont work
                        # and needs to be redesigned!
                        if next(t2views, None) is not None:
                            sys.exit('ElasticcClassPublisher cannot parse multiple configs.')
                        # Skip if unit running or waiting to be run
                        if t2view.code == -1 or t2view.code == -4:
                            continue
                        # If run failed, make a note and let report go
                        if t2view.code == -3 or t2view.code == -6:
                            gen.send((
                                tran_view.id,
                                StockAttributes(
                                    journal=JournalAttributes(
                                        extra={
                                            "t1State": t1_link,
                                            "failedElasticcUnit": t2unit
                                            },
                                        ),
                                        tag="ELASTICC_FAIL",
                                )
                            ))
                            self.logger.info('ELAsTICC FAIL',
                                            extra={'unit':t2unit, 't1State': t1_link})
                            classifications[t2unit]=[]
                        elif t2view.code == 0:
                            if t2unit=='T2RunParsnip':
                                classifications['parsnip'] = self._get_parsnip_classification(t2view)
                            # Add the others... Rapid, TDE?, ...

                # Did we find all?
                if len(classifications)!=len(self.t2classifiers):
                    self.logger.info('Not all classifiers done', extra={'t1':t1_link})
                    continue

                # Time to submit!
                # Convert classifications to elasticc format
                class_report = {
                    'alertId': alertinfo['alertId'],
                    'diaSourceId': alertinfo['diaSourceId'],
                    'elasticcPublishTimestamp': alertinfo['elasticcPublishTimestamp'],
                    'brokerIngestTimestamp': alertinfo['brokerIngestTimestamp'],
                    'brokerName': self.broker_name,
                    'brokerVersion': self.broker_version,
                    'classifications': []
                    }
                for cname, clist in classifications.items():
                    class_report['classifications'].extend(clist)

                self.logger.debug('', extra={'classification_report': class_report})

                # use the ElasticcTomClient
                desc_response = self.tomclient.tom_post(class_report)

                # Check output:
                # if as expected store to journal that transfer is complete.
                # if not as expected, log what data is available and possible
                # a t3 document with this content??
                if desc_response['success']:
                    gen.send((
                        tran_view.id,
                        StockAttributes(
                            journal=JournalAttributes(
                                extra={
                                    "t1State": t1_link,
                                    "descPutResponse": desc_response,
                                    "descPutComplete": True,
                                    "descPutUnits": tuple(self.t2classifiers),
                                    "descReport": class_report,
                                    },
                                    ),
                                    )
                                ))
                else:
                    gen.send((
                        tran_view.id,
                        StockAttributes(
                            journal=JournalAttributes(
                                extra={
                                    "t1State": t1_link,
                                    "descPutResponse": desc_response,
                                    "descPutComplete": False,
                                    "descPutUnits": tuple(self.t2classifiers),
                                    "descReport": class_report,
                                    },
                                    ),
                                    )
                                ))
                    self.logger.info('desc post failed', extra={
                        "descResponse":desc_response,
                        "descReport": class_report, })
