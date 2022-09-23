#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ElasticcClassPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                11.04.2022
# Last Modified Date:  24.09.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from itertools import islice
from typing import Iterable, TYPE_CHECKING
from collections.abc import Generator

from ampel.struct.StockAttributes import StockAttributes
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.enum.DocumentCode import DocumentCode
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit, T3Send
from ampel.secret.NamedSecret import NamedSecret
from ampel.view.TransientView import TransientView
from ampel.view.T2DocView import T2DocView
from ampel.view.T3Store import T3Store

from ampel.contrib.hu.t3.ElasticcTomClient import ElasticcTomClient

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord

def chunks(l: Iterable, n: int) -> Generator[list, None, None]:
    source = iter(l)
    while True:
        chunk = list(islice(source, n))
        if chunk:
            yield chunk
        if len(chunk) < n:
            break

class ElasticcClassPublisher(AbsT3ReviewUnit):
    """

    This unit is intended to submit classifications to the DESC TOM db during
    the ELAsTICC LSST alert simulation.

    This will have to proceed through the following stages for each stock:
    - Retrieve all of the states and associate these to elasticc alertId.
    - Check logs for which of these states and t2configs a report was already (successfully) sent.
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

    desc_user: NamedSecret[str]
    desc_password: NamedSecret[str]

    #: prepare report, but do not submit to TOM
    dry_run: bool = False
    #: submit reports in batches
    batch_size: int = 1000

    unit: str = "T2ElasticcReport"

    def post_init(self) -> None:
        self.tomclient = ElasticcTomClient(self.desc_user.get(), self.desc_password.get(), self.logger)

    def search_journal_elasticc(
        self, tran_view: TransientView
    ) -> dict[int,list]:
        """
        Look through the journal for mapping between alert ID, timestampe and state id.

        Assuming journal entries from this unit has the following layout
        extra = {
            "t1State": t1_link,
            "t2Config": t2_config,
            "descPutResponse": response,
            "descPutComplete": True,
            "descPutUnit": self.unit,
            }

        Returns dict:
        {state:
            list(t2config)      # List of t2config for which report is done
              },
         ...}

        """

        # Create a list of states for which the list of units
        def select_submitted(entry: "JournalRecord") -> bool:
            return bool(
                (entry.get("extra") is not None and ("descPutComplete" in entry["extra"])
                and (entry["extra"]["descPutComplete"]) )
                and (entry["extra"]["descPutUnit"]==self.unit)
                and entry["unit"] and entry["unit"] == self.__class__.__name__
            )

        done_t1states: dict[int,list] = {}

        # All entries which we found should correspond to correctly sumitted classifications
        if jentries := list(tran_view.get_journal_entries(tier=3, filter_func=select_submitted)):
            for entry in jentries:
                if entry['extra']['t1State'] not in done_t1states.keys():
                    done_t1states[entry['extra']['t1State']] = []
                done_t1states[entry['extra']['t1State']].append( entry['extra']['t2Config'] )

        # Next section would look through the journal and find the elasticc alert
        # data needed. Here we are doing some short version of it
        # Perform another journal search and check whether this unit was run
        # for this state
        state_map = {}
        def select_alerts(entry: "JournalRecord") -> bool:
            return bool(
                entry.get("alert") is not None and entry.get("link") is not None
            )
        if jentries := list(tran_view.get_journal_entries(tier=0, filter_func=select_alerts)):
            for entry in jentries:
                assert isinstance(link := entry.get("link"), int)
                state_map[link] = done_t1states.get(link, [])
        return state_map


    def _get_reports(self, gen: Generator[TransientView, T3Send, None]) -> Generator[tuple[TransientView,int,dict,dict],None,None]:

        for tran_view in gen:
            # Check journal for state/alertId combos and whether already
            # submitted (for this t2classifiers list).
            state_alert = self.search_journal_elasticc(tran_view)

            for t1_link, submitted_configs in state_alert.items():
#                if submitted:
                if t2views := tran_view.get_t2_views(unit=self.unit, link=t1_link, code=DocumentCode.OK):
                    for t2view in t2views:
                        # Check whether report for t1link/state & t2config done
                        if t2view.config in submitted_configs:
                            self.logger.debug('submitted', extra={'t1':t1_link, 'config':t2view.config})
                            continue
                        if not isinstance((body := t2view.get_payload()), dict):
                            continue
                        yield tran_view, t1_link, t2view.config, body["report"]

    def process(self, gen: Generator[TransientView, T3Send, None], t3s: T3Store) -> None:
        """

        """

        for chunk in chunks(self._get_reports(gen), self.batch_size):
            tran_views, t1_links, t2_configs, class_reports = zip(*chunk)

            if self.dry_run:
                continue

            # use the ElasticcTomClient
            desc_response = self.tomclient.tom_post(class_reports)

            # Check output:
            # if as expected store to journal that transfer is complete.
            # if not as expected, log what data is available and possible
            # a t3 document with this content??
            for tran_view, t1_link, t2_config, class_report in zip(tran_views, t1_links, t2_configs, class_reports):
                if desc_response['success']:
                    gen.send((
                        tran_view.id,
                        StockAttributes(
                            journal=JournalAttributes(
                                extra={
                                    "t1State": t1_link,
                                    "t2Config": t2_config,
                                    "descPutResponse": desc_response,
                                    "descPutComplete": True,
                                    "descPutUnit": self.unit,
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
                                    "descPutUnits": self.unit,
                                    },
                                    ),
                                    )
                                ))
                    self.logger.info('desc post failed', extra={
                        "descResponse":desc_response,
                        "descReport": class_report, })
