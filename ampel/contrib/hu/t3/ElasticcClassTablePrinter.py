#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ElasticcClassTablePrinter.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                11.04.2022
# Last Modified Date:  14.09.2023
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

import os
from collections import defaultdict
from collections.abc import Generator
from typing import TYPE_CHECKING

import pandas as pd

from ampel.abstract.AbsT3Unit import AbsT3Unit, T3Send
from ampel.contrib.hu.t3.ElasticcTomClient import Elasticc2ClassificationDict
from ampel.enum.DocumentCode import DocumentCode
from ampel.log import LogFlag
from ampel.struct.T3Store import T3Store
from ampel.view.TransientView import TransientView

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord


class ElasticcClassTablePrinter(AbsT3Unit):
    """

    Unit which will print Elasticc unit classifications to STDOUT and/or an output csv file

    """

    # Map of classifier internal names to plot, to output names to use
    classifier_map: dict = {
        "ElasticcLiveSNGuess": "SNGuess",
        "ElasticcLiveSNGuessParsnip": "FollowMe",
        "ElasticcLiveSNGuessParsnipPrior": "Final",
    }
    # Only include probabilities above this limit
    min_prob: float = 0.05

    #: Output settings for saving pandas *csv* file
    print_stdout: bool = False
    output_dir: str = "."
    output_file: None | str = None  # Set to name to generate

    # Name of unit where reports were collected
    unit: str = "T2ElasticcReport"

    def search_journal_elasticc(self, tran_view: TransientView) -> dict[int, bool]:
        """
        Look through the journal for mapping between alert ID, timestampe and state id.

        Assuming journal entries from this unit has the following layout
        extra = {
            "t1State": t1_link,
            "descPutResponse": response,
            "descPutComplete": True,
            "descPutUnit": self.unit,
            }

        Returns dict:
        {state:
            bool      # report submitted for this state
              },
         ...}

        """

        # Create a list of states for which the list of units
        def select_submitted(entry: "JournalRecord") -> bool:
            return bool(
                (
                    entry.get("extra") is not None
                    and ("descPutComplete" in entry["extra"])
                    and (entry["extra"]["descPutComplete"])
                )
                and (entry["extra"]["descPutUnit"] == self.unit)
                and entry["unit"]
                and entry["unit"] == self.__class__.__name__
            )

        done_t1states = set()

        # All entries which we found should correspond to correctly sumitted classifications
        if jentries := list(
            tran_view.get_journal_entries(tier=3, filter_func=select_submitted)
        ):
            for entry in jentries:
                done_t1states.add(entry["extra"]["t1State"])

        # Next section would look through the journal and find the elasticc alert
        # data needed. Here we are doing some short version of it
        # Perform another journal search and check whether this unit was run
        # for this state
        state_map = {}

        def select_alerts(entry: "JournalRecord") -> bool:
            return bool(
                entry.get("alert") is not None and entry.get("link") is not None
            )

        if jentries := list(
            tran_view.get_journal_entries(tier=0, filter_func=select_alerts)
        ):
            for entry in jentries:
                link = entry.get("link")
                assert isinstance(link, int)
                if link in done_t1states:
                    state_map[link] = True
                else:
                    state_map[link] = False
        return state_map

    def _parse_classifier_reports(self, t2ElasticcReport: dict):
        """
        Reconstruct T2ElasticcReport output
        """
        flatreport = {
            "alertId": t2ElasticcReport["alertId"],
            "diaSourceId": t2ElasticcReport["diaSourceId"],
        }
        # Sort classifications
        by_classifier: defaultdict[
            tuple[str, str], list[Elasticc2ClassificationDict]
        ] = defaultdict(list)

        # Elasticc1 allowed multiple classifiers in each report, storing names with the probabilities, while e2 assumed a single named stored in the report body
        if (cname := t2ElasticcReport.get("classifierName")) is not None:
            if cname not in self.classifier_map:
                # elasticc 2, and this classifier not requested
                return []
            by_classifier[cname, t2ElasticcReport["classifierParams"]] = []
            for c in t2ElasticcReport["classifications"]:
                if c["probability"] < self.min_prob:
                    continue
                by_classifier[cname, t2ElasticcReport["classifierParams"]].append(
                    {"classId": c["classId"], "probability": c["probability"]}
                )
        else:
            for c in t2ElasticcReport["classifications"]:
                if c["probability"] < self.min_prob:
                    continue
                if c["classifierName"] not in self.classifier_map:
                    continue
                by_classifier[(c["classifierName"], c["classifierParams"])].append(
                    {"classId": c["classId"], "probability": c["probability"]}
                )

        # Parse classifications
        for (name, _), classifications in by_classifier.items():
            flatreport.update(
                {
                    "{}_P({})".format(self.classifier_map[name], c["classId"]): c[
                        "probability"
                    ]
                    for c in classifications
                }
            )
        return flatreport

    def _get_reports(
        self, gen: Generator[TransientView, T3Send, None]
    ) -> Generator[tuple[TransientView, int, dict], None, None]:
        stats = {
            "stocks": 0,
            "views": 0,
            "pending": 0,
            "submitted": 0,
        }

        for tran_view in gen:
            state_alert = self.search_journal_elasticc(tran_view)

            stats["stocks"] += 1

            for t1_link, submitted in state_alert.items():
                stats["views"] += 1
                if submitted:
                    self.logger.debug("submitted", extra={"t1": t1_link})
                    stats["submitted"] += 1
                    continue
                if t2views := tran_view.get_t2_views(
                    unit=self.unit, link=t1_link, code=DocumentCode.OK
                ):
                    t2view = next(t2views, None)
                    if t2view is None:
                        self.logger.debug("No T2Doc found", extra={"unit": self.unit})
                        stats["pending"] += 1
                        continue  # No T2 ticket found
                    # Only reason there could be multiple views here is if we
                    # are running different configs... if so this unit wont work
                    # and needs to be redesigned!
                    if (
                        t2_view_extra := next(t2views, None)
                    ) is not None and t2_view_extra.config != t2view.config:
                        raise RuntimeError(
                            "ElasticcClassPublisher cannot parse multiple configs. "
                            f"Got configs={t2view.config},{t2_view_extra.config} for "
                            f"stock:{t2view.stock},link:{t2view.link},unit:{t2view.unit}"  # type: ignore[str-bytes-safe]
                        )
                    if not isinstance((body := t2view.get_payload()), dict):
                        continue
                    # Either a singular report or a list of reports
                    if (report := body.get("report", None)) is not None:
                        yield tran_view, t1_link, report
                    elif isinstance((reports := body["reports"]), tuple):
                        for report in reports:
                            yield tran_view, t1_link, report
                    else:
                        yield tran_view, t1_link, body["report"]
        self.logger.log(LogFlag.SHOUT, "filtered states", extra=stats)

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: T3Store
    ) -> None:
        """ """

        reports = []
        for tran_view, t1_link, class_report in self._get_reports(gen):
            flatreport = self._parse_classifier_reports(class_report)
            flatreport["stock"] = tran_view.id
            flatreport["t1_link"] = t1_link

            if self.print_stdout:
                print(flatreport)  # noqa: T201
            reports.append(flatreport)

        # Save dataframe
        if self.output_file is not None:
            df = pd.DataFrame.from_records(reports)
            full_path = os.path.join(self.output_dir, self.output_file)
            df.to_csv(full_path)
