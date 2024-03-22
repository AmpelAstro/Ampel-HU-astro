#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/SubmitTNS.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                1.03.2024
# Last Modified Date:  1.03.2024
# Last Modified By:    jnordin@physik.hu-berlin.de

import asyncio
import time
from collections.abc import Generator
from typing import Any

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.contrib.hu.t3.tns.tns_ampel_util import get_tns_t2remarks, ztfdps_to_tnsdict
from ampel.contrib.hu.t3.tns.TNSClient import TNSClient
from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.contrib.hu.util.TNSMirrorSearcher import TNSMirrorSearcher
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import StockId, T3Send, UBson
from ampel.view.TransientView import TransientView


class SubmitTNS(AbsPhotoT3Unit, TNSMirrorSearcher):
    """
    Submit candidates to TNS (unless already submitted).

    Note that it is assumed that all selected transients are to be submitted.
    """

    # AT report config
    base_at_dict: dict = {
        "reporting_group_id": "82",  # Should be ampel
        "discovery_data_source_id": "48",
        "reporter": "J. Nordin, V. Brinnel, J. van Santen (HU Berlin), A. Gal-Yam, O. Yaron (Weizmann) on behalf of ZTF",
        "at_type": "1",
    }
    baseremark: str = "See arXiv:1904.05922 for selection criteria."

    # Connect information
    tns_key: NamedSecret[dict]
    timeout: float = 120.0
    max_parallel_requests: int = 8
    maxdist: float = 2.0  # max squared dist, in arcsec.
    tns_doublecheck: bool = True  # Also do a TNS name search - is this needed?
    tns_submit: bool = False  # Also do a TNS name search - is this needed?

    def post_init(self) -> None:
        self.client = TNSClient(
            TNSToken(**self.tns_key.get()),
            self.timeout,
            self.max_parallel_requests,
            self.logger,
        )

    async def get_tns_names(self, ra, dec):
        names = []
        async for doc in self.client.search(
            ra=ra, dec=dec, radius=self.maxdist, units="arcsec"
        ):
            names.extend(doc["internal_names"].split(", "))
        return names

    def sendReports(self, reports: list[dict]) -> dict:
        """
        Based on a lists of reportlists, send to TNS.
        Return results for journal entries
        """
        MAX_LOOP = 25
        SLEEP = 2

        reportresult = {}
        for atreport in reports:
            # Submit a report
            for _ in range(MAX_LOOP):
                reportid = asyncio.run(self.client.sendReport(atreport))
                if reportid:
                    print("TNS report ID %s" % (reportid))
                    break
                time.sleep(SLEEP)
            else:
                self.logger.info("TNS Report sending failed")
                print("TNS bulk report failed")
                continue

            # Try to read reply
            for _ in range(MAX_LOOP):
                time.sleep(SLEEP)
                response = asyncio.run(self.client.reportReply(reportid))
                if isinstance(response, list):
                    break
            else:
                self.logger.info("TNS Report reading failed")
                continue

            print("responseZZ", response)

            # Check whether request was bad. In this case TNS looks to return a list with dicts
            # of failed objects which does not correspond to the order of input atdicts.
            # In any case, nothing in the submit is posted.
            # Hence only checking first element
            bad_request = None
            for key_atprop in ["ra", "decl", "discovery_datetime"]:
                if key_atprop in response[0]:
                    try:
                        bad_request = response[0][key_atprop]["value"]["5"]["message"]
                        break
                    except KeyError:
                        pass
            if bad_request is not None:
                self.logger.info("bad TNS request", extra=bad_request)
                continue

            # Parse reply for evaluation
            for k, v in atreport["at_report"].items():
                if "100" in response[k]:
                    self.logger.info(
                        "TNS Inserted with name %s" % (response[k]["100"]["objname"])
                    )
                    reportresult[v["internal_name"]] = [
                        "TNS inserted",
                        {"TNSName": response[k]["100"]["objname"]},
                    ]
                elif "101" in response[k]:
                    self.logger.info(
                        "Already existing with name %s"
                        % (response[k]["101"]["objname"])
                    )
                    reportresult[v["internal_name"]] = [
                        "TNS pre-existing",
                        {"TNSName": response[k]["101"]["objname"]},
                    ]

        return reportresult

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson | UnitResult:
        # Reports to be sent, indexed by the transient view IDs (so that we can check in the replies)
        atreports: dict[StockId, dict[str, Any]] = {}

        for tran_view in gen:
            # Base information
            atdict = ztfdps_to_tnsdict(tran_view.get_photopoints())
            if atdict is None:
                self.logger.debug("Not enough info for TNS submission")
                print("... not enough info for submission")
                continue
            atdict.update(self.base_at_dict)

            # from T2s
            catremarks = get_tns_t2remarks(tran_view)
            if catremarks is not None:
                atdict.update(catremarks)

            # Was this already submitted to TNS, and if so also based on ZTF?
            (tnsname, internal_names) = self.get_tns_name_internal(
                atdict["ra"]["value"], atdict["dec"]["value"]
            )
            print("found tnsname", tnsname, internal_names)
            if atdict["internal_name"] in internal_names:
                print("... found it submitted in mirror db, skip")
                continue

            # and yet another check, directly from TNS... unnecessary?
            if self.tns_doublecheck:
                tnsmatch = asyncio.run(
                    self.get_tns_names(
                        ra=atdict["ra"]["value"], dec=atdict["dec"]["value"]
                    )
                )
                print("direct tns check", tnsmatch)
                if atdict["internal_name"] in tnsmatch:
                    print("... found it submitted at TNS!")
                    continue

            # Collected necessary data, not already published - add to submission list
            print("final atdict", atdict)
            atreports[tran_view.id] = atdict

        if len(atreports) == 0:
            # Nothing to submit
            self.logger.info("Nothing to report.")
            return None

        # atreports is now a dict with tran_id as keys and atreport as keys
        # what we need is a list of dicts with form {'at_report':atreport }
        # where an atreport is a dictionary with increasing integer as keys and atreports as values
        atreportlist = [
            {"at_report": {i: report for i, report in enumerate(atreports.values())}}
        ]

        print("final list of atreports")
        print(atreportlist)
        # Submit the reports
        if not self.tns_submit:
            return None
        tnsreply = self.sendReports(atreportlist)
        print("tnsreply", tnsreply)

        # Store report reply into t3 collection.
        return tnsreply
