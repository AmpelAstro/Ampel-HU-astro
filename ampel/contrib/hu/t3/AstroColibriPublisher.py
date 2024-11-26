#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/AstroColibriPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                12.11.2022
# Last Modified Date:  13.11.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

# FIXME: restore mypy when this is actually ready
# type: ignore

# from itertools import islice
import os
from collections.abc import Generator

import numpy as np
from astropy.time import Time

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.content.JournalRecord import JournalRecord
from ampel.contrib.hu.t3.AstroColibriClient import (
    AstroColibriClient,
    AstroColibriPhot,
)
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.struct.T3Store import T3Store
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id


def dps_to_astrocolobri(dps_det: list[dict]) -> dict[str, AstroColibriPhot]:
    """
    From a list of datapoints, construct a dict with the photometry
    information expected by AstroColibri
    {
        'first': {'mag','band','datetime'},
        'peak': {'mag','band','datetime'},
        'last': {'mag','band','datetime'},
    }

    """

    def convert_dp(dp):
        bandmap = {1: "ZTFg", 2: "ZTFR", 3: "ZTFi"}
        return {
            "mag": dp["body"]["magpsf"],
            "band": bandmap[dp["body"]["fid"]],
            "datetime": Time(dp["body"]["jd"], format="jd").iso,
        }

    outphot = {}

    # Get the first and last
    dps_sort = sorted(dps_det, key=lambda pp: pp["body"]["jd"])
    outphot["first"] = convert_dp(dps_sort[0])
    outphot["last"] = convert_dp(dps_sort[-1])

    # Get the peak
    dps_sort = sorted(dps_det, key=lambda pp: pp["body"]["magpsf"])
    outphot["peak"] = convert_dp(dps_sort[0])

    return outphot


class AstroColibriPublisher(AbsPhotoT3Unit):
    """

    Publish results to AstroColibri. This demo version will:
    - Find the first, brightest and last photometry.
    - Get the position.
    Collect attributes:
    - "Nearby" if AmpelZ<0.02
    - "ProbSNIa" if ParsnipP(SNIa)>0.5
    - "ProbSN" if SNGuess=True at any phase.
    - Kilonovaness if available.

    Will update if new obs was made after last posting.

    AC requests:


    - Can you make sure that the event time corresponds to the one submitted to TNS?
    - We suggest a slight renaming of the events:
     x   * source_name: use the name of the event given by TNS and add "(Ampel)" to it. Example: TNS name "AT 2024edy" => "AT 2024edy (Ampel)"
     x   * trigger_id: keep as it is
     x   * we suggest you fill the ZTF identifier (e.g. ZTF24aahbwis) that you currently use a source_name into the field "discoverer_internal_name"
    x Type: for events that are not submitted as classified to TNS (i.e. listed as AT in TNS), please change the value of the "type" parameter to "ot" (instead of "ot_sn")
    * Classification: you can add the Ampel classification into the field "classification". E.g. if these are supernova candidates use "SN" or any classification that one could also find on TNS






    """

    # Limits for attributes
    nearby_z: float = 0.02
    snia_minprob: float = 0.7  # Parsnip class prob to call something SNIa
    min_kilonovaness: float = 0.0  # Parsnip class prob to call something SNIa

    # Image upload
    # ZTFNAME will be replaced with the ZTF interpreted stock id
    image_path: None | str = None

    # TODO: What is the use here, when selecting matching journal entries.
    process_name: None | str = None

    # Testcollect (will not try to post)
    testcollect: bool = False
    # randname (for repeat submission to dev AstroColibri, not checking TNS)
    randname: bool = False

    #    user: NamedSecret[str]
    #    password: NamedSecret[str]
    user: str = "AMPEL"
    password: str = "AMPEL_in_Astro-COLIBRI"

    def post_init(self) -> None:
        #        self.colibriclient = AstroColibriClient(self.user.get(), self.password.get(), self.logger)
        self.colibriclient = AstroColibriClient(self.user, self.password, self.logger)

    def _filter_journal_astrocolobri(self, jentry: JournalRecord, after: float):
        """
        Select journal entries from AstroColibriPublisher newer than last update
        Does not work when passed as a function to get_journal_entries...
        """
        return (
            jentry["unit"] == "AstroColibriPublisher"
            and (self.process_name is None or jentry["process"] == self.process_name)
            and jentry["ts"] >= after
            and jentry.get("extra") is not None
            and jentry["extra"]["success"]
            and not jentry["extra"].get("testcollect", False)
            and not jentry["extra"].get("randname", True)  # Allow submission of same
            # Possibly add check for whether submit was ok: sucess=True. But see how this looks in journal.
            # and (jentry.get("extra") is not None and ("descPutComplete" in entry["extra"]) )
            # and (entry["extra"]["descPutComplete"]) )
        )

    def requires_update(self, view: TransientView) -> bool:
        if not view.stock:
            return False
        # find latest activity activity at lower tiers
        latest_activity = max(
            (
                jentry["ts"]
                for jentry in view.get_journal_entries() or []
                if jentry.get("tier") in {0, 1, 2}
            ),
            default=float("inf"),
        )
        # Manual
        t3journals = [
            je
            for je in view.get_journal_entries(tier=3)
            if self._filter_journal_astrocolobri(je, latest_activity)
        ]
        return bool(t3journals)

    def submitted(self, view: "TransientView") -> bool:
        # Was transient (successfully pushed)
        if not view.stock:
            return False
        return bool(
            [
                je
                for je in view.get_journal_entries(tier=3)
                if self._filter_journal_astrocolobri(je, 0)
            ]
        )

    def process(
        self,
        tviews: Generator[TransientView, JournalAttributes, None],
        t3s: "None | T3Store" = None,
    ) -> None:
        """
        Iterate through TransientView and check which should be
        submitted / updated.
        """

        for tview in tviews:
            # Gather general information, including coordinate
            payload = {
                "type": "ot",
                "observatory": "ztf",
                # "source_name": to_ztf_id(int(tview.id)),
                "discoverer_internal_name": to_ztf_id(int(tview.id)),
                #                'trigger_id': self.trigger_id+':'+str(tview.id),  # How do these work?
                # "trigger_id": "TNS" + tns_name,
                #                'ivorn': self.trigger_id+':'+str(tview.id),       # Need ivorn schema
                # "timestamp": Time.now().iso,
            }

            # Gather photometry based information
            assert tview.t0 is not None
            # Get subset of real detections.
            # Could be complemented with further restrictions, e.g. filter / RB
            dps_det = [
                pp
                for pp in tview.t0
                if (pp["id"] > 0 and pp["body"].get("isdiffpos", False))
            ]
            payload["photometry"] = dps_to_astrocolobri(dps_det)
            payload["ra"] = np.mean([pp["body"]["ra"] for pp in dps_det])
            payload["dec"] = np.mean([pp["body"]["dec"] for pp in dps_det])
            payload["err"] = 1.0 / 3600  # position err ~1 arcsec in dec

            # Find TNS name
            if self.randname:
                import random

                # Generate random name (assuming publishing to dev AC)
                tns_name = f"AmpelRand{random.randint(1, 999)}"
                tns_submission_time = Time.now().iso
            elif isinstance(tview.extra, dict) and "TNSReports" in tview.extra:
                # A tns name is required, here obtained from the mirror DB through a T3 complement
                tnsreport = next(iter(tview.extra["TNSReports"]))
                tns_name = tnsreport["objname"]
                tns_submission_time = tnsreport["discoverydate"]
            else:
                self.logger.debug("No TNS name", extra={"tnsName": None})
                continue

            payload["trigger_id"] = "TNS" + tns_name
            payload["source_name"] = tns_name + " (AMPEL)"
            payload["time"] = tns_submission_time

            # Check if this was submitted
            # TODO: How should the first submit differ from updates?
            if self.submitted(tview):
                # Check if it needs an update
                if self.requires_update(tview):
                    post_update = True
                else:
                    continue
            else:
                post_update = False

            # If part of random testing, perturb coordinates
            if self.randname:
                payload["ra"] += random.randrange(-1, 1)
                payload["dec"] += random.randrange(-1, 1)
                payload["source_name"] += payload["trigger_id"][12:]

            # Gather attributes
            attributes = {"classification": {}, "ampelProp": []}
            # Nearby attribute
            t2res = tview.get_t2_body(unit="T2DigestRedshifts")
            if isinstance(t2res, dict) and t2res.get("ampel_z", 999) < self.nearby_z:
                attributes["ampelProp"].append("Nearby")
                attributes["ampelProp"].append("AmpelZ{:.2f}".format(t2res["ampel_z"]))
            # Infant attribute
            t2res = tview.get_t2_body(unit="T2InfantCatalogEval")
            if isinstance(t2res, dict) and t2res.get("action", False):
                attributes["ampelProp"].append("Young")
            # SNIa
            t2res = tview.get_t2_body(unit="T2RunParsnip")
            if (
                isinstance(t2res, dict)
                and "classification" in t2res
                and t2res["classification"]["SNIa"] > self.snia_minprob
            ):
                attributes["ampelProp"].append("ProbSNIa")
                attributes["classification"] = {
                    "class": ["SNIa", "Other"],
                    "prob": [
                        t2res["classification"]["SNIa"],
                        1 - t2res["classification"]["SNIa"],
                    ],
                }
            # Kilonovaness
            t2res = tview.get_t2_body(unit="T2KilonovaEval")
            if (
                isinstance(t2res, dict)
                and t2res.get("kilonovaness", -99) > self.min_kilonovaness
            ):
                attributes["ampelProp"].append(
                    "Kilonovaness{}".format(t2res["kilonovaness"])
                )
                attributes["ampelProp"].append("LVKmap{}".format(t2res["map_name"]))

            # Check whether we have a figure to upload.
            # Assuming this exists locally under {stock}.png
            if self.image_path is not None:
                ipath = self.image_path.replace("ZTFNAME", to_ztf_id(int(tview.id)))
                ipath = self.image_path.replace("STOCK", str(tview.id))
                # Only upload if it actually exists:
                if not os.path.isfile(ipath):
                    ipath = None
            else:
                ipath = None

            payload["broker_attributes"] = attributes
            self.logger.debug("reacting", extra={"payload": payload})

            # Ok, so we have a transient to react to
            if self.testcollect:
                if post_update:
                    self.logger.debug("Should send an update")
                else:
                    self.logger.debug("New submission")
                self.logger.debug(payload)

                # Journal submission
                jcontent = {
                    "reaction": "fake submission",
                    "success": True,
                    "testcollect": True,
                }

            else:
                self.logger.debug(
                    "colibri submitting",
                    extra={"payload": payload, "image_path": ipath},
                )
                jcontent = self.colibriclient.firestore_post(payload, image_path=ipath)

            if jcontent:
                tviews.send(JournalAttributes(extra=jcontent))
