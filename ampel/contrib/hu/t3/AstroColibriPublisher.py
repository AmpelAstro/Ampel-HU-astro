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
import re
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

    def get_tnsname(self, view: "TransientView") -> bool:
        # Was transient (successfully pushed)
        if not view.stock:
            return False
        # Check whether the name is found in the name collection
        if len(names := view.stock.get("name", [])) > 0:
            # Will only be able to require TNS name through format
            # dddd
            for name in names:
                if re.search(r"\d{4}\D{3}\D?", name):
                    return name

        # Should we look through the Journal for entries from the TNSTalker?
        # It *should* also save these entries to name so should not be needed...

        return None

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
            if self.randname:
                import random

                # Generate random name (assuming publishing to dev AC)
                tns_name = f"AmpelRand{random.randint(1, 999)}"
            else:
                # Find TNS name (required for AstroColibri posting)
                # Currently assumes that this is stored either in the
                # stock name list or can be found in the T3 journal
                # (probably from the TNSTalker)
                tns_name = self.get_tnsname(tview)
                if not tns_name:
                    self.logger.info("No TNS.", extra={"tnsName": None})
                    continue

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

            # Gather general information
            payload = {
                "type": "ot_sn",  # for optical? when to change to ot_sn?
                "observatory": "ztf",
                "source_name": to_ztf_id(int(tview.id)),
                #                'trigger_id': self.trigger_id+':'+str(tview.id),  # How do these work?
                "trigger_id": "TNS" + tns_name,
                #                'ivorn': self.trigger_id+':'+str(tview.id),       # Need ivorn schema
                "timestamp": Time.now().iso,
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

            # If part of random testing, perturb coordinates
            if self.randname:
                payload["ra"] += random.randrange(-1, 1)
                payload["dec"] += random.randrange(-1, 1)
                payload["source_name"] += payload["trigger_id"][12:]

            # Gather attributes
            attributes = []
            # Nearby attribute
            t2res = tview.get_t2_body(unit="T2DigestRedshifts")
            if isinstance(t2res, dict) and t2res.get("ampel_z", 999) < self.nearby_z:
                attributes.append("Nearby")
                attributes.append("AmpelZ{:.2f}".format(t2res["ampel_z"]))
            # Infant attribute
            t2res = tview.get_t2_body(unit="T2InfantCatalogEval")
            if isinstance(t2res, dict) and t2res.get("action", False):
                attributes.append("Young")
            # SNIa
            t2res = tview.get_t2_body(unit="T2RunParsnip")
            if isinstance(t2res, dict) and "classification" in t2res.keys():
                if t2res["classification"]["SNIa"] > self.snia_minprob:
                    attributes.append("ProbSNIa")
            # Kilonovaness
            t2res = tview.get_t2_body(unit="T2KilonovaEval")
            if (
                isinstance(t2res, dict)
                and t2res.get("kilonovaness", -99) > self.min_kilonovaness
            ):
                attributes.append("Kilonovaness{}".format(t2res["kilonovaness"]))
                attributes.append("LVKmap{}".format(t2res["map_name"]))

            # Check whether we have a figure to upload.
            # Assuming this exists locally under {stock}.png
            if self.image_path is not None:
                ipath = self.image_path.replace("ZTFNAME", to_ztf_id(int(tview.id)))
                # Only upload if it actually exists:
                if not os.path.isfile(ipath):
                    ipath = None
            else:
                ipath = None

            payload["ampel_attributes"] = attributes
            self.logger.debug("reacting", extra={"payload": payload})

            # Ok, so we have a transient to react to
            if self.testcollect:
                if post_update:
                    print("Should send an update")
                else:
                    print("New submission")
                print("*** Posting ***")
                print(payload)
                print("*** Done ***")

                # Journal submission
                jcontent = {
                    "reaction": "fake submission",
                    "success": True,
                    "testcollect": True,
                }

            else:
                jcontent = self.colibriclient.firestore_post(payload, image_path=ipath)

            if jcontent:
                tviews.send(JournalAttributes(extra=jcontent))
        return None
