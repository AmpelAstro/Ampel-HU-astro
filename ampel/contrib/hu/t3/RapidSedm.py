#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/RapidSedm.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 05.08.2019
# Last Modified Date: 06.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import datetime
from typing import Any, Dict, List, Optional, Tuple

import requests

from ampel.contrib.hu.t3.RapidBase import RapidBase, Secret
from ampel.struct.JournalTweak import JournalTweak
from ampel.type import StockId
from ampel.view.TransientView import TransientView
from ampel.ztf.utils import to_ampel_id, to_ztf_id


class RapidSedm(RapidBase):
    """
    Select transients for rapid reactions. Intended as base class where the react method can be
    implemented as wished and a testreact method posts test reactions to Slack.

    This version reacts by setting a target for SEDM observatoins
    """

    # Base SEDM trigger info
    sedm_url: str = "http://pharos.caltech.edu/request"
    sedm_payload: Dict[str, Any] = {
        "obj_ra": None,
        "obj_dec": None,
        "obj_epoch": 2000,
        "obj_mag": None,
        "obj_name": None,
        "allocation": 20180319205302725,
        "status": "PENDING",
        "inidate": None,
        "enddate": None,
        "priority": 5.0,
        "ifu": True,
        "ab": "n",
        "ifu_use_mag": "y",
        "rc": False,
        "rc_use_mag": "y",
        "do_r": "y",
        "r_exptime": 0,
        "r_repeats": 1,
        "do_g": "y",
        "g_exptime": 0,
        "g_repeats": 1,
        "do_i": "y",
        "i_exptime": 0,
        "i_repeats": 1,
        "do_u": "y",
        "u_exptime": 0,
        "u_repeats": 1,
        "maxairmass": 2.5,
        "min_moon_dist": 30,
        "user_id": 284,
    }

    sedm_username: str
    sedm_password: Secret[str]

    # maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    max_redshift: float = 0.05
    # arcsec, maximum distance
    max_dist: float = 30.0
    # range of peak magnitudes for submission
    min_peak_mag: float = 19.25
    # Minimal median RB.
    drb_minmed: float = 0.995
    # Limiting magnitude to consider upper limits as 'significant'
    maglim_min: float = 19.25

    def post_init(self):
        """"""
        self.name = "RapidSedm"
        super().post_init()

    def react(
        self, tran_view: TransientView, info: Dict[str, Any]
    ) -> Tuple[bool, Optional[JournalTweak]]:
        """
        Send a trigger to the SEDM. Note that we have no good way of investigating the queue at this time
        """
        assert isinstance(tran_view.id, int)
        # Assemble required information. These *should* already be present in the default info
        # provided by the info dict returned for a sucessfull accept_tview
        react_dict = {}
        react_dict.update(self.sedm_payload)
        react_dict["obj_ra"] = info["ra"]
        react_dict["obj_dec"] = info["dec"]
        react_dict["obj_mag"] = info[
            "latest_mag"
        ]  # Assuming that the object is not declining?
        react_dict["obj_name"] = to_ztf_id(tran_view.id)
        react_dict["inidate"] = datetime.datetime.utcnow()
        react_dict["enddate"] = datetime.datetime.utcnow() + datetime.timedelta(days=2)

        # We are still in debug stage, turn down priority
        # react_dict['priority'] = 1

        self.logger.debug(
            "SEDM trigger for %s w dict %s" % (to_ztf_id(tran_view.id), react_dict)
        )

        # Make the post
        response = requests.post(
            self.sedm_url,
            data=react_dict,
            auth=(self.sedm_username, self.sedm_password.get()),
        )

        # Check result
        success = response.status_code == 200

        # Document what we did
        jcontent = {"t3unit": self.name, "reactDict": react_dict, "success": success}
        jup = JournalTweak(extra=jcontent)

        return success, jup

    def add(self, transients: Tuple[TransientView, ...]) -> Dict[StockId, JournalTweak]:
        """
        Loop through transients and check for TNS names and/or candidates to submit
        """

        journal_updates: Dict[StockId, JournalTweak] = {}
        # We will here loop through transients and react individually
        for tv in transients:
            matchinfo = self.accept_tview(tv)

            # Check sumission criteria
            if not matchinfo:
                continue

            self.logger.info("Passed reaction threshold", extra={"tranId": tv.id})

            # Ok, so we have a transient to react to
            if self.do_react:
                success, jup = self.react(tv, matchinfo)
                if jup is not None:
                    journal_updates[tv.id] = jup
                if success:
                    self.logger.info(
                        "React success",
                        extra={"tranId": tv.id, "success": success},
                    )
                else:
                    self.logger.info(
                        "React failure",
                        extra={"tranId": tv.id, "success": success},
                    )
            else:
                success = False
                jup = None

            # Otherwise, test
            matchinfo["SEDM trigger success"] = success
            if self.do_testreact:
                test_success, jup = self.test_react(tv, matchinfo)
                if jup is not None:
                    journal_updates[tv.id] = jup

        return journal_updates

    def done(self):
        """ """

        # Should possibly do some accounting or verification
        self.logger.info("done running T3")
