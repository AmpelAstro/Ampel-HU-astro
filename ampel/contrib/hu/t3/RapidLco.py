#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/rapidBase
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                05.08.2019
# Last Modified Date:  20.08.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import datetime
from typing import Any

import requests

from ampel.contrib.hu.t3.RapidBase import RapidBase
from ampel.secret.NamedSecret import NamedSecret
from ampel.util.freeze import recursive_unfreeze
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id


class RapidLco(RapidBase):
    """
    Submit LCO triggers for candidates passing criteria.
    """

    # LCO trigger info
    lco_api: NamedSecret[str] = NamedSecret(label="lco/jnordin")
    # A dict of LCO API triggers to be sent for each SN that fulfills all
    # criteria. Assumed to have the following key content:
    # 'trigger_name': {'start_delay':X (days), 'end_delay':Y (days), 'api_form':Z}
    # Where the start and end delays define the allowed LCO time range and the api_form provides the request to be submitted.
    # The following keys of the api_form will be changed: name, target:name, target:ra, target:dec, windows:end, windows:start
    lco_payload: dict[str, Any] = {
        "lco_u_rapid": {
            "start_delay": 0,
            "end_delay": 1,
            "api_form": {
                "group_id": "ZTF_rapid_sample",
                "proposal": "SUPA2021A-002",
                "ipp_value": 1.05,
                "operator": "SINGLE",
                "observation_type": "RAPID_RESPONSE",
                "requests": [
                    {
                        "acceptability_threshold": 90,
                        "configurations": [
                            {
                                "type": "EXPOSE",
                                "instrument_type": "1M0-SCICAM-SINISTRO",
                                "instrument_configs": [
                                    {
                                        "bin_x": 1,
                                        "bin_y": 1,
                                        "exposure_count": "1",
                                        "exposure_time": "750",
                                        "mode": "full_frame",
                                        "rotator_mode": "",
                                        "extra_params": {"defocus": 0},
                                        "optical_elements": {"filter": "up"},
                                    },
                                    {
                                        "bin_x": 1,
                                        "bin_y": 1,
                                        "exposure_count": "1",
                                        "exposure_time": "60",
                                        "mode": "full_frame",
                                        "rotator_mode": "",
                                        "extra_params": {"defocus": 0},
                                        "optical_elements": {"filter": "gp"},
                                    },
                                ],
                                "acquisition_config": {
                                    "mode": "OFF",
                                    "extra_params": {},
                                },
                                "guiding_config": {
                                    "mode": "ON",
                                    "optional": True,
                                    "extra_params": {},
                                },
                                "target": {
                                    "name": "ZTF_rapid_sample",
                                    "type": "ICRS",
                                    "ra": "x",
                                    "dec": "y",
                                    "proper_motion_ra": 0,
                                    "proper_motion_dec": 0,
                                    "epoch": 2000,
                                    "parallax": 0,
                                },
                                "constraints": {
                                    "max_airmass": 1.6,
                                    "min_lunar_distance": 30,
                                },
                            }
                        ],
                        "windows": [{"end": "x"}],
                        "location": {"telescope_class": "1m0"},
                    }
                ],
            },
        },
        "lco_u_queue": {
            "start_delay": 1,
            "end_delay": 3,
            "api_form": {
                "group_id": "ZTF_rapid_follow",
                "proposal": "SUPA2021A-002",
                "ipp_value": 1.05,
                "operator": "SINGLE",
                "observation_type": "NORMAL",
                "requests": [
                    {
                        "acceptability_threshold": 90,
                        "configurations": [
                            {
                                "type": "EXPOSE",
                                "instrument_type": "1M0-SCICAM-SINISTRO",
                                "instrument_configs": [
                                    {
                                        "bin_x": 1,
                                        "bin_y": 1,
                                        "exposure_count": "2",
                                        "exposure_time": "500",
                                        "mode": "full_frame",
                                        "rotator_mode": "",
                                        "extra_params": {"defocus": 0},
                                        "optical_elements": {"filter": "up"},
                                    },
                                    {
                                        "bin_x": 1,
                                        "bin_y": 1,
                                        "exposure_count": "1",
                                        "exposure_time": "30",
                                        "mode": "full_frame",
                                        "rotator_mode": "",
                                        "extra_params": {"defocus": 0},
                                        "optical_elements": {"filter": "gp"},
                                    },
                                ],
                                "acquisition_config": {
                                    "mode": "OFF",
                                    "extra_params": {},
                                },
                                "guiding_config": {
                                    "mode": "ON",
                                    "optional": True,
                                    "extra_params": {},
                                },
                                "target": {
                                    "name": "ZTF_rapid_follow",
                                    "type": "ICRS",
                                    "ra": "x",
                                    "dec": "x",
                                    "proper_motion_ra": 0,
                                    "proper_motion_dec": 0,
                                    "epoch": 2000,
                                    "parallax": 0,
                                },
                                "constraints": {
                                    "max_airmass": 1.6,
                                    "min_lunar_distance": 30,
                                },
                            }
                        ],
                        "windows": [{"end": "x", "start": "x"}],
                        "location": {"telescope_class": "1m0"},
                    }
                ],
            },
        },
    }

    def react(
        self, tran_view: TransientView, info: None | dict[str, Any]
    ) -> tuple[bool, dict[str, Any]]:
        """
        Send a trigger to the LCO
        """

        if not isinstance(info, dict):
            return False, {"success": False}
        assert isinstance(tran_view.id, int)
        transient_name = to_ztf_id(tran_view.id)

        # Look for coordinates in the T2 info dicts.
        # Assuming ra and dec exists in there
        ra, dec = None, None
        for t2info in info.values():
            if "ra" in t2info.keys():
                ra = t2info["ra"]
            if "dec" in t2info.keys():
                dec = t2info["dec"]
        if ra is None or dec is None:
            # Look at the response
            self.logger.info(
                "No LCO trigger: Could not find ra/dec",
                extra={
                    "target": transient_name,
                },
            )

        # Step through all LCO submit forms
        success = True  # Will be set to false if any submit fails
        submitted = []
        responses = []

        for submit_name, submit_info in self.lco_payload.items():
            # Create submit dictionary
            react_dict = recursive_unfreeze(submit_info["api_form"])

            # Update with information
            react_dict["name"] = submit_name + "_" + transient_name
            react_dict["requests"][0]["configurations"][0]["target"][
                "name"
            ] = transient_name
            react_dict["requests"][0]["configurations"][0]["target"]["ra"] = str(ra)
            react_dict["requests"][0]["configurations"][0]["target"]["dec"] = str(dec)
            # Some keys are not necessarily there
            timenow = datetime.datetime.utcnow()
            if "start" in react_dict["requests"][0]["windows"][0].keys():
                dtime = datetime.timedelta(days=submit_info["start_delay"])
                react_dict["requests"][0]["windows"][0]["start"] = "%s" % (
                    timenow + dtime
                )
            if "end" in react_dict["requests"][0]["windows"][0].keys():
                dtime = datetime.timedelta(days=submit_info["end_delay"])
                react_dict["requests"][0]["windows"][0]["end"] = "%s" % (
                    timenow + dtime
                )

            self.logger.debug(
                "Starting LCO trigger",
                extra={"target": transient_name, "react_dict": react_dict},
            )
            # Make a test to validate
            testreply = requests.post(
                "https://observe.lco.global/api/requestgroups/validate/",
                headers={"Authorization": f"Token {self.lco_api.get()}"},
                json=react_dict,
            )

            # Abort if we have errors
            if len(testreply.json()["errors"]) > 0:
                self.logger.info(
                    "Validating LCO trigger fails for for %s" % (transient_name),
                    extra={"target": transient_name, "react_dict": react_dict},
                )
                success = False
                continue

            # Submit full trigger
            response = requests.post(
                "https://observe.lco.global/api/requestgroups/",
                headers={"Authorization": f"Token {self.lco_api.get()}"},
                json=react_dict,
            )

            # Check whether this was successful
            try:
                response.raise_for_status()
            except requests.exceptions.HTTPError:
                self.logger.info(
                    "Submit LCO fails for %s" % (transient_name),
                    extra={
                        "target": transient_name,
                        "react_dict": react_dict,
                        "response": response.content,
                    },
                )
                success = False

            # Look at the response
            self.logger.info(
                "Submit LCO succeeds for %s" % (transient_name),
                extra={
                    "target": transient_name,
                    "react_dict": react_dict,
                    "response": response.json(),
                },
            )

            submitted.append(react_dict)
            responses.append(response.json())

        # Document what we did
        jcontent = {
            "reactDicts": submitted,
            "success": success,
            "lcoResponses": responses,
        }

        # Note: out-commented because unused
        # jup = JournalAttributes(extra=jcontent)

        return success, jcontent
