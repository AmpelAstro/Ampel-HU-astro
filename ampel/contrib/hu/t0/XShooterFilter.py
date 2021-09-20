#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/XShooterFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 28.08.2018
# Last Modified Date: 19.04.2021
# Last Modified By  : jnordin 


from typing import Optional, Sequence, Union

from astropy.time import Time
from numpy import array

from ampel.alert.PhotoAlert import PhotoAlert
from ampel.ztf.t0.DecentFilter import DecentFilter


class XShooterFilter(DecentFilter):
    """
    Filter derived from the DecentFilter, in addition selecting very new
    transients which are visible from the South. In particular, the transient
    are accepted if they are detected during the last 6h, at least one non-detection
    during the last 5 days (and no detection during this time).
    """

    max_dec: float  # maximum allowed value for the declination
    det_within: float  # the transient must have been detected within the last 'DET_WITHIN' days
    ul_within: float  # the transient must AT LEAST one ulim within the last 'UL_WITHIN' days
    # Updated parameters based on infant detections spring 2021. Defaults conservative
    max_chipsf: float = 4        # Best guess value 2
    max_seeratio: float = 2      # Best guess value 1.3
    min_sumrat: float = 0.6      # Best guess value 0.8

    def post_init(self):

        super().post_init()
        self.keys_to_check += ("jd",)

    # Override
    def process(self, alert: PhotoAlert) -> Optional[Union[bool, int]]:
        """
        run the decent filter on the alert
        """

        # cut on declination
        latest = alert.pps[0]
        if latest["dec"] > self.max_dec:
            self.logger.debug(
                f"Rejected: declination {latest['dec']:.2f} deg is "
                f"above maximum allowed of {self.max_dec:.2f} deg"
            )
            return None

        # CUT ON LATEST SUBTRACTION QUALITY
        #################################

        if latest["chipsf"] > self.max_chipsf:
            self.logger.debug(
                f"Rejected: chipsf {latest['chipsf']:.2f}  "
                f"above maximum allowed of {self.max_chipsf:.2f}"
            )
            return None
        if latest["seeratio"] > self.max_seeratio:
            self.logger.debug(
                f"Rejected: seeratio {latest['seeratio']:.2f}  "
                f"above maximum allowed of {self.max_seeratio:.2f}"
            )
            return None
        if latest["sumrat"] < self.min_sumrat:
            self.logger.debug(
                f"Rejected: sumrat {latest['sumrat']:.2f}  "
                f"below min allowed of {self.min_sumrat:.2f}"
            )
            return None


        # CUT ON THE HISTORY OF THE ALERT
        #################################

        now_jd = latest["jd"]
        self.logger.debug(f"Setting 'now' to JD {now_jd:.4f} to cut on alert history")

        # check on history 1: detected in the last 6h
        detection_jds = array(alert.get_values("jd"))
        recent_detections = detection_jds > (now_jd - self.det_within)
        if not any(recent_detections):
            self.logger.debug(
                f"Rejected: no detection within the last {self.det_within:.3f} days "
                f"(latest one {(now_jd - max(detection_jds)):.3f} days ago)"
            )
            return None

        # check on the history 2: at least one upper limit in the last 5 days
        ulim_jds = alert.get_values("jd", data="uls")
        if ulim_jds is None:
            self.logger.debug("Rejected: this alert has no upper limits")
            return None

        if not any(array(ulim_jds) > (now_jd - self.ul_within)):
            self.logger.debug(
                f"Rejected: no upper limit in the last {self.ul_within:.3f} days"
            )
            return None

        # check on history 3: no detection within the last 5 days
        not_so_recent_detections = detection_jds[~recent_detections]
        if any(not_so_recent_detections > (now_jd - self.ul_within)):
            self.logger.debug(
                f"Rejected: at least one detection within the last {self.ul_within:.3f} days "
                f"(but not within {self.det_within:.3f} days)."
            )
            return None

        # now apply the DecentFilter
        return super().process(alert)
