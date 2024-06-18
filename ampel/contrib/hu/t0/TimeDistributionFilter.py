#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t0/TimeDistributionFilter.py
# License:             BSD-3-Clause
# Author:              jnordin <jnordin@desy.de>
# Date:                28.08.2018
# Last Modified Date:  24.11.2021
# Last Modified By:    jnordin

import numpy as np

from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.t0.DecentFilter import DecentFilter


class TimeDistributionFilter(DecentFilter):
    """
    Filter derived from the DecentFilter, in addition:

    For alerts accepted by the DecentFilter, will also look at the distribution of times for detections.

    Main goal is to reject spurioius early / late single detections and then check that the remaining
    detections fulfill criteria.

    """

    min_magfromlim: float = 0.1  # Min distance between mag lim and (ap) photometry to include in time discussion
    jd_rej_sigma: float = 5.0  # Reject dates outside this range

    min_masked_duration: float = 5.0
    max_masked_duration: float = 360.0

    def post_init(self):
        super().post_init()

    # Override
    def process(self, alert: AmpelAlertProtocol) -> None | bool | int:
        """
        run  filter on provided alert

        Return 0 / False if alert is to be rejected.
        Return True / >0 if alert eccepted
        """

        # Apply the DecentFilter
        if (is_decent := super().process(alert)) is False:
            self.logger.info("rejected by decent")
            return is_decent

        # Determine photopoints to use for time study
        def pp_significant(pp):
            # Has to be a candidate
            if pp.get("candid") is None:
                return False
            # Has to be positive
            if pp.get("isdiffpos") in ["f", "0"]:
                return False
            if pp.get("magfromlim") < self.min_magfromlim:
                return False

            # PP ok
            return True

        # Mask based on unique dates - will prevent domination by high cadence observations. 
        jds = np.unique( [int(pp.get("jd")) for pp in alert.datapoints if pp_significant(pp)] )        
        mask = jds > 0
        t_median = np.median( jds[mask] )

        # Iteratively reject datapoints until we reach a minimum (or find no more)
        # Using the median absolute deviation
        while sum(mask) >= 1:
            sig_est = 1.48 * np.median(np.abs(jds[mask] - t_median))
            new_mask = np.abs(jds - t_median) < self.jd_rej_sigma * sig_est
            if sum(new_mask) < sum(mask):
                mask = new_mask
                t_median = np.median(jds[mask])
            else:
                break

        # We are here looking at dates of observatiion, so would be confusing to treat that as detections.
        #if sum(mask) < self.min_ndet:
        #    return False

        t_masked_duration = (
            np.max(jds[mask]) - np.min(jds[mask]) if sum(mask) > 0 else 0
        )
        if (
            t_masked_duration < self.min_masked_duration
            or t_masked_duration > self.max_masked_duration
        ):
            return False

        return True
