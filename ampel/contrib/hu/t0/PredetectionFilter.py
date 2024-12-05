#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t0/XShooterFilter.py
# License:             BSD-3-Clause
# Author:              S. Reusch <simeon.reusch@desy.de>
# Date:                07.05.2023
# Last Modified Date:  06.12.2023
# Last Modified By:    A. Ernst <ernstand@physik.hu-berlin.de>


from numpy import array

from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.t0.DecentFilter import DecentFilter


class PredetectionFilter(DecentFilter):
    """
    Filter derived from the DecentFilter.
    It selects new transients without detections prior to a certain date
    (trigger_jd), e.g. for kilonova searches in GW follow up.
    """

    trigger_jd: float | None = 0
    map_dir: str | None = "./"
    map_name: str | None = None

    def post_init(self):
        super().post_init()
        # self.keys_to_check += ("jd",)
        # Is none not working, now doing this manually
        # self.select_upper_limits = [{'attribute': 'magpsf', 'operator': 'is', 'value': None}]
        self.select_photopoints = [
            {"attribute": "magpsf", "operator": "is not", "value": None}
        ]

        if self.map_name is not None:
            ah = AmpelHealpix(map_name=self.map_name, map_url="", save_dir=self.map_dir)
            # map_hash = ah.process_map()
            self.trigger_jd = ah.get_triggertime()

    def process(self, alert: AmpelAlertProtocol) -> None | bool | int:
        """
        Run the predetection filter on the alert, followed by DecentFilter
        """
        alert_jds = array(alert.get_values("jd", filters=self.select_photopoints))

        predetection_jds = [jd for jd in alert_jds if jd < self.trigger_jd]

        if len(predetection_jds) > 0:
            self.logger.debug(
                "Transient is too old. There are detections prior to trigger time"
            )
            return None

        # now apply the DecentFilter
        return super().process(alert)
