#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t0/XShooterFilter.py
# License:             BSD-3-Clause
# Author:              S. Reusch <simeon.reusch@desy.de>
# Date:                07.05.2023
# Last Modified Date:  07.05.2023
# Last Modified By:    S. Reusch

from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.t0.DecentFilter import DecentFilter
from numpy import array


class PredetectionFilter(DecentFilter):
    """
    Filter derived from the DecentFilter.
    It selects new transients without detections prior to a certain date
    (trigger_jd), e.g. for kilonova searches in GW follow up.
    """

    trigger_jd: float

    def post_init(self):
        super().post_init()
        # self.keys_to_check += ("jd",)
        # Is none not working, now doing this manually
        # self.select_upper_limits = [{'attribute': 'magpsf', 'operator': 'is', 'value': None}]
        self.select_photopoints = [
            {"attribute": "magpsf", "operator": "is not", "value": None}
        ]

    def process(self, alert: AmpelAlertProtocol) -> None | bool | int:
        """
        Run the predetection filter on the alert, followed by DecentFilter
        """
        alert_jds = array(alert.get_values("jd", filters=self.select_photopoints))

        predetection_jds = [jd for jd in alert_jds if jd < self.trigger_jd]

        if len(predetection_jds) > 0:
            self.logger.info(
                f"Transient is too old. There are detections prior to trigger time"
            )
            return None

        # now apply the DecentFilter
        return super().process(alert)
