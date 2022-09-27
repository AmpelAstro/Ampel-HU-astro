#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2FastDecliner.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                28.06.2022
# Last Modified Date:  26.09.2022
# Last Modified By:    Jakob Nordin <jnordin@physik.hu-berlin.de>

from astropy.table import Table
from typing import Any, TYPE_CHECKING
from ampel.types import UBson
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.struct.UnitResult import UnitResult
from ampel.view.LightCurve import LightCurve


class T2FastDecliner(AbsLightCurveT2Unit):
    """
    Determine decline rate in two last obs.
    Check whether above some limit.

    :param float min_declinerate: Min decline rate for alarm
    :param float min_trange: Min time diff for consideration (days)
    :param float max_trange: Max time diff for consideration (days)
    :param dict[filter_id, filter_name] filter_keys: Filters to check

    :return: dict, decline rate in each filter and
            bool: fast_decliner
    """

    min_declinerate: float = 0.1
    min_trange: float = 0.1
    max_trange: float = 2.
    filter_keys: dict[str, str] = {"1": "g", "2": "R"}

    def process(self, light_curve: LightCurve) -> UBson | UnitResult:

        # For each filter, check datapoints for latest datapoints
        t2_output = {'fast_decliner': False}
        for filt_id, filt_name in self.filter_keys.items():
            pps = light_curve.get_ntuples( ["jd", "magpsf", "sigmapsf"],
                 filters = [
                    {"attribute": "fid", "operator": "==", "value": int(filt_id)},
                ]
            )
            if pps is None or len(pps)<2:
                continue
            # Simple diffs of two most resent observations
            pps = sorted(pps, key=lambda d: d[0])
            delta_t = pps[-1][0] - pps[-2][0]
            delta_mag = pps[-1][1] - pps[-2][1]
            decline_rate = delta_mag / delta_t

            # Create output structure, annotating if lc "fast declining"
            t2_output[filt_name] = {'delta_t': delta_t, 'decline_rate': decline_rate}
            if (delta_t>self.min_trange and delta_t<self.max_trange and
                    decline_rate>self.min_declinerate):
                t2_output['fast_decliner'] = True

        return t2_output
