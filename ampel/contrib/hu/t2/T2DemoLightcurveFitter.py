#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2DemoLightcurveFitter.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 24.09.2021
# Last Modified Date: 06.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

from collections.abc import Sequence

# The following three only used if correcting for MW dust
import numpy as np

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2BaseLightcurveFitter import T2BaseLightcurveFitter
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson

# from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView


class T2DemoLightcurveFitter(T2BaseLightcurveFitter):
    """

    Demonstration class showing how methods of T2BaseLightcurveFitter can be used
    develop a specific classifier.

    Variables of T2BaseLightcurveFitter determines how to load data and whether to
    try to load a redshift.

    """

    fit_order: int = 3

    def post_init(self) -> None:
        """
        Retrieve models and potentially dustmaps.
        """

        super().post_init()

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """

        Fit a model to the lightcurve of this transient.
        Called for each state of each transient.

        Returns
        -------
        dict
        """

        # Initialize output dict
        t2_output: dict[str, UBson] = {
            "model": "demoFitter",
            "order": self.fit_order,
        }

        # Load photometry as table
        # fitdatainfo contains redshift, if requested
        (sncosmo_table, fitdatainfo) = self.get_fitdata(datapoints, t2_views)
        t2_output.update(fitdatainfo)
        if sncosmo_table is None:
            return t2_output

        # Lightcurve encoded in sncosmo_table can now be fit to model.
        # Fit parameters stored in t2_output and returned.

        for band in set(sncosmo_table["band"]):
            i = sncosmo_table["band"] == band
            if sum(i) <= (self.fit_order + 1):  # Not enough data in band
                continue
            t2_output["polyfit" + band] = list(
                np.polyfit(
                    sncosmo_table[i]["time"] - np.min(sncosmo_table[i]["time"]),
                    sncosmo_table[i]["flux"],
                    deg=self.fit_order,
                )
            )

        return t2_output
