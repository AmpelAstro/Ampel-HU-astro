#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2GuessSN.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                06.04.2020
# Last Modified Date:  03.08.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import numpy as np
from typing import Union
from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.view.LightCurve import LightCurve
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.contrib.hu.t2.T2RiseDeclineStat import T2RiseDeclineBase
import ampel.contrib.hu.t2.xgb_trees as xgb_trees


class T2BrightSNProb(AbsLightCurveT2Unit, T2RiseDeclineBase):
    """
    Derive a number of simple metrics describing the rise, peak and decline of a lc.
    Run a XGB tree trained to check whether this transient are likely to be an RCF SN.
    """

    def post_init(self):
        # Files and selection parameters for tree. Dict keys correspond to nbr of detections.
        # The other parameters were requirements for inclusion in training (and thus usage)
        self.model_params = [
            "mag_det",
            "mag_last",
            "t_lc",
            "rb_med",
            "col_det",
            "t_predetect",
            "distnr_med",
            "magnr_med",
            "classtar_med",
            "sgscore1_med",
            "distpsnr1_med",
            "neargaia_med",
            "maggaia_med",
        ]
        self.xgb_tree_param = {
            2: {"max_duration": 3.5, "max_predetect": 3.5, "min_detmag": 16},
            3: {"max_duration": 6.5, "max_predetect": 3.5, "min_detmag": 16},
            4: {"max_duration": 6.5, "max_predetect": 3.5, "min_detmag": 16},
            5: {"max_duration": 10, "max_predetect": 3.5, "min_detmag": 16},
            6: {"max_duration": 10, "max_predetect": 3.5, "min_detmag": 16},
            100: {"max_duration": 90, "max_predetect": 10, "min_detmag": 16},
        }
        # Load the (large) set of trees
        self.xgb_tree = xgb_trees.xgboost_tree()


    def process(self, light_curve: LightCurve) -> Union[UBson, UnitResult]:

        # Output dict that we will start to populate
        o = self.compute_stats(light_curve)

        if not o["success"]:
            return o

        # This is where the data collection stops and evaluation starts.
        # Even though not all of the properties above are used we keep them for future compatibility.

        # Did not train for gap objects or objects with intervening upper limits
        if o["bool_hasgaps"] or not o["bool_pure"]:
            self.logger.info(
                "Abort due to LC quality",
                extra={"bool_hasgaps": o["bool_hasgaps"], "bool_pure": o["bool_pure"]},
            )
            o["SNGuess"] = 0.5  # Model could not run
            o["SNGuessBool"] = None
            o["success"] = False
            return o

        # Determine which lightcurve requirements to use
        if o["ndet"] < 2 or o["ndet"] > 100:
            self.logger.info("LC size not covered", extra={"nbr_det": o["ndet"]})
            o["SNGuess"] = 0.5  # Model could not run
            o["SNGuessBool"] = None
            o["success"] = False
            return o

        # Which tree param key to use?
        if o["ndet"] < 7:
            treekey = o["ndet"]
        elif o["ndet"] <= 100:
            treekey = 100

        # One large bin for all longer lightcurves at this time
        max_lc_time = self.xgb_tree_param[treekey]["max_duration"]
        max_predetect_time = self.xgb_tree_param[treekey]["max_predetect"]
        min_det_mag = self.xgb_tree_param[treekey]["min_detmag"]

        # Verify lightcurve properties
        if (
            o["t_lc"] > max_lc_time
            or o["t_predetect"] is None
            or o["t_predetect"] > max_predetect_time
            or o["mag_det"] < min_det_mag
        ):
            self.logger.info(
                "LC prop outside model training",
                extra={
                    "t_lc": o["t_lc"],
                    "t_predetect": o["t_predetect"],
                    "mag_det": o["mag_det"],
                },
            )
            o["SNGuess"] = 0.5  # Model could not run
            o["SNGuessBool"] = None
            o["success"] = False
            return o

        # The xgb tree cannot handle None, but np.nan works
        nandict = {k: (np.nan if v is None else v) for (k, v) in o.items()}
        fitresult = self.xgb_tree.predict_risedecline(nandict, treekey)
        o["SNGuess"] = fitresult

        if fitresult > 0.5:
            o["SNGuessBool"] = 1
        else:
            o["SNGuessBool"] = 0

        # In either case the fit succeeded.
        o["success"] = True
        return o
