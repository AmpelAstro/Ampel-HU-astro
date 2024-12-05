#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2ElasticcKnClassifier.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                22.08.2023
# Last Modified Date:  22.08.2023
# Last Modified By:    jnordin@physik.hu-berlin.de

import copy

from ampel.contrib.hu.t2.T2Elasticc2Classifier import BaseElasticc2Classifier


class T2ElasticcKnClassifier(BaseElasticc2Classifier):
    """

    Classifier module implementing a version of the ZTF/LSST rest
    feature based identification.

    Tries to make use of features available in the RiseDecline feature base set


    Emulate the KN selection made by the Andreoni Rest pipeline, using the features calculated by RiseDecline.
    We assume this to be:
        - t_rise < 2
        - t_fall < 10
       # t-rise cannot always be found among features, will allow this to be skipped for now.

       # Converting to flux limits
       10**( 1 / 2.5) = 2.52
       10**(0.3 / 2.5) = 1.32


        So have to say:
        - require total duration < 10
        - any band, fall_slope < - 1.32
        if bool_peaked:
           any band, rise_slope > 2.5



    """

    # Setting for report to construct
    classifier_name: str = "ElasticcAtRest"
    classifier_version: str = "v231026"

    # Selection limits
    duration_max: float = 10.0
    slope_fall_max: float = -1.32
    slope_rise_min: float = 2.5
    require_rise: bool = False
    rest_bands: list[str] = ["lsstu", "lsstg", "lsstr", "lssti", "lssty", "lsstz"]

    def classify(
        self, base_class_dict, features, flux_table
    ) -> tuple[list[dict], dict]:
        """

        For now a binary classifier, will return:
        - 2232: 1.0 (if alert passes KN selection)
        - 2000: 1.0 (otherwise, saying nothing)

        """

        kn_class = copy.deepcopy(base_class_dict)
        #        kn_class['classifierName'] += 'AtRest'

        # Need minimal nbr of detections but short duration
        # Case I: either no detection OR too long duration
        if features["ndet"] <= 0 or features.get("t_lc", -1) > self.duration_max:
            kn_class["classifications"].append({"classId": 2000, "probability": 1.0})
            return ([kn_class], {"atRest": {"2000": 1.0}})

        falls = [
            features.get("fall_slope_" + band)
            for band in self.rest_bands
            if "fall_slope_" + band in features
        ]
        rise = [
            features.get("rise_slope_" + band)
            for band in self.rest_bands
            if "rise_slope_" + band in features
        ]

        # Case II: either no fall slope measurement, or it is too shallow
        if len(falls) == 0 or min(falls) > self.slope_fall_max:
            # No fall info, or not falling fast enough
            kn_class["classifications"].append({"classId": 2000, "probability": 1.0})
            return ([kn_class], {"atRest": {"2000": 1.0}})

        # Examining rise-time.
        # Case III: we have a measurement of the rise-time slope
        if len(rise) > 0:
            # Case IIIa: rise slope too slow, skip event.
            if self.slope_rise_min > max(rise):
                # Rise measured, too slopw
                kn_class["classifications"].append(
                    {"classId": 2000, "probability": 1.0}
                )
                return ([kn_class], {"atRest": {"2000": 1.0}})
            # Case IIIb: rise slope fast enough, go on the KN classification
        # Case IV: no constraint on rise-time and we require this to exist: non-KN class
        elif self.require_rise:
            # Require rise-time measurement, not available here
            kn_class["classifications"].append({"classId": 2000, "probability": 1.0})
            return ([kn_class], {"atRest": {"2000": 1.0}})

        # Case V: either fast falltime, rise constraint does not exist but this is OK OR case IIIb as above
        # Made it here, Fulfilling both rise and fall requirements
        kn_class["classifications"].append({"classId": 2232, "probability": 1.0})
        return ([kn_class], {"atRest": {"2232": 1.0}})
