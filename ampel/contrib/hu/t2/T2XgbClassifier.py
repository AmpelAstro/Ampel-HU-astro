#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2XgbClassifier.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                03.07.2022
# Last Modified Date:  03.07.2022
# Last Modified By:    jnordin@physik.hu-berlin.de

import os
from collections.abc import Sequence
from typing import Any, Literal, cast

import numpy as np
import xgboost as xgb

from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


# Q: Could this be run as an AbsTiedT2Unit as we do not need the dps?
# class T2XgbClassifier(AbsTiedStateT2Unit):
class T2XgbClassifier(AbsTiedStateT2Unit):
    """
    Load a series of xgboost classifier models (distinguished by number
    of detections) and return a classification.

    Will test whether the first (null) model of the loaded classifier is
    more likely, P(model 0)>0.5. What this means is determined by the training.

    E.g. For the elasticc1v2 model, the null model corresponds to Elasticc class 1
    (and model 1 to class 2.)

    """

    # Detection bins. Note that these reflect _significant_ detections
    # in the sense of TabulatorRiseDecline
    det_ranges: Sequence = [
        [1, 1],
        [2, 2],
        [3, 4],
        [5, 6],
        [7, 9],
        [10, 14],
        [15, 20],
        [21, 30],
        [31, 45],
        [46, 75],
        [76, 1000],
    ]
    model_folder: str = "/home/jnordin/github/zweilasticc/ml_workshop/elasticc/iter100"
    model_prefix: str = "xgboost_elasticc1v2_"
    model_type: str = ".json"

    _classifiers: list[xgb.XGBClassifier]

    # Columns (in order) used for training
    # Which columsn to use for training. But these seem not to align with
    # parquet file!!!!
    #    use_cols: list[str] = ['bool_rise', 'bool_fall', 'bool_peaked', 'bool_pure',
    #       'bool_fastrise', 'bool_fastfall', 'bool_hasgaps', 'mag_det',
    #       'mag_last', 'det_bands', 'peak_bands', 'last_bands', 't_predetect',
    #       't_lc', 't_rise', 't_fall', 'rise_slope_lsstu',
    #       'rise_slopesig_lsstu', 'fall_slope_lsstu', 'fall_slopesig_lsstu',
    #       'rise_slope_lsstg', 'rise_slopesig_lsstg', 'fall_slope_lsstg',
    #       'fall_slopesig_lsstg', 'rise_slope_lsstr', 'rise_slopesig_lsstr',
    #       'fall_slope_lsstr', 'fall_slopesig_lsstr', 'rise_slope_lssti',
    #       'rise_slopesig_lssti', 'fall_slope_lssti', 'fall_slopesig_lssti',
    #       'rise_slope_lsstz', 'rise_slopesig_lsstz', 'fall_slope_lsstz',
    #       'fall_slopesig_lsstz', 'rise_slope_lssty', 'rise_slopesig_lssty',
    #       'fall_slope_lssty', 'fall_slopesig_lssty', 'lsstu-lsstg_det',
    #       'lsstg-lsstr_det', 'lsstr-lssti_det', 'lssti-lsstz_det',
    #       'lsstz-lssty_det', 'lsstu-lsstg_peak', 'lsstg-lsstr_peak',
    #       'lsstr-lssti_peak', 'lssti-lsstz_peak', 'lsstz-lssty_peak',
    #       'lsstu-lsstg_last', 'lsstg-lsstr_last', 'lsstr-lssti_last',
    #       'lssti-lsstz_last', 'lsstz-lssty_last', 'host_sep', 'z', 'z_err',
    #       'band_det_id', 'band_last_id']
    use_cols: list[str] = [
        "ndet",
        "jd_det",
        "jd_last",
        "t_predetect",
        "mag_det",
        "band_det_id",
        "mag_last",
        "band_last_id",
        "t_lc",
        "bool_pure",
        "jd_peak_lsstz",
        "bool_rise",
        "bool_fall",
        "bool_peaked",
        "bool_fastrise",
        "bool_fastfall",
        "bool_hasgaps",
        "det_bands",
        "last_bands",
        "jd_peak_lsstr",
        "rise_slope_lsstz",
        "rise_slopesig_lsstz",
        "fall_slope_lsstz",
        "fall_slopesig_lsstz",
        "jd_peak",
        "t_rise",
        "t_fall",
        "peak_bands",
        "jd_peak_lssty",
        "lsstz-lssty_last",
        "lsstz-lssty_peak",
        "fall_slope_lssty",
        "fall_slopesig_lssty",
        "jd_peak_lssti",
        "fall_slope_lssti",
        "fall_slopesig_lssti",
        "jd_peak_lsstu",
        "fall_slope_lsstr",
        "fall_slopesig_lsstr",
        "jd_peak_lsstg",
        "rise_slope_lsstg",
        "rise_slopesig_lsstg",
        "rise_slope_lsstu",
        "rise_slopesig_lsstu",
        "rise_slope_lssti",
        "rise_slopesig_lssti",
        "lsstu-lsstg_last",
        "fall_slope_lsstu",
        "fall_slopesig_lsstu",
        "lssti-lsstz_last",
        "rise_slope_lssty",
        "rise_slopesig_lssty",
        "fall_slope_lsstg",
        "fall_slopesig_lsstg",
        "lsstg-lsstr_last",
        "lsstr-lssti_last",
        "lsstu-lsstg_peak",
        "rise_slope_lsstr",
        "rise_slopesig_lsstr",
        "lsstg-lsstr_peak",
        "lsstz-lssty_det",
        "lssti-lsstz_peak",
        "lsstr-lssti_peak",
        "lssti-lsstz_det",
        "lsstg-lsstr_det",
        "lsstr-lssti_det",
        "lsstu-lsstg_det",
        "z",
        "z_err",
        "host_sep",
    ]

    # Which units should this be changed to
    t2_dependency: Sequence[
        StateT2Dependency[
            Literal["T2TabulatorRiseDecline", "T2ElasticcRedshiftSampler"]
        ]
    ]

    def post_init(self) -> None:
        """
        Load models
        """

        self._classifiers = []

        for dbin in self.det_ranges:
            fname = os.path.join(
                self.model_folder,
                self.model_prefix + f"ndet{dbin[0]}_{dbin[1]}" + self.model_type,
            )
            model = xgb.XGBClassifier()
            model.load_model(fname=fname)
            self._classifiers.append(model)

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> dict[str, UBson]:
        """

        Extract results from TabulatorRiseDecline and apply suitable model.

        Parameters
        -----------

        t2_records: List of T2Records from (here) TabulatorRiseDecline.

        Returns
        -------
        dict
        """

        t2data: dict[str, Any] = {}
        # Parse t2views - should not be more than one.
        for t2_view in t2_views:
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            # So far only knows how to parse phases from T2TabulatorRiseDecline
            if t2_view.unit == "T2TabulatorRiseDecline":
                t2_res = cast(dict[str, Any], t2_view.get_payload())
                t2data.update(t2_res)
            if t2_view.unit == "T2ElasticcRedshiftSampler":
                t2_res = cast(dict[str, Any], t2_view.get_payload())
                # For some reason we trained xgb using z, zerr and host_sep
                zdata = {"z": None, "z_err": None, "host_sep": None}
                if t2_res["z_source"] in [
                    "HOSTGAL2_ZQUANT",
                    "HOSTGAL_ZQUANT",
                    "HOSTGAL_ZSPEC",
                    "default",
                ]:
                    if len(t2_res["z_samples"]) == 3:
                        # This was the sampling used for job_training
                        zdata["z"] = t2_res["z_samples"][1]
                        zdata["z_err"] = t2_res["z_samples"][1] - t2_res["z_samples"][0]
                        zdata["host_sep"] = t2_res["host_sep"]
                    elif (
                        len(t2_res["z_samples"]) == 4 and t2_res["z_samples"][0] == 0.01
                    ):
                        # ... unless the event was hostless, in which case
                        # xgboost use None values.
                        pass
                    else:
                        self.logger.info(
                            "Do not know how to handle z info", extra=t2_res
                        )
                        print(t2_res)
                        return {"model": self.model_prefix, "xgbsuccess": False}
                t2data.update(zdata)

        if t2data["ndet"] == 0:
            # Explore whether negative detections can be used to parse type
            if "nnegdet" in t2data and t2data["nnegdet"] > 0:
                if t2data["z"] is not None and t2data["z"] > 0.001:
                    # High redshifts are typically AGNs
                    return {
                        "model": self.model_prefix,
                        "xgbsuccess": False,
                        "cause": "No pos det",
                        "direct_eval": "AGN",
                    }
                # Low redshifts are (typically) either uLens or EB.
                # The former have a large fraction of neg det (if any)
                nfrac = float(t2data["nnegdet"]) / t2data["alldet"]
                if nfrac > 0.2:
                    return {
                        "model": self.model_prefix,
                        "xgbsuccess": False,
                        "cause": "No pos det",
                        "direct_eval": "uLens",
                    }
                return {
                    "model": self.model_prefix,
                    "xgbsuccess": False,
                    "cause": "No pos det",
                    "direct_eval": "EB",
                }
            # Otherwise, not much to go by
            return {
                "model": self.model_prefix,
                "xgbsuccess": False,
                "imodel": -1,
                "cause": "No sig. det",
            }
        if not t2data["success"]:
            return {
                "model": self.model_prefix,
                "xgbsuccess": False,
                "imodel": -1,
                "cause": "RiseDecline error.",
            }

        # Find which tree to use
        for b, dbin in enumerate(self.det_ranges):  # noqa: B007
            if t2data["ndet"] >= dbin[0] and t2data["ndet"] <= dbin[1]:
                break

        t2out: dict[str, UBson] = {"imodel": b, "model": self.model_prefix}

        # Can we find away round creating a numpy array?
        prob = self._classifiers[b].predict_proba(
            np.array(
                [[t2data[col] if col in t2data else None for col in self.use_cols]]
            )
        )
        t2out["prob0"] = float(prob[0][0])
        t2out["is_0"] = float(prob[0][0]) > 0.5

        return t2out
