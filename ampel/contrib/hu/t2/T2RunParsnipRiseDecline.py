#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2RunParsnipRiseDecline.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 12.09.2024
# Last Modified Date: 12.09.2024
# Last Modified By  : jnordin@physik.hu-berlin.de

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Iterable	

import os
import gc
import numpy as np
import timeout_decorator
import backoff
from astropy.table import Table
from scipy.stats import chi2
import matplotlib.pyplot as plt

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.view.T2DocView import T2DocView
from ampel.contrib.hu.t2.T2BaseLightcurveFitter import T2BaseLightcurveFitter
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import T2TabulatorRiseDeclineBase
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import BaseLightCurveFeatures
from ampel.contrib.hu.t2.T2BaseClassifier import BaseClassifier

import parsnip
import lcdata

# Correct types for parsnip properties
dcast_pred = {
    "object_id": str,
    "type": str,
    "count": int,
    "count_s2n_3": int,
    "count_s2n_5": int,
    "count_s2n_3_pre": int,
    "count_s2n_3_rise": int,
    "count_s2n_3_post": int,
    "model_dof": int,
}
dcast_class = {
    "object_id": str,
}



class T2RunParsnipRiseDecline(
    T2BaseLightcurveFitter,
    AbsTabulatedT2Unit,
    T2TabulatorRiseDeclineBase,
    BaseLightCurveFeatures,
    BaseClassifier
):





    def post_init(self) -> None:
        """
        Retrieve models and potentially dustmaps.
        """

        # Load model and classifiers (from T2BaseClassifier)
        self.read_class_models()
        
        # Load Lightcurve extract model (from T2TabulatorRiseDecline)    
        self.init_lightcurve_extractor()

        super().post_init()



    def process(
        self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> dict[str, Any]:
        """
        Process datapoints belonging to one state of one transient.
        A commong Table is generated which is used as input
        to the feature generator.

    Derive features from:
    - TabulatorRiseDecline
    - Parsnip
    
    combine features and use tese as input for XGBoost classification.


    BaseLightcurveFitter inherits from T2DigestRedshifts which is a AbsTiedLightCurveT2Unit unit, accepting a LightCurve as process input.
    T2TabulatorRiseDecline inheirts from AbsStateT2Unit, expecting a T1Document as input. 
    Logical might be to change the latter unit, but this is now incorporated into e.g. T2RunSncosmo, so would require all of these to be changed. 
    So fastest solution is probably to rewrite and BaseLightcurveFitter...


        """

        # Convert input datapoints to standardized Astropy Table
        # Using standard tabulators
        # Will also correct for MW extinction, if chosen.
        # fitdatainfo contains redshift, if requested
        (flux_table, fitdatainfo) = self.get_fitdata(datapoints, t2_views)
        
        if flux_table is None:
            return {"fitdatainfo": fitdatainfo}


        ## Calculate lightcurve features
        # Obtain RiseDecline features
        features = self.compute_stats(flux_table)
        # Calculate light_curve features
        lcfeat = self.extract_lightcurve_features(flux_table)
        features.update(lcfeat)
        # Create averaged values
        avgfeat = self.average_filtervalues(features)
        features.update(avgfeat)
        t2_output = { "risedeclinefeatures": features, "fitdatainfo": fitdatainfo }
        
        ## Run the classifiers
        t2_output["classifications"] = self.classify( 
            features, 
            flux_table,
            fitdatainfo["z"], 
            redshift_weights=fitdatainfo["z_weights"], 
            transient_name=str(compound.get("stock")) )


        return t2_output




