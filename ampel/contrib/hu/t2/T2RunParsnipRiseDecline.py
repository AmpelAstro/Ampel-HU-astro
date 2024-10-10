#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2RunParsnipRiseDecline.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 12.09.2024
# Last Modified Date: 12.09.2024
# Last Modified By  : jnordin@physik.hu-berlin.de

from collections.abc import Iterable, Sequence
from typing import Any

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2BaseClassifier import BaseClassifier
from ampel.contrib.hu.t2.T2BaseLightcurveFitter import T2BaseLightcurveFitter
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import (
    BaseLightCurveFeatures,
    T2TabulatorRiseDeclineBase,
)
from ampel.view.T2DocView import T2DocView

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


def get_probability_evolution(classouts, classtype, classifier, classlabel=None):
    """
    Helpfunction to sort through the output of T2RunParsnipRiseDecline and extract
    evolution of probability fo be of one class.

    Obtain evolution of classified type {classtype} for classifier {classifier}
    (with the optional subkey {lasslabel}) from list of classification outputs {classouts}
    as a function of time.

    Sample usage:
    classlabeling = {
    'Parsnip': {'targetname': 'SNIa', 'model': 'parsnip', 'training': 'snlong' },
    'XRDSampling': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'sample1' },
    'XRDPostpeak': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'later' },
    'XRDEarly': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'early' },
    'XParsnipLast': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'last_parsnip' },
    'XRDLast': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'last_risedec' },
    'XRDPLast': {'targetname': 'snia', 'model': 'xgbmulti', 'training': 'last_prd' },
    }

    results = { modellabel: get_pevo( outs, modelkeys['targetname'], modelkeys['model'], classlabel=modelkeys['training'] )
           for modellabel, modelkeys in classlabeling.items() }

    """
    t, c = [], []
    for allresult in classouts:
        if allresult["classifications"]["features"]["ndet"] == 0:
            continue
        time = allresult["classifications"]["features"]["jd_last"]
        classresult = allresult["classifications"][classifier]
        # Parsnip specific - can this destroy?
        if "Failed" in classresult:
            continue
        if classlabel:
            classresult = classresult[classlabel]
        if "classification" in classresult:
            classresult = classresult["classification"]
        if "classifications" in classresult:
            classresult = classresult["classifications"]
        try:
            c.append(classresult[classtype])
            t.append(time)
        except KeyError:
            pass

    return t, c


class T2RunParsnipRiseDecline(
    T2BaseLightcurveFitter,
    AbsTabulatedT2Unit,
    T2TabulatorRiseDeclineBase,
    BaseLightCurveFeatures,
    BaseClassifier,
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
        (flux_table, fitdatainfo) = self.get_fitdata(list(datapoints), t2_views)

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
        t2_output: dict[str, Any] = {
            "risedeclinefeatures": features,
            "fitdatainfo": fitdatainfo,
        }

        ## Run the classifiers
        t2_output["classifications"] = self.classify(
            features,
            flux_table,
            fitdatainfo["z"],
            redshift_weights=fitdatainfo["z_weights"],
            transient_name=str(compound.get("stock")),
        )

        return t2_output
