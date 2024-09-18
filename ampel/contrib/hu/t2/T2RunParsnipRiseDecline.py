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
):


    # Name (in case standard) or path to parsnip model to load
    parsnip_model: None | str
    # Path to classifier to apply to lightcurve fit. If not set, no classification will be done.
    parsnip_classifier: None | str 
    # Zeropoint parameters
    # These are separately set in the Parsnip model settings. The zeropoint
    # can does vary between input data, training data and the model.
    # Try to adjust this relative to the 'zp' field of the input tabulated lc
    training_zeropoint: float = 27.5  # Used in Elasticc training sample
    default_zeropoint: float = 25.0  # Default parsnip value
    # Save / plot parameters
    parsnipplot_suffix: None | str = None
    parsnipplot_dir: None | str = None



    def post_init(self) -> None:
        """
        Retrieve models and potentially dustmaps.
        """

        # Load model and classifier
        self.model = parsnip.load_model(self.parsnip_model, threads=1)
        self.classifier = None
        if self.parsnip_classifier:
            self.classifier = parsnip.Classifier.load(self.parsnip_classifier)
        
        # Load Lightcurve extract model (from T2TabulatorRiseDecline)    
        self.init_lightcurve_extractor()

        super().post_init()

    @backoff.on_exception(
        backoff.constant,
        timeout_decorator.timeout_decorator.TimeoutError,
        max_tries=3,
    )
    @timeout_decorator.timeout(5, use_signals=True)
    def _classify_parsnip(self, predictions):
        """
        Carry out the parsnip classification.
        """
        if self.classifier is not None:
            return self.classifier.classify(predictions)
        raise RuntimeError("No classifier configured")

    @backoff.on_exception(
        backoff.constant,
        timeout_decorator.timeout_decorator.TimeoutError,
        max_tries=3,
    )
    @timeout_decorator.timeout(5, use_signals=True)
    def _predict_parsnip(self, dataset):
        """
        Carry out the parsnip predictions.
        """
        return self.model.predict_dataset(dataset)

    def run_parsnip(
        self,
        sncosmo_table: Table,
        z: list[float],
        z_weights: None | list[float],
        stockname: str,
    ) -> dict[str, Any]:
        """
        Run parsnip based on input redshifts.
        
        Will do this based on the list of potential redshifts provided
        by get_redshift in T2DigestRedshifts. These are typically either a catalog redshift, AmpelZ, 
        from catalog matching, or a sample of redshifts.
        fitdatainfo["z"]
        
        """ 
        
        # Initialize result dict
        parsnip_output: dict[str, UBson] = {
            "model": self.parsnip_model,
            "classifier": self.parsnip_classifier,
        }

        

        # Adjust zeropoint
        # Required if model was trained assuming a different zeropoint compared with the input table.
        run_zeropoints = set(sncosmo_table["zp"])
        if len(run_zeropoints) > 1:
            self.logger.info("Warning, multiple zeropoints, using avg.")
            run_zeropoint = np.mean(list(run_zeropoints))
        else:
            run_zeropoint = run_zeropoints.pop()
        self.model.settings["zeropoint"] = (
            self.default_zeropoint + run_zeropoint - self.training_zeropoint
        )
        
        
        # Create a list of lightcurves, each at a discrete redshift
        lcs = []
        for redshift in z:
            use_lc = sncosmo_table.copy()
            use_lc.meta["object_id"] = f"parsnip_z{redshift:4f}"
            use_lc.meta["redshift"] = redshift
            lcs.append(use_lc)
        lc_dataset = lcdata.from_light_curves(lcs)

        lc_predictions = self._predict_parsnip(lc_dataset)
        lc_classifications = self._classify_parsnip(lc_predictions)

        # Cast result for storage and look at relative probabilities
        predictions = {}
        classifcations = {}
        for i, redshift in enumerate(z):
            foo = dict(lc_predictions[i][lc_predictions.colnames[1:]])
            predictions[str(redshift)] = {
                k: dcast_pred[k](v) if k in dcast_pred and v is not None else float(v)
                for k, v in foo.items()
            }
            # Not sure whether the dof could change? Normalizing now
            if foo["model_dof"] > 0:
                predictions[str(redshift)]["chi2pdf"] = chi2.pdf(
                    foo["model_chisq"], foo["model_dof"]
                )
            else:
                # Not enough data - earlier check?
                predictions[str(redshift)]["chi2pdf"] = 0.0
            foo = dict(lc_classifications[i][lc_classifications.colnames[1:]])
            classifcations[str(redshift)] = {
                k: dcast_class[k](v) if k in dcast_class and v is not None else float(v)
                for k, v in foo.items()
            }

        # Marginalize over the redshift
        # p(c|y) = Integral[p(c|z,y) p(z|y) dz]
        types = lc_classifications.colnames[1:]
        dtype = lc_classifications[types[0]].dtype
        probabilities = lc_classifications[types].as_array().view((dtype, len(types)))
        # Now we could normalize these z prob and normalize types over redshifts
        z_probabilities = np.array(
            [lcfit["chi2pdf"] for redshift, lcfit in predictions.items()]
        )

        parsnip_output["predictions"] = predictions
        parsnip_output["classifications"] = classifcations

        if np.sum(z_probabilities) > 0:
            # Take redshift probabilities into account, if available
            if z_weights is not None:
                z_probabilities *= z_weights
            integrated_probabilities = z_probabilities.dot(probabilities)
            integrated_probabilities /= np.sum(integrated_probabilities)
            parsnip_output["marginal_lc_classifications"] = dict(
                zip(types, integrated_probabilities, strict=False)
            )
            # Find the best redshifts
            parsnip_output["z_at_minchi"] = z[np.argmax(z_probabilities)]
            # Map these to singular value predictions/lc_classifications
            # (wastes DB space, but possible to filter based on)
            parsnip_output["prediction"] = predictions[str(parsnip_output["z_at_minchi"])]
            parsnip_output["classification"] = classifcations[str(parsnip_output["z_at_minchi"])]
        else:
            # Not enough data for a chi2 estimate
            parsnip_output["Failed"] = "NoDOF"
            return parsnip_output

        # Plot
        if self.parsnipplot_suffix and self.parsnipplot_dir:

            fig = plt.figure()
            ax = plt.gca()

            # Set redshift to best value and plot this fit
            lc_dataset.light_curves[0].meta["redshift"] = parsnip_output["z_at_minchi"]

            parsnip.plot_light_curve(lc_dataset.light_curves[0], self.model, ax=ax)
            plt.tight_layout()
            plt.savefig(
                os.path.join(self.parsnipplot_dir, f"t2parsnip_{stockname}.{self.parsnipplot_suffix}")
            )

            plt.close("fig")
            plt.close("all")
            del fig
            gc.collect()
            
        return parsnip_output




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
        
        ## Calculate Parnsip features
        z = fitdatainfo["z"]
        z_source = fitdatainfo["z_source"]
        if self.parsnip_model is not None:
            parsnip_features = self.run_parsnip(
                flux_table,
                fitdatainfo["z"],
                fitdatainfo["z_weights"],
                str( compound.get("stock") )
                )
            t2_output['parsnipfeatures'] = parsnip_features

        return t2_output




