#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2BaseClassifier.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                25.08.2022
# Last Modified Date:  18.09.2024
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Any

import joblib
import lcdata
import matplotlib.pyplot as plt
import numpy as np

# Parsnip models - avoid importing these if running in the light mode?
import parsnip
from astropy.table import Table
from typing_extensions import TypedDict

from ampel.content.NewSVGRecord import NewSVGRecord
from ampel.contrib.hu.model.ParsnipRiseDeclineResult import (
    ClassificationEntry,
    Classifications,
    ParsnipFailure,
    ParsnipPrediction,
    ParsnipResult,
    PredictionEntry,
)
from ampel.model.PlotProperties import PlotProperties
from ampel.model.UnitModel import UnitModel
from ampel.plot.create import create_plot_record

# All parsnip predictions that are not floats
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


# METHODDATA
parsnip_taxonomy = {
    "SNII": 2224,
    "SNIa": 2222,
    "SNibc": 2223,
    "SNIbc": 2223,
    "TDE": 2243,
    "CART": 2245,
    "ILOT": 2244,
    "Mdwarf-flare": 2233,
    "PISN": 2246,
    "KN": 2232,
    "SLSN-I": 2242,
    "SLSN": 2242,
    "SNIa91bg": 2226,
    "SNIax": 2225,
    "dwarf-nova": 2234,
    "uLens": 2235,
}


class ParsnipModelFiles(TypedDict):
    model: str
    classifier: str


class XgbMultiModelFiles(TypedDict):
    path: str
    classes: list[int | str]


def create_parsnip_plot(
    lc_dataset: lcdata.Dataset,
    model: parsnip.ParsnipModel,
    redshift: float,
    properties: PlotProperties,
    extra: dict[str, Any],
) -> NewSVGRecord:
    """
    Create a parsnip fit plot for the provided lightcurve table and model.
    """

    fig, ax = plt.subplots()
    try:
        lc_dataset.light_curves[0].meta["redshift"] = redshift

        parsnip.plot_light_curve(lc_dataset.light_curves[0], model, ax=ax)
        fig.tight_layout()
        return create_plot_record(
            fig,
            properties,
            close=False,
            extra=extra,
        )
    finally:
        fig.clear()
        plt.close(fig)


def run_parsnip_zsample(
    sncosmo_table: Table,
    zs: list[float],
    zweights: list[float] | None,
    *,
    label: str,
    model: parsnip.ParsnipModel,
    classifier: parsnip.Classifier,
    delta_zp: float = 0.0,
    transient_name: str = "noname",
    plot: PlotProperties | None = None,
) -> ParsnipResult | ParsnipFailure:
    """
    Fit a parsnip model for multiple redshifts provided in the `zs` list,
    providing a weighted average result combining `zweights` with the fit chi results.

    delta_zp: the parsnip model zp is set to what is found in the lc table
        with this offset added (for example to reflect different training zp)
    """

    # Adjust zeropoint - needed when training done with inconsitent zeropoints.
    run_zeropoints = set(sncosmo_table["zp"])
    if len(run_zeropoints) > 1:
        raise ValueError("Multiple zeropoints.")
    model.settings["zeropoint"] = run_zeropoints.pop() + delta_zp

    # Create a list of lightcurves, each at a discrete redshift
    lcs = []
    for redshift in zs:
        use_lc = sncosmo_table.copy()
        use_lc.meta["object_id"] = f"parsnip_z{redshift:4f}"
        use_lc.meta["redshift"] = redshift
        lcs.append(use_lc)
    lc_dataset = lcdata.from_light_curves(lcs)
    try:
        lc_predictions = model.predict_dataset(lc_dataset)
    except ValueError:
        # This is to catch a rare error when the phase guess of parsnip fails. As a starting point
        # it looks for the most common 0.1 decimal. Sometimes all datapoints have the same decimal, causing an error ...
        # Looks like a parsnip bug which should be solved there.
        # Should be rare, when only few datapoints, where parsnip anyway does not work well
        return ParsnipFailure(model=label, failed="PhaseGuess")
    lc_classifications = classifier.classify(lc_predictions)

    # Cast result for storage and look at relative probabilities
    predictions: list[PredictionEntry] = []
    classifications: list[ClassificationEntry] = []
    for i, redshift in enumerate(zs):
        prediction = ParsnipPrediction(**lc_predictions[i][lc_predictions.colnames[1:]])
        predictions.append(PredictionEntry(z=redshift, prediction=prediction))
        classifications.append(
            ClassificationEntry(
                z=redshift,
                classification=dict(
                    lc_classifications[i][lc_classifications.colnames[1:]]
                ),
            )
        )

    # Marginalize over the redshift
    # p(c|y) = Integral[p(c|z,y) p(z|y) dz]
    types = lc_classifications.colnames[1:]
    dtype = lc_classifications[types[0]].dtype
    probabilities = lc_classifications[types].as_array().view((dtype, len(types)))
    # Now we could normalize these z prob and normalize types over redshifts
    z_probabilities = np.array([pred.prediction.chi2pdf for pred in predictions])

    if np.sum(z_probabilities) > 0:
        # Take redshift probabilities into account, if available
        if zweights is not None:
            z_probabilities *= zweights
        integrated_probabilities = z_probabilities.dot(probabilities)
        integrated_probabilities /= np.sum(integrated_probabilities)

        # Find the best redshift
        best = np.argmax(z_probabilities)

        return ParsnipResult(
            model=label,
            predictions=predictions,
            classifications=classifications,
            marginal_lc_classifications=dict(
                zip(types, integrated_probabilities, strict=False)
            ),
            z_at_minchi=zs[best],
            # Map result at best redshift to singular value predictions/lc_classifications
            # (wastes DB space, but possible to filter based on)
            prediction=predictions[best].prediction,
            classification=classifications[best].classification,
            plot=create_parsnip_plot(
                lc_dataset,
                model,
                zs[best],
                plot,
                extra={"stock": transient_name, "model": label},
            )
            if plot is not None
            else None,
        )

    # Not enough data for a chi2 estimate or fits with negligable probability
    return ParsnipFailure(model=label, failed="NoFit")


class BaseClassifier:
    """

    Base class for carrying out parsnip and/or xgb classifications.


    Classifiers are expected to be one of:
    - xgboost binary classifier dumped jusing joblib, together with columns to use:
        {'columns':['ndet','t_lc',...],'model':...}
    - xgboost multi classifier...
    - parsnip

    Each are configured by label:path. The label will be used when calling a classifier to use.

    """

    # Setting for report to construct
    classifier_name: str
    classifier_version: str

    ## Paths to classifiers to load.
    paths_xgbbinary: dict[str, str] = {}
    paths_xgbmulti: dict[str, XgbMultiModelFiles] = {}
    paths_parsnip: dict[str, ParsnipModelFiles] = {}

    # Parameters controlling any parsnip fit
    # The parsnip model zeropoint will by default be set to that of
    # the input lightcurve. This can be adjusted by zp offset.
    # Offset can e.g. appear if training performed with wrong zeropoint...
    parsnip_zeropoint_offset: float = 0.0
    # Save / plot parameters
    parsnipplot_properties: PlotProperties | None = None

    add_parsnip_from: None | str = None

    # Include features and ML results in t2_record
    return_features: bool = False

    result_adapter: None | UnitModel = None

    def read_class_models(self) -> None:
        self._class_xgbbinary = {
            label: joblib.load(path) for label, path in self.paths_xgbbinary.items()
        }
        self._class_xgbmulti: dict[str, Any] = {
            label: {**joblib.load(pathdir["path"]), "classes": pathdir["classes"]}
            for label, pathdir in self.paths_xgbmulti.items()
        }
        self._class_parsnip = {
            label: {
                "classifier": parsnip.Classifier.load(paths["classifier"]),
                "model": parsnip.load_model(paths["model"], threads=1),
            }
            for label, paths in self.paths_parsnip.items()
        }

    def post_init(self) -> None:
        self.read_class_models()

    def get_xgb_class(self, classlabel, features):
        """
        Return the classification for the labeled xgb model
        Classify based on features ordered according to columns
        """

        return self._class_xgbbinary[classlabel]["model"].predict_proba(
            np.array(
                [
                    features.get(col, np.nan)
                    for col in self._class_xgbbinary[classlabel]["columns"]
                ]
            ).reshape(1, -1)
        )

    def get_multixgb_class(self, classlabel, features):
        """
        Return the classification for the labeled xgb model
        Classify based on features ordered according to columns
        """

        if self._class_xgbmulti[classlabel]["model"].objective == "multi:softmax":
            # Single type classifier
            modelguess = self._class_xgbmulti[classlabel]["model"].predict(
                np.array(
                    [
                        features.get(col, np.nan)
                        for col in self._class_xgbmulti[classlabel]["columns"]
                    ]
                ).reshape(1, -1)
            )
            pvals = [np.zeros(len(self._class_xgbmulti[classlabel]["classes"]))]
            pvals[0][modelguess] = 1
        elif self._class_xgbmulti[classlabel]["model"].objective == "multi:softprob":
            # Multiple
            pvals = self._class_xgbmulti[classlabel]["model"].predict_proba(
                np.array(
                    [
                        features.get(col, np.nan)
                        for col in self._class_xgbmulti[classlabel]["columns"]
                    ]
                ).reshape(1, -1)
            )

        return {
            str(self._class_xgbmulti[classlabel]["classes"][k]): float(prob)
            for k, prob in enumerate(list(pvals[0]))
        }

    def get_parsnip_class(
        self,
        model_label: str,
        lctable: Table,
        zsample: list[float],
        zweights: list[float] | None = None,
        transient_name: str = "noname",
    ) -> ParsnipResult | ParsnipFailure:
        """
        Return the classification for the labeled parsnip model
        at the redshift list provided by zsample.

        Final classification is weighted mean of fits at sample redshifts,
        weighted by zweights if existing.

        Also limit to max model redshift.
        """

        zw = [1.0 / len(zsample) for z in zsample] if zweights is None else zweights

        return run_parsnip_zsample(
            lctable,
            zsample,
            zw,
            model=self._class_parsnip[model_label]["model"],
            classifier=self._class_parsnip[model_label]["classifier"],
            label=model_label,
            plot=self.parsnipplot_properties,
            transient_name=transient_name,
            delta_zp=self.parsnip_zeropoint_offset,
        )

    def classify(
        self,
        features,
        flux_table,
        redshift_samples,
        transient_name="noname",
        xgb_binarymodels: list | None = None,
        xgb_multimodels: list | None = None,
        parsnip_models: list | None = None,
        redshift_weights=None,
    ) -> list[Classifications]:
        """
        Run the loaded classifiers on the input features / lightcurve.
        'xgb_binarymodels', 'xgb_multimodels', 'parsnip_models' all allow to select which models to use (selected by labels).
        If left empty, all models will be run.

        'add_parsnip_from' points to one parsnip model, for which the (combined) output features are added to the xgb feature collection.
        """

        ### Parsnip
        parsnip_class: dict[str, ParsnipResult | ParsnipFailure] = {}
        if len(flux_table) >= 4:
            models = parsnip_models if parsnip_models else self._class_parsnip.keys()
            if (
                self.add_parsnip_from is not None
                and self.add_parsnip_from not in models
            ):
                raise ValueError("Feature parsnip model not requested to run.")
            parsnip_class = {
                model_name: self.get_parsnip_class(
                    model_name,
                    flux_table,
                    redshift_samples,
                    redshift_weights,
                    transient_name,
                )
                for model_name in models
            }

            ### Update feature dictionary with parsnip output if requested.
            if self.add_parsnip_from is not None and isinstance(
                add_parsnip_from := parsnip_class[self.add_parsnip_from], ParsnipResult
            ):
                # Previously we added parsnip_ as prefix to all features, but this was not done in training.
                #            features.update( {'parsnip_'+featname : featval for featname, featval in  parsnip_class[self.add_parsnip_from]['prediction'].items() } )
                #            features['parsnip_modelchidof'] = features['parsnip_model_chisq'] / features['parsnip_model_dof']
                # Check for feature overlap
                if (
                    len(
                        duplicate := {
                            key: features[key]
                            for key in features
                            if key in add_parsnip_from.prediction.model_dump()
                        }
                    )
                    > 0
                ):
                    raise ValueError(
                        "Found duplicate keys, cannot safely merge feature dictionaries: ",
                        duplicate,
                    )
                features.update(add_parsnip_from.prediction.model_dump())
                features["parsnip_modelchidof"] = (
                    features["model_chisq"] / features["model_dof"]
                )
        else:
            parsnip_class = {}

        ### Run binary XGB classifiers
        # Either run all loaded models or the specified list
        models = xgb_binarymodels if xgb_binarymodels else self._class_xgbbinary.keys()
        binary_class: dict[str, Any] = {
            model_name: self.get_xgb_class(model_name, features)
            for model_name in models
        }

        ### Run multi XGB classifiers
        # Either run all loaded models or the specified list
        models = xgb_multimodels if xgb_multimodels else self._class_xgbmulti.keys()
        multi_class: dict[str, Any] = {
            model_name: self.get_multixgb_class(model_name, features)
            for model_name in models
        }

        ### Prepare output
        return [
            Classifications(
                name=self.classifier_name,
                version=self.classifier_version,
                parsnip=list(parsnip_class.values()),
                xgbbinary=binary_class,
                xgbmulti=multi_class,
                features=features if self.return_features else None,
            )
        ]
