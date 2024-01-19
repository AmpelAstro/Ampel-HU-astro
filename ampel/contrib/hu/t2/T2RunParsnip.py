#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2RunParsnip.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 24.09.2021
# Last Modified Date: 06.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

import gc
import os
import warnings
from collections.abc import Sequence
from typing import Literal

import backoff
import matplotlib.pyplot as plt
import numpy as np
import packaging
import scipy
import timeout_decorator
from astropy.table import Table
from scipy.stats import chi2

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView

# do not warning about scipy.stats.mode(keepdims=None)
if packaging.version.parse(scipy.__version__) < packaging.version.parse("1.11"):
    warnings.filterwarnings(
        "ignore", category=FutureWarning, module="parsnip.light_curve", lineno=31
    )

import extinction  # type: ignore[import]
import lcdata
import parsnip
import sncosmo  # type: ignore[import]

# The following three only used if correcting for MW dust
from sfdmap2.sfdmap import SFDMap  # type: ignore[import]

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


class T2RunParsnip(AbsTiedStateT2Unit, AbsTabulatedT2Unit):
    """
    Gathers information and runs the parsnip model and classifier.
    - Obtain model (read from file unless not in sncosmo registry)
    - Parse previous (chained) T2results for redshift or phase limits.
    - Creates suitable photometry, using the converter provided and  phase limits.
    - Defines model appropritately, including fit ranges and fixed parameters.
    - Run fit, potentially iterative in case of non-convergence.
    - Plot output if requested

    TODO:
    - Add option for redoing fits with disturbed initial conditions to avoid local minima
    - Add option for masking data?
    """

    # Name (in case standard) or path to parsnip model to load
    parsnip_model: str
    # Path to classifier to apply to lightcurve fit. If not set, no classification will be done.
    parsnip_classifier: None | str

    # Redshift usage options. Current options
    # T2MatchBTS : Use the redshift published by BTS and  synced by that T2.
    # T2DigestRedshifts : Use the best redshift as parsed by DigestRedshift.
    # T2ElasticcRedshiftSampler: Use a list of redshifts and weights from the sampler.
    # None : run sncosmo template fit with redshift as free parameter OR use backup_z if set
    redshift_kind: None | str
    # If loading redshifts from DigestRedshifts, provide the max ampel z group to make use of.
    # (note that filtering based on this can also be done for a potential t3)
    max_ampelz_group: int = 3
    # It is also possible to use fixed redshift whenever a dynamic redshift kind is not possible
    # This could be either a single value or a list
    fixed_z: None | float | Sequence[float]
    # Finally, the provided lens redshift might be multiplied with a scale
    # Useful for lensing studies, or when trying multiple values
    scale_z: None | float

    max_fit_z: None | float = None

    # Remove MW dust absorption.
    # MWEBV is either derived from SFD maps using the position from light_curve
    # (assuming the SFD_DIR env var is set) OR retrieved from stock (ELASTICHOW?)
    # The default value of Rv will be used.
    apply_mwcorrection: bool = False

    # Further fit parameters
    # Bounds - not yet implemented
    # sncosmo_bounds : Dict[ str, List[float] ] = {}
    # Remove MW dust absorption using SFD maps.
    # assumes that the position can be retrieved from the light_curve and
    # that the SFD_DIR env var is set to allow them to be found.
    # The default value of Rv will be used.
    # apply_mwcorrection : bool = False

    # Phase range usage. Current option:
    # T2PhaseLimit : use the jdmin jdmax provided in this unit output
    # None : use full datapoint range
    # (T2BayesianBlocks should be added)
    phaseselect_kind: None | str

    # Abort veto (if fulfilled, skip run)
    abort_map: None | dict[str, list]

    # Zeropoint parameters
    # These are separately set in the Parsnip model settings. The zeropoint
    # can does vary between input data, training data and the model.
    # Try to adjust this relative to the 'zp' field of the input tabulated lc
    training_zeropoint: float = 27.5  # Used in Elasticc training sample
    default_zeropoint: float = 25.0  # Default parsnip value

    # Save / plot parameters
    plot_suffix: None | str
    plot_dir: None | str

    # Which units should this be changed to
    t2_dependency: Sequence[
        StateT2Dependency[
            Literal[
                "T2ElasticcRedshiftSampler",
                "T2DigestRedshifts",
                "T2MatchBTS",
                "T2PhaseLimit",
                "T2XgbClassifier",
            ]
        ]
    ]

    def post_init(self) -> None:
        """
        Retrieve models.
        """

        # Load model and classifier
        self.model = parsnip.load_model(self.parsnip_model, threads=1)
        self.classifier = None
        if self.parsnip_classifier:
            self.classifier = parsnip.Classifier.load(self.parsnip_classifier)

        if self.apply_mwcorrection:
            self.dustmap = SFDMap()

    def _get_redshift(
        self, t2_views
    ) -> tuple[None | list[float], None | str, None | list[float]]:
        """
        Can potentially also be replaced with some sort of T2DigestRershift tabulator?

        Assuming that only one instance of redshift sources exist
        """

        # Examine T2s for eventual information
        z: None | list[float] = None
        z_source: None | str = None
        z_weights: None | list[float] = None

        if self.redshift_kind in [
            "T2MatchBTS",
            "T2DigestRedshifts",
            "T2ElasticcRedshiftSampler",
        ]:
            for t2_view in t2_views:
                if not t2_view.unit == self.redshift_kind:
                    continue
                self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )
                # Parse this
                if self.redshift_kind == "T2MatchBTS":
                    if (
                        "bts_redshift" in t2_res.keys()
                        and not t2_res["bts_redshift"] == "-"
                    ):
                        z = [float(t2_res["bts_redshift"])]
                        z_source = "BTS"
                elif self.redshift_kind == "T2DigestRedshifts":
                    if (
                        "ampel_z" in t2_res.keys()
                        and t2_res["ampel_z"] is not None
                        and t2_res["group_z_nbr"] <= self.max_ampelz_group
                    ):
                        z = [float(t2_res["ampel_z"])]
                        z_source = "AMPELz_group" + str(t2_res["group_z_nbr"])
                elif self.redshift_kind == "T2ElasticcRedshiftSampler":
                    z = t2_res["z_samples"]
                    z_source = t2_res["z_source"]
                    z_weights = t2_res["z_weights"]
        else:
            # Check if there is a fixed z set for this run, otherwise keep as free parameter
            if self.fixed_z is not None:
                if isinstance(self.fixed_z, float):
                    z = [self.fixed_z]
                else:
                    z = list(self.fixed_z)
                z_source = "Fixed"
            else:
                z = None
                z_source = "Fitted"

        if (z is not None) and (z_source is not None) and self.scale_z:
            z = [onez * self.scale_z for onez in z]
            z_source += f" + scaled {self.scale_z}"

        return z, z_source, z_weights

    def _get_abort(self, t2_views) -> tuple[bool, dict]:
        """
        Check potential previous t2s for whether the run should be aborted.
        (For perfomance reasons).

        Implemented case is concerns T2XgbClassifier.
        """

        if not self.abort_map or len(self.abort_map) == 0:
            # Not looking for any
            return (False, {})

        abort, abort_maps = False, {}
        for t2_view in t2_views:
            if t2_view.unit not in self.abort_map.keys():
                continue
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
            abort_maps.update(t2_res)

            for abort_map in self.abort_map[t2_view.unit]:
                if all(t2_res.get(key, None) == val for key, val in abort_map.items()):
                    abort = True

        return (abort, abort_maps)

    def _get_phaselimit(self, t2_views) -> tuple[None | float, None | float]:
        """
        Can potentially also be replaced with some sort of tabulator?

        """

        # Examine T2s for eventual information
        jdstart: None | float = None
        jdend: None | float = None

        if self.phaseselect_kind is None:
            jdstart = -np.inf
            jdend = np.inf
        else:
            for t2_view in t2_views:
                # So far only knows how to parse phases from T2PhaseLimit
                if not t2_view.unit == "T2PhaseLimit":
                    continue
                self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )
                jdstart = t2_res["t_start"]
                jdend = t2_res["t_end"]

        return jdstart, jdend

    def _deredden_mw_extinction(self, ebv, phot_tab, rv=3.1) -> Table:
        """
        For an input photometric table, try to correct for mw extinction.
        """

        # Find effective wavelength for all filters in phot_tab
        filterlist = set(phot_tab["band"])
        eff_wave = [sncosmo.get_bandpass(f).wave_eff for f in filterlist]

        # Determine flux correction (dereddening) factors
        flux_corr = 10 ** (0.4 * extinction.ccm89(np.array(eff_wave), ebv * rv, rv))

        # Assign this appropritately to Table
        phot_tab["flux_original"] = phot_tab["flux"]
        phot_tab["fluxerr_original"] = phot_tab["fluxerr"]
        for k, band in enumerate(filterlist):
            phot_tab["flux"][(phot_tab["band"] == band)] *= flux_corr[k]
            phot_tab["fluxerr"][(phot_tab["band"] == band)] *= flux_corr[k]

        return phot_tab

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

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """

        Fit the parameters of the initiated snocmo_model to the light_curve
        provided. Depending on the configuration, the provided T2DovViews
        are used to look for redshift information and any phase (time)
        limits for the fit.


        Parameters
        -----------
        light_curve: "ampel.view.LightCurve" instance.
        See the LightCurve docstring for more info.

        t2_records: List of T2Records from the following units (if available)
        T2DigestRedshifts (redshift parsed from catalogs)
        T2MatchBTS (redshifts synced from BTS page)
        T2PhaseLimit (fit time-limits as determined from lightcurve)

        Returns
        -------
        dict
        """

        # Initialize output dict
        t2_output: dict[str, UBson] = {
            "model": self.parsnip_model,
            "classifier": self.parsnip_classifier,
        }

        # Check whether no computation should be done (due to previous fit)
        (abort_run, abort_info) = self._get_abort(t2_views)
        t2_output["abort_maps"] = abort_info
        if abort_run:
            t2_output["aborted"] = True
            return t2_output

        # Check for phase limits
        (jdstart, jdend) = self._get_phaselimit(t2_views)
        t2_output["jdstart"] = jdstart
        t2_output["jdend"] = jdend
        if t2_output["jdstart"] is None:
            return t2_output

        # Obtain photometric table
        sncosmo_table = self.get_flux_table(datapoints)
        sncosmo_table = sncosmo_table[
            (sncosmo_table["time"] >= jdstart) & (sncosmo_table["time"] <= jdend)
        ]
        self.logger.debug(f"Sncosmo table {sncosmo_table}")

        # Adjust zeropoint - does this matter? and should we have changed it?
        run_zeropoints = set(sncosmo_table["zp"])
        if len(run_zeropoints) > 1:
            self.logger.info("Warning, multiple zeropoints, using avg.")
            run_zeropoint = np.mean(list(run_zeropoints))
        else:
            run_zeropoint = run_zeropoints.pop()
        self.model.settings["zeropoint"] = (
            self.default_zeropoint + run_zeropoint - self.training_zeropoint
        )

        # Potentially correct for dust absorption
        if self.apply_mwcorrection:
            # Get ebv from coordiantes.
            # Here there should be some option to read it from journal/stock etc
            mwebv = self.dustmap.ebv(*self.get_pos(datapoints, which="mean"))
            t2_output["mwebv"] = mwebv
            sncosmo_table = self._deredden_mw_extinction(mwebv, sncosmo_table)

        ## Obtain redshift(s) from catalog fit or a RedshiftSample
        z, z_source, z_weights = self._get_redshift(t2_views)
        t2_output["z"] = z
        t2_output["z_source"] = z_source
        t2_output["z_weights"] = z_weights
        # A source class of None indicates that a redshift source was required, but not found.
        if z is None or z_source is None:
            return t2_output

        # If redshift should be fitted, we start with getting samples
        if z_source == "Fitted":
            if not hasattr(self.model, "predict_redshift_distribution"):
                self.logger.warn(
                    "Redshift fitting is not supported in that version of parsnip"
                )
                return t2_output
            z, z_probabilities = self.model.predict_redshift_distribution(
                sncosmo_table, max_redshift=self.max_fit_z
            )
        assert z is not None

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

        t2_output["predictions"] = predictions
        t2_output["classifications"] = classifcations

        if np.sum(z_probabilities) > 0:
            # Take redshift probabilities into account, if available
            if z_weights is not None:
                z_probabilities *= z_weights
            integrated_probabilities = z_probabilities.dot(probabilities)
            integrated_probabilities /= np.sum(integrated_probabilities)
            t2_output["marginal_lc_classifications"] = dict(
                zip(types, integrated_probabilities, strict=False)
            )
            # Find the best redshifts
            t2_output["z_at_minchi"] = z[np.argmax(z_probabilities)]
            # Map these to singular value predictions/lc_classifications
            # (wastes DB space, but possible to filter based on)
            t2_output["prediction"] = predictions[str(t2_output["z_at_minchi"])]
            t2_output["classification"] = classifcations[str(t2_output["z_at_minchi"])]
        else:
            # Not enough data for a chi2 estimate
            t2_output["Failed"] = "NoDOF"
            return t2_output

        # Plot
        if self.plot_suffix and self.plot_dir:
            # How to construct name?
            tname = compound.get("stock")
            # Need plotting tools to define id mapper
            # tname = ZTFIdMapper.to_ext_id(light_curve.stock_id)

            fig = plt.figure()
            ax = plt.gca()

            # Set redshift to best value and plot this fit
            lc_dataset.light_curves[0].meta["redshift"] = t2_output["z_at_minchi"]

            parsnip.plot_light_curve(lc_dataset.light_curves[0], self.model, ax=ax)
            plt.tight_layout()
            plt.savefig(
                os.path.join(self.plot_dir, f"t2parsnip_{tname}.{self.plot_suffix}")
            )

            plt.close("fig")
            plt.close("all")
            del fig
            gc.collect()

        return t2_output
