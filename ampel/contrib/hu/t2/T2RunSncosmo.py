#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2RunSncosmo.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                11.05.2021
# Last Modified Date:  18.02.2022
# Last Modified By:    simeon.reusch@desy.de


import copy
import errno
import os
from collections.abc import Sequence
from typing import Literal

import backoff
import matplotlib.pyplot as plt
import numpy as np
import sncosmo  # type: ignore[import]
from astropy.table import Table
from sncosmo.fitting import DataQualityError

# from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
# from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2BaseLightcurveFitter import T2BaseLightcurveFitter
from ampel.model.PlotProperties import PlotProperties
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView
from ampel.ztf.util.ZTFNoisifiedIdMapper import ZTFNoisifiedIdMapper


class T2RunSncosmo(T2BaseLightcurveFitter):
    """
    Gathers information and runs Sncosmo. Steps include:
    - Obtain model (read from file unless not in sncosmo registry)
    - Parse previous (chained) T2results for redshift or phase limits.
    - Creates suitable photometry, using the converter provided and  phase limits.
    - Defines model appropritately, including fit ranges and fixed parameters.
    - Run fit, potentially iterative in case of non-convergence.
    - Plot output if requested

    Fits lightcurves using SNCOSMO (using SALT2 defaultwise)
    which is assumed to be chained to other T2units that provide redshift and
    fit limits.

    Redshift selection parameters can be found in T2DigestRedshifts.
    Phase selection and dereddening parameter can be found in T2BaseLightcurveFitter.


    TODO:
    - Add option for redoing fits with disturbed initial conditions to avoid local minima
    - Add option for masking data?
    - Directly limit to bands which cover redshift model (and exit if none exists)
    - Save error message to log/t2 output

    Plot props sample use:
    plot_props: Optional[PlotProperties] = PlotProperties(
        tags=["SALT", "SNCOSMO"],
        file_name={
            "format_str": "%s_%s_%s.svg",
            "arg_keys": ["stock", "model", "redshift_kind"]
        },
        title={
            "format_str": "%s %s %s",
            "arg_keys": ["stock", "model", "redshift_kind"]
        },
        fig_text={
            "format_str": "%s %s \nz-source %s \nchisq %.2f ndof %s",
            "arg_keys": ["stock", "model", "redshift_kind", "chisq", "ndof"]
        },
        width=10,
        height=6,
        id_mapper="ZTFIdMapper"
    )


    """

    # Parameters regulating model
    # This unit requires that the model either exists in the current registry, or can be retrieved.
    # Non-standard models thus first have to be stored into the registry as part of the AMPEL init.
    sncosmo_model_name: str

    # Sncosmo parameters
    # Bounds - This is propagated directly into sncosmo. Beware e.g. clashed with catalog redshifts
    # When fitting redshift this needs to be included here, e.g. "sncosmo_bounds": {"z":(0.001,0.3)}
    sncosmo_bounds: dict[str, list[float]] = {}

    noisified: bool = False

    # Plot parameters
    plot_db: bool = False
    plot_props: None | PlotProperties = None  # Plot properties for SvgRecord creation
    plot_suffix: None | str = None  # Suffix if stored (locally) through matplotlib (e.g. _crayzmodel.png). Will add transient name
    plot_dir: str = "."  # Suffix if stored (locally) through matplotlib (e.g. _crayzmodel.png). Will add transient name

    # These are the units through which we look for redshifts
    # Which units should this be changed to
    t2_dependency: Sequence[  # type: ignore[assignment]
        StateT2Dependency[
            Literal[
                "T2CatalogMatch",
                "T2LSPhotoZTap",
                "T2CatalogMatchLocal",
                "T2MatchBTS",
                "T2ElasticcRedshiftSampler",
            ]
        ]
    ]

    def post_init(self) -> None:
        """
        Retrieve models.
        """

        # Setup model. Currently only parsnip treated separatly
        # If possible, use T2RunParnsip as the parsnip
        # sncosmo model is very slow.
        if self.sncosmo_model_name == "parsnip_plasticc":
            import parsnip  # type: ignore[import]

            source = parsnip.ParsnipSncosmoSource("plasticc")
        else:
            source = self.sncosmo_model_name  # Directly loaded

        self.sncosmo_model = sncosmo.Model(source=source)
        self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)

        # If redshift _should_ be provided we remove this from fit parameters
        if self.redshift_kind is not None or self.fixed_z is not None:
            self.fit_params.remove("z")

        self.default_param_vals = self.sncosmo_model.parameters

        # retry on with exponential backoff on "too many open files"
        self.process = backoff.on_exception(  # type: ignore[assignment]
            backoff.expo,
            OSError,
            giveup=lambda exc: not isinstance(exc, OSError)
            or exc.errno != errno.EMFILE,
            logger=self.logger,  # type: ignore[arg-type]
            max_time=300,
        )(self.process)

        # Load potential MW extinction map
        super().post_init()

    def _get_fit_metrics(self, sncosmo_result, sncosmo_table, sncosmo_model) -> dict:
        """
        Obtain metrics such as peak magnitude and lc residuals.
        Assumes that all models have at least z and t0 parameters.
        """

        # Fixing method parameters here, to avoid overloading unit params.
        detection_sigma = (
            3  # Detection sigma threshold to look for phase of first detection
        )
        pull_range = [-10, 20]  # Phase range used when calculating uniform chi2/dof

        z = sncosmo_model.get("z")

        lc_metrics = {}
        lc_metrics["restpeak_model_absmag_B"] = sncosmo_model.source_peakabsmag(
            "bessellb", "ab"
        )
        # Assuming all models have t0 as peak time parameter
        try:
            lc_metrics["obspeak_model_B"] = sncosmo_model.bandmag(
                "bessellb", "ab", sncosmo_model.get("t0")
            )
        except ValueError:
            # Likely too high redshift for predicting mag
            lc_metrics["obspeak_model_B"] = None

        sncosmo_table["phase"] = (sncosmo_table["time"] - sncosmo_model.get("t0")) / (
            1 + z
        )
        # Determine the phase of the first 3 sigma detection
        i_first = np.where(
            (sncosmo_table["flux"] / sncosmo_table["fluxerr"]) > detection_sigma
        )[0]
        # table might not be ordered
        lc_metrics[f"phase_{detection_sigma}sigma"] = np.min(
            sncosmo_table["phase"][i_first]
        )

        # Determine the explosion time (JD) according to the model
        # i.e. first time when model was defined.
        lc_metrics["jd_model_start"] = (
            sncosmo_model.source.minphase() + sncosmo_model.get("t0")
        )

        # Determine the chi/dof and dof for observations around peak light
        pulls = []
        for band in np.unique(sncosmo_table["band"]):
            band_tab = sncosmo_table[
                (sncosmo_table["band"] == band)
                & (sncosmo_table["phase"] >= pull_range[0])
                & (sncosmo_table["phase"] <= pull_range[1])
            ]
            # band_pulls = (band_tab["flux"] - sncosmo_model.bandflux(
            #    band, band_tab["jd"], zp=25., zpsys='ab'))
            try:
                # Using the same zeropoint / sys as when creating the table above
                band_pulls = (
                    band_tab["flux"]
                    - sncosmo_model.bandflux(
                        band, band_tab["time"], zp=25.0, zpsys="ab"
                    )
                ) / band_tab["fluxerr"]
                pulls.extend(list(band_pulls))
            except ValueError:
                self.logger.info("Sncosmo get fit metric error")
                lc_metrics["pull_retrieval_error"] = True

        lc_metrics["nbr_peak_pulls"] = len(pulls)
        lc_metrics["absmean_peak_pull"] = np.mean(np.abs(pulls))

        return lc_metrics

    def fit_sncosmo(
        self,
        phottab: Table,
        z: float | None,
    ) -> dict:
        """
        Fit the loaded sncosmo model at the redshift provided.

        Return collected info:

        """

        # Fitting section
        # To handle multiple redshifts, fitting section below should be put into function.
        # Do we keep the best fit, or all fits (as in t2parsnip)

        self.sncosmo_model.parameters = (
            self.default_param_vals
        )  # Reset models (if previoius runs done)
        if z is not None:
            self.sncosmo_model.set(z=z)

        # Carry out fit. Bounds are directly carried from parameters
        # todo: gravefully check which observed bands cover redshifted model
        try:
            sncosmo_result, fitted_model = sncosmo.fit_lc(
                phottab,
                self.sncosmo_model,
                self.fit_params,
                bounds=self.sncosmo_bounds,
            )
        except ValueError:
            self.logger.info("Sncosmo fit error")
            return {"run_error": True}
        except RuntimeError:
            # Might have worked with different initial conditions?
            self.logger.info("Sncosmo fit error")
            return {"run_error": True}
        except DataQualityError:
            self.logger.info("Sncosmo fit error")
            return {"run_error": True}

        self.logger.debug(f"Run results {sncosmo_result}")

        # Derive model metrics
        sncosmo_result["fit_metrics"] = self._get_fit_metrics(
            sncosmo_result, phottab, fitted_model
        )

        # How to best serialize these for mongo storage?
        sncosmo_result["parameters"] = sncosmo_result["parameters"].tolist()
        sncosmo_result["data_mask"] = sncosmo_result["data_mask"].tolist()
        try:
            sncosmo_result["covariance"] = sncosmo_result["covariance"].tolist()
        except KeyError:
            sncosmo_result["covariance"] = []

        # For filtering purposes we want a proper dict
        sncosmo_result["paramdict"] = {}
        for ix, pname in enumerate(sncosmo_result["param_names"]):
            sncosmo_result["paramdict"][pname] = sncosmo_result["parameters"][ix]

        # Return also model
        return {"fitresult": sncosmo_result, "fitted_model": fitted_model}

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
        T2LoadRedshift (redshift from .csv file)
        T2PhaseLimit (fit time-limits as determined from lightcurve)

        Returns
        -------
        dict
        """

        # Initialize output dict
        t2_output: dict[str, UBson] = {"model_name": self.sncosmo_model_name}

        # Load photometry as table
        # fitdatainfo contains redshift, if requested
        (sncosmo_table, fitdatainfo) = self.get_fitdata(datapoints, t2_views)
        t2_output.update(fitdatainfo)
        if sncosmo_table is None:
            return t2_output

        # Ok, sncosmo cannot handle negative fluxes?

        # Fitting section
        if fitdatainfo["z_source"] == "Fitted":
            # Run fit with z as free parameter.
            # Add this...
            fitout = self.fit_sncosmo(sncosmo_table, None)
            fitted_model = fitout["fitted_model"]
            t2_output["sncosmo_result"] = fitout["fitresult"]
        else:
            fit_allz = {
                str(z): self.fit_sncosmo(sncosmo_table, z) for z in fitdatainfo["z"]
            }

            # Get z for best results (min chi2)
            chisqs: dict[str, float] = {
                z: float(fitatz["fitresult"]["chisq"] / fitatz["fitresult"]["ndof"])
                for z, fitatz in fit_allz.items()
                if "fitresult" in fitatz and fitatz["fitresult"]["success"]
            }
            if len(chisqs) > 0:
                # At least one good fit
                bestz = min(chisqs, key=chisqs.get)  # type: ignore[arg-type]
                t2_output["sncosmo_result"] = fit_allz[bestz]["fitresult"]
                t2_output["best_z"] = bestz
                fitted_model = fit_allz[bestz]["fitted_model"]
                t2_output["result_at_z"] = {
                    z: fitatz["fitresult"]
                    for z, fitatz in fit_allz.items()
                    if "fitresult" in fitatz and fitatz["fitresult"]["success"]
                }
            else:
                fitted_model = None
                t2_output["success"] = False
                t2_output["sncosmo_result"] = {"chisq": -1, "ndof": -1}

        chisq: float = -1.0
        ndof: int = 0
        if isinstance(t2_output["sncosmo_result"], dict):
            chisq = float(t2_output["sncosmo_result"]["chisq"])
            ndof = t2_output["sncosmo_result"]["ndof"]

        # Save plot
        if fitted_model and chisq > 0 and (self.plot_props or self.plot_suffix):
            # Construct name JN: What are the standards for noisified?
            stock_id = "-".join([str(v) for v in self.get_stock_id(datapoints)])
            tname = "-".join([str(v) for v in self.get_stock_name(datapoints)])

            if self.noisified:
                stock_id = "-".join([str(v) for v in self.get_stock_id(datapoints)])
                tname = ZTFNoisifiedIdMapper().to_ext_id(stock_id)

            # Add some info
            plot_fig_text = "{} {} {} \nchisq {:.2f}\nndof {}".format(
                tname,
                self.sncosmo_model_name,
                self.redshift_kind,
                chisq,
                ndof,
            )
            plot_extra = {
                "model": self.sncosmo_model_name,
                "redshift_kind": self.redshift_kind,
                "chisq": chisq,
                "ndof": ndof,
                "stock": self.get_stock_id(datapoints)[
                    0
                ],  # Only using first name, assuming this is from ZTF
            }

            fig = sncosmo.plot_lc(
                sncosmo_table,
                model=fitted_model,
                pulls=True,
                figtext=plot_fig_text,
                ncol=3,
                # sncosmo using old syntax, crashing if cmap not given
                cmap=plt.cm.jet_r,  #  type: ignore[attr-defined]
                # fill_data_marker = self.fit_mask, # Activate if add fit mask
            )

            if self.plot_props:
                plots = [
                    create_plot_record(
                        fig, self.plot_props, plot_extra, logger=self.logger
                    )
                ]
                # Also store to DB if requested
                if self.plot_db:
                    t2_output["plots"] = plots
            if self.plot_suffix:
                fpath = os.path.join(self.plot_dir, tname + self.plot_suffix)
                fig.savefig(fpath)

        return t2_output
