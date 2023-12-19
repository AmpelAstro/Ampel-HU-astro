#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2RunSncosmo.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                11.05.2021
# Last Modified Date:  18.02.2022
# Last Modified By:    simeon.reusch@desy.de


import os
import copy
import errno
from collections.abc import Sequence
from typing import Literal, Optional

import backoff
import numpy as np
import sncosmo  # type: ignore[import]
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.model.PlotProperties import PlotProperties
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from ampel.ztf.util.ZTFNoisifiedIdMapper import ZTFNoisifiedIdMapper
from astropy.table import Table
from sfdmap2.sfdmap import SFDMap  # type: ignore[import]
from sncosmo.fitting import DataQualityError


class T2RunSncosmo(AbsTiedStateT2Unit, AbsTabulatedT2Unit):
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

    # Redshift usage options. Current options
    # T2MatchBTS : Use the redshift published by BTS and  synced by that T2.
    # T2DigestRedshifts : Use the best redshift as parsed by DigestRedshift.
    # None : run sncosmo template fit with redshift as free parameter OR use backup_z if set
    redshift_kind: None | str
    # If loading redshifts from DigestRedshifts, provide the max ampel z group to make use of.
    # (note that filtering based on this can also be done for a potential t3)
    max_ampelz_group: int = 3
    # If none of the above is selected, a fixed redshift can be provided.
    # If this is None, the redshift will be included as a free parameter in the fit
    fixed_z: None | float
    # It is also possible to use fixed redshift whenever a dynamic redshift kind is not possible
    backup_z: None | float
    # Finally, the provided lens redshift might be multiplied with a scale
    # Useful for lensing studies, or when trying multiple values
    scale_z: None | float

    # Sncosmo parameters
    # Bounds - This is propagated directly into sncosmo. Beware e.g. clashed with catalog redshifts
    # When fitting redshift this needs to be included here, e.g. "sncosmo_bounds": {"z":(0.001,0.3)}
    sncosmo_bounds: dict[str, list[float]] = {}
    # Remove MW dust absorption using SFD maps. This assumes that the position
    # can be retrieved from the light_curve and that the SFD_DIR env var is set
    # to allow them to be found. The default value of Rv will be used.
    apply_mwcorrection: bool = False

    # Phase range usage. Current option:
    # T2PhaseLimit : use the jdmin jdmax provided in this unit output
    # None : use full datapoint range
    # (T2BayesianBlocks should be added)
    phaseselect_kind: None | str

    noisified: bool = False


    # Plot parameters
    plot_db: bool = False
    plot_props: None | PlotProperties = None   # Plot properties for SvgRecord creation
    plot_matplotlib_suffix: None | str = None    # Suffix if stored (locally) through matplotlib (e.g. _crayzmodel.png). Will add transient name 
    plot_matplotlib_dir: str = '.'    # Suffix if stored (locally) through matplotlib (e.g. _crayzmodel.png). Will add transient name 


    # Units from which time limits to use or redshifts can be picked. 
    t2_dependency: Sequence[StateT2Dependency[Literal[
        "T2ElasticcRedshiftSampler",
        "T2DigestRedshifts",
        "T2MatchBTS",
        "T2PhaseLimit",
        "T2XgbClassifier"]]]



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

        # Setup model, with or without MW correction
        if self.apply_mwcorrection:
            dust = sncosmo.models.CCM89Dust()
            self.sncosmo_model = sncosmo.Model(
                source=source,
                effects=[dust],
                effect_names=["mw"],
                effect_frames=["obs"],
            )
            self.dustmap = SFDMap()
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)
            self.fit_params.remove("mwebv")
            self.fit_params.remove("mwr_v")
        else:
            self.sncosmo_model = sncosmo.Model(source=source)
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)

        # If redshift _should_ be provided we remove this from fit parameters
        if self.redshift_kind is not None or self.backup_z is not None:
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

    # Should be rearranged to provide a list of redshifts a la T2RunParsnip with signature as this 
#    def _get_redshift(self, t2_views) -> tuple[Optional[list[float]], Optional[str], Optional[list[float]]]:
    def _get_redshift(self, t2_views) -> tuple[None | float, None | str]:
        """
        Can potentially also be replaced with some sort of T2DigestRershift tabulator?

        Assuming that only one instance of redshift sources exist
        """

        # Examine T2s for eventual information
        z: Optional[list[float]] = None
        z_source: Optional[str] = None
        z_weights: Optional[list[float]] = None


        if self.redshift_kind in ['T2MatchBTS', 'T2DigestRedshifts', 'T2ElasticcRedshiftSampler']:
            for t2_view in t2_views:
                if not t2_view.unit == self.redshift_kind:
                    continue
                self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                # Parse this
                if self.redshift_kind == 'T2MatchBTS':
                    if 'bts_redshift' in t2_res.keys() and not t2_res['bts_redshift'] == '-':
                        z = [float(t2_res['bts_redshift'])]
                        z_source = "BTS"
                elif self.redshift_kind == 'T2DigestRedshifts':
                    if ('ampel_z' in t2_res.keys() and t2_res['ampel_z'] is not None
                            and t2_res['group_z_nbr'] <= self.max_ampelz_group):
                        z = [float(t2_res['ampel_z'])]
                        z_source = "AMPELz_group" + str(t2_res['group_z_nbr'])
                elif self.redshift_kind == 'T2ElasticcRedshiftSampler':
                    z = t2_res['z_samples']
                    z_source = t2_res['z_source']
                    z_weights = t2_res['z_weights']
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
            z = [onez*self.scale_z for onez in z]
            z_source += " + scaled {}".format(self.scale_z)

        # TODO: return the list instead of this
        # return z, z_source, z_weights
        # We now simply pick the middle number
        if isinstance(z_weights, list):
            z = z[ z_weights.index( max(z_weights) ) ]  # type: ignore
            print('INPUT')
            print(z, z_weights)
            print('SNCOSMOS z', z)
        elif isinstance(z, list):
            if len(z) % 2 != 0:
                z = z[int(len(z) / 2)]  # type: ignore
            else:
                z = (((z[int(len(z) / 2)]) + (z[int(len(z) / 2) - 1])) / 2)  # type: ignore
                
        return z, z_source  # type: ignore


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
                self.logger.debug("Parsing t2 results from {}".format(t2_view.unit))
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )
                jdstart = t2_res["t_start"]
                jdend = t2_res["t_end"]

        return jdstart, jdend

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
        lc_metrics["phase_{}sigma".format(detection_sigma)] = np.min(
            sncosmo_table["phase"][i_first]
        )

        # Determine the explosion time (JD) according to the model
        # i.e. first time when model was defined.
        lc_metrics[
            "jd_model_start"
        ] = sncosmo_model.source.minphase() + sncosmo_model.get("t0")

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
            except ValueError as e:
                self.logger.info("Sncosmo get fit metric error")
                lc_metrics["pull_retrieval_error"] = True

        lc_metrics["nbr_peak_pulls"] = len(pulls)
        lc_metrics["absmean_peak_pull"] = np.mean(np.abs(pulls))

        return lc_metrics

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
        t2_output: dict[str, UBson] = {"model_name": self.sncosmo_model_name}

        # Obtain redshift
        z, z_source = self._get_redshift(t2_views)
        t2_output["z"] = z
        t2_output["z_source"] = z_source
        # A source class of None indicates that a redshift source was required, but not found.
        if t2_output["z_source"] is None:
            return t2_output

        # Check for phase limits
        (jdstart, jdend) = self._get_phaselimit(t2_views)
        t2_output["jdstart"] = jdstart
        t2_output["jdend"] = jdend
        if t2_output["jdstart"] is None:
            return t2_output

        # Obtain photometric table
        sncosmo_table = self.get_flux_table(datapoints)
        print("T2RUNSNCOSMO:: ", sncosmo_table)
        sncosmo_table = sncosmo_table[
            (sncosmo_table["time"] >= jdstart) & (sncosmo_table["time"] <= jdend)
        ]

        self.logger.debug("Sncosmo table {}".format(sncosmo_table))

        # Fitting section
        # To handle multiple redshifts, fitting section below should be put into function.
        # Do we keep the best fit, or all fits (as in t2parsnip)
                
        self.sncosmo_model.parameters = self.default_param_vals  # Reset

        # Define fit parameter and ranges
        if self.apply_mwcorrection:
            transient_mwebv = self.dustmap.ebv(*self.get_pos(datapoints, which="mean"))
            self.sncosmo_model.set(mwebv=transient_mwebv)

        # Set redshift if provided
        if isinstance(t2_output["z"], float):
            self.sncosmo_model.set(z=t2_output["z"])

        self.logger.debug(
            "Starting fit with fit params {}, all parameters {} and start values {}".format(
                self.fit_params,
                self.sncosmo_model.param_names,
                self.sncosmo_model.parameters,
            )
        )

        # Carry out fit. Bounds are directly carried from parameters
        # todo: gravefully check which observed bands cover redshifted model
        try:
            sncosmo_result, fitted_model = sncosmo.fit_lc(
                sncosmo_table,
                self.sncosmo_model,
                self.fit_params,
                bounds=self.sncosmo_bounds,
            )
        except ValueError as e:
            self.logger.info("Sncosmo fit error")
            print("value error", e)
            t2_output["run_error"] = True
            return t2_output
        except RuntimeError as e:
            # Might have worked with different initial conditions?
            print("value error", e)
            self.logger.info("Sncosmo fit error")
            t2_output["run_error"] = True
            return t2_output
        except DataQualityError as e:
            print("value error", e)
            self.logger.info("Sncosmo fit error")
            t2_output["run_error"] = True
            return t2_output

        self.logger.debug("Run results {}".format(sncosmo_result))

        # Derive model metrics
        t2_output["fit_metrics"] = self._get_fit_metrics(
            sncosmo_result, sncosmo_table, fitted_model
        )

        # How to best serialize these for mongo storage?
        sncosmo_result["parameters"] = sncosmo_result["parameters"].tolist()
        sncosmo_result["data_mask"] = sncosmo_result["data_mask"].tolist()
        try:
            sncosmo_result["covariance"] = sncosmo_result["covariance"].tolist()
        except:
            sncosmo_result["covariance"] = []

        # For filtering purposes we want a proper dict
        sncosmo_result["paramdict"] = {}
        for ix, pname in enumerate(sncosmo_result["param_names"]):
            sncosmo_result["paramdict"][pname] = sncosmo_result["parameters"][ix]

        # Finish up
        t2_output["sncosmo_result"] = sncosmo_result

        # Save plot
        if self.plot_props or self.plot_matplotlib_suffix:
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
                sncosmo_result["chisq"],
                sncosmo_result["ndof"],
            )
            plot_extra = {
                "model": self.sncosmo_model_name,
                "redshift_kind": self.redshift_kind,
                "chisq": sncosmo_result["chisq"],
                "ndof": sncosmo_result["ndof"],
                "stock": self.get_stock_id(datapoints)[0]   # Only using first name, assuming this is from ZTF
            }
            
            fig = sncosmo.plot_lc(
                sncosmo_table,
                model=fitted_model,
                pulls=True,
                figtext=plot_fig_text,
                ncol=3,
                # fill_data_marker = self.fit_mask, # Activate if add fit mask
            )


            if self.plot_props:
                plots = [
                    create_plot_record(fig, self.plot_props, plot_extra, logger=self.logger)
                ]
                # Also store to DB if requested
                if self.plot_db:
                    t2_output["plots"] = plots
            if self.plot_matplotlib_suffix:
                fpath = os.path.join(self.plot_matplotlib_dir,tname+self.plot_matplotlib_suffix)
                fig.savefig(fpath)

        return t2_output
