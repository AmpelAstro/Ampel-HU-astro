#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2TabulatorRiseDecline.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                28.12.2018
# Last Modified Date:  03.08.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import os
import re
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, Any

import light_curve
import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit

from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.base.LogicalUnit import LogicalUnit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.protocol.LoggerProtocol import LoggerProtocol


def getMag(tab: Table, err=False, time=False):
    """
    Shorthand to getting the magnitude from flux.
    """

    iPos = tab["flux"] > 0
    m = -2.5 * np.log10(tab["flux"][iPos]) + tab["zp"][iPos]
    if err:
        # Simple symmetric method
        merr = 2.5 / np.log(10) * (tab["fluxerr"][iPos] / tab["flux"][iPos])
        if time:
            return m, merr, tab["time"][iPos]
        return m, merr
    if time:
        return m, tab["time"][iPos]
    return m


def getMeanflux(tab: Table, jdstart, jdend):
    """
    Return dict with the mean flux for all bands
    with data in the selected range.
    """

    means = {}
    used_bands = []
    tt = tab[(tab["time"] >= jdstart) & (tab["time"] <= jdend)]
    for band in set(tt["band"]):
        btt = tt[tt["band"] == band]
        means[band] = np.average(btt["flux"], weights=1.0 / btt["fluxerr"] ** 2)
        #        means[band] = 1
        used_bands.append(band)
    return means, used_bands


def getBandBits(bands: Sequence):
    """
    Return number quantifying which bands are included.
    """
    bandval = {
        "lsstu": 1,
        "lsstg": 2,
        "lsstr": 4,
        "lssti": 8,
        "lsstz": 16,
        "lssty": 32,
        "ztfg": 64,
        "ztfr": 128,
        "ztfi": 256,
    }
    index = 0
    for band in bands:
        index += bandval[band]
    return index


class FitFailed(RuntimeError):
    ...


class T2TabulatorRiseDeclineBase(AmpelBaseModel):
    """
     Derive a number of simple metrics describing the rise, peak and decline of a lc.

     This version assumes input provided by flux table.

     Derived values:
     * t_predetect : time between first detection and previous non-detection
                     in "significant bands". Todo: add min error?
     * t_lc : duration (time between first and most recent detection)
     * jd_max : jd of peak light. None unless bool_peaked
     * jd_det : jd of first detection
     * jd_last : jd of last detection
     * ndet : number of significant detections
     * bool_peaked : is the lc estimated to be declining?
     * bool_pure : has there been no significant non-detections after first detection?
     * bool_rise : was the peak light within "cadence" days of the most recent detection?
    # Not yet implemented    * bool_norise : was the first detection NOT significantly fainter than
         mag_peak IF bool_peaked, ELSE mag_lst
     * bool_hasgaps : The lc has a gap between detections of at least 30 days,
         indicating either a recurrent event or a chance coincidental detection.
     * mag_peak : magnitude at peak light (significant band). Only calculated if bool_peaked
     * mag_det : detection magnitude (significant band)
     * mag_last : magnitude of last detection (significant band)
     Following done for each "color_list" i-ii
     * i-ii_peak : color at peak. None unless bool_peaked AND i+ii obs made within
                  "cadence"  days of jd_max
     * i-ii_det : color at detection. None unless i+ii obs made within "cadence"
                   days of jd_det
     * i-ii_last : color at last detection. None unless i+ii obs made within
                   "cadence" days of jd_last
     * slope_rise_{i,ii} : magnitude slope between jd_det and jd_max. None if bool_norise
     * slope_decline_{i,ii} : magnitude slope between jd_max and jd_lst. None unless bool_peaked
     * postpeak_fluxevo_{band} : flux ratio between peak and {dt_fluxevo}

     Additionally, we would like to access host properties like distance and host mags and sizes.
     Would that have to be a different T2?.
     Rather make a different T2 which is chained to the redshift sampler and this.

     "t_cadence" is meant to approximate the rough practical cadence, i.e. how often repeated
     observations should come and what can be seen as "simulateneous"

     "significant_bands" is a list of "deep bands". Collected data from these will be
     used to calculate ages, peak magnitude etc.

     "sigma_det" quantifies what should be seen as a detection.

     "color_list" contains a list of lists of colors which should be considered

    """

    t_cadence: float = 5.0
    significant_bands: Sequence[str] = ["lsstg", "lsstr", "lssti", "lsstz"]
    sigma_det: float = 5.0
    sigma_slope: float = 3.0  # Threshold for having detected a slope (flux/day)
    dt_fluxevo: float = 20.0  # Calculating flux ratio between peak and this time (roughly modeled after SNIa 2nd peak)
    color_list: Sequence[Sequence[str]] = [
        ["lsstu", "lsstg"],
        ["lsstg", "lsstr"],
        ["lsstr", "lssti"],
        ["lssti", "lsstz"],
        ["lsstz", "lssty"],
        ["ztfg", "ztfr"],
        ["ztfr", "ztfi"],
    ]
    max_tgap: int = 30

    # Cut the lightcurve if longer than this limit.
    # Motivated by worse classification for longer (inbalanced training?)
    # max_ndet: int = 20
    # For new training round, change this
    max_ndet: int = 200000

    if TYPE_CHECKING:
        logger: LoggerProtocol

    def get_bandfeatures(self, ftable):
        """
        Go through all bands of input table and:
        - Try to find a peak, together with rise and fall slopes.
        - Based on this try to determine if it is rising, falling and have peaked.
        - Storing slopes where possible.
        """

        def linearFunc(x, intercept, slope):
            return intercept + slope * x

        banddata = {}
        tscale = np.mean(ftable["time"])

        for band in set(ftable["band"]):
            bt = ftable[ftable["band"] == band]

            max_flux = bt["flux"].max()
            max_flux_time = bt[bt["flux"] == max_flux]["time"][0]
            banddata["jd_peak_" + band] = max_flux_time

            # Divide Table
            riset = bt[bt["time"] <= max_flux_time]
            fallt = bt[bt["time"] >= max_flux_time]

            # Examine rise
            if len(riset) > 1:
                try:
                    fit, cov = curve_fit(
                        linearFunc,
                        riset["time"] - tscale,
                        riset["flux"],
                        sigma=riset["fluxerr"],
                        absolute_sigma=True,
                    )
                    banddata["rise_slope_" + band] = fit[1]
                    banddata["rise_slopesig_" + band] = fit[1] / np.sqrt(cov[1][1])
                except RuntimeError:
                    self.logger.info("Risetime curve fit failed.")
            if len(fallt) > 1:
                try:
                    fit, cov = curve_fit(
                        linearFunc,
                        fallt["time"] - tscale,
                        fallt["flux"],
                        sigma=fallt["fluxerr"],
                        absolute_sigma=True,
                    )
                    banddata["fall_slope_" + band] = fit[1]
                    banddata["fall_slopesig_" + band] = fit[1] / np.sqrt(cov[1][1])
                except RuntimeError:
                    self.logger.info("Falltime curve fit failed.")
                    
            # Check for flux decline until the time for a possible second bump 
            if sum( 
                    (
                    isecpeak:=(np.abs(fallt["time"] - max_flux_time - self.dt_fluxevo) < self.t_cadence/2)
                    ) 
                )>0:
                banddata["fluxevo_ratio_" + band] = np.mean(fallt["flux"][isecpeak]) / max_flux




        # Check whether we have a significant rise detected in any band.
        risepulls = [
            banddata.get("rise_slopesig_" + band, 0) for band in set(ftable["band"])
        ]
        if sum(risepulls) > self.sigma_slope:
            banddata["bool_rise"] = True
        else:
            banddata["bool_rise"] = False
        # Check whether we see a decline
        decpulls = [
            banddata.get("fall_slopesig_" + band, 0) for band in set(ftable["band"])
        ]
        if sum(decpulls) < -self.sigma_slope:
            banddata["bool_fall"] = True
        else:
            banddata["bool_fall"] = False

        # If the transient has both a rise and a fall we can
        # define a central peak
        if (
            banddata["bool_rise"]
            and banddata["bool_fall"]
            and (
                len(
                    peakjds := [
                        banddata["jd_peak_" + band]
                        for band in set(ftable["band"])
                        if "rise_slope_" + band in banddata
                        and "fall_slope_" + band in banddata
                    ]
                )
                > 0
            )
        ):
            banddata["bool_peaked"] = True
            # Use jd of all bands for which we could estimate rise+fall
            banddata["jd_peak"] = np.median(peakjds)
        else:
            banddata["bool_peaked"] = False

        # Could include the abs mag at peak, but we argued this would not
        # be as useful?

        return banddata

    def cut_flux_table(self, flux_table: Table) -> Table:
        """
        Limit to some _total_ set of significant detections.
        """
        sig_mask = np.abs((flux_table["flux"]) / flux_table["fluxerr"]) > self.sigma_det
        sig_time = list(flux_table["time"][sig_mask])
        if len(sig_time) > self.max_ndet:
            max_time = sorted(sig_time)[self.max_ndet]
            flux_table = flux_table[flux_table["time"] <= max_time]

        return flux_table

    def average_filtervalues(
        self, features: dict[str, Any], matchlist: Iterable[str] = tuple()
    ) -> dict[str, Any]:
        """
        Many (most) features calculated per band or color, with available filter/color varying from object
        to object.

        This method will attempt to average over any filter or color-specific fields available and create an averaged entry.

        A base matchlist is created based on significant bands. (plus potentially an input matchlist)

        """

        mymatchlist = [f"_{band}_" for band in self.significant_bands]
        # Ignore the colors for now, would need to rework how these averages are described.
        # mymatchlist.extend(
        # ['{}-{}_'.format(col[0],col[1]) for col in self.color_list]
        #    )
        mymatchlist.extend(matchlist)

        matchedvalues: dict[str, list] = {}
        # Can this be redone to avoid nested loop?
        for match in mymatchlist:
            p = re.compile(match)
            for key, val in features.items():
                stub = p.sub("_", key)  # Name of averaged entry in return dict
                if stub == key:  # Match string not found
                    continue
                if stub not in matchedvalues:
                    matchedvalues[stub] = []
                matchedvalues[stub].append(val)

        return {k: np.nanmean(v) for k, v in matchedvalues.items()}

    def compute_stats(self, flux_table: Table) -> dict[str, Any]:
        # Output dict that we will start to populate
        o: dict[str, Any] = {}

        # Create subset of table with significant detections in significant bands
        band_mask = [
            bandobs in self.significant_bands for bandobs in flux_table["band"]
        ]
        sig_mask = np.abs((flux_table["flux"]) / flux_table["fluxerr"]) > self.sigma_det
            
        det_table = flux_table[band_mask & sig_mask]
        # Calculate fraction negative detection (we no longer cut only because of this)
        o["ndet"] = len(det_table)

        if o["ndet"] == 0:
            o["success"] = False
            o["cause"] = "No data survive significance criteria."
            # Gather additional information for evaluation
            o["alldet"] = len(flux_table)
            neg_mask = (-flux_table["flux"]) / flux_table["fluxerr"] > self.sigma_det
            o["nnegdet"] = len(flux_table[band_mask & neg_mask])

            return o

        o["frac_pos"] = np.sum(det_table["flux"] > 0) / o["ndet"]

        o["jd_det"] = det_table["time"].min()
        o["jd_last"] = det_table["time"].max()
        o["t_lc"] = o["jd_last"] - o["jd_det"]

        # Get the max time of obs in signifant bands prior to jd_det
        if flux_table[band_mask]["time"].min() < o["jd_det"]:
            o["t_predetect"] = (
                o["jd_det"]
                - flux_table[band_mask][flux_table["time"][band_mask] < o["jd_det"]][
                    "time"
                ].max()
            )
        else:
            o["t_predetect"] = None

        # Magnitude and color calculations below assume detections are positive
        if (
            flux_table[flux_table["time"] == o["jd_det"]]["flux"] < 0
            or flux_table[flux_table["time"] == o["jd_last"]]["flux"] < 0
        ):
            # Did we succeed with feature calculation or not?
            o["success"] = True
            return o

        o["mag_det"] = float(getMag(flux_table[flux_table["time"] == o["jd_det"]]))
        o["band_det_id"] = getBandBits(
            [flux_table[flux_table["time"] == o["jd_det"]]["band"][0]]
        )

        o["mag_last"] = float(getMag(flux_table[flux_table["time"] == o["jd_last"]]))
        o["band_last_id"] = getBandBits(
            [flux_table[flux_table["time"] == o["jd_last"]]["band"][0]]
        )

        # Check fails irregularly
        try:
            o["mag_min"] = float(
                getMag(flux_table[flux_table["flux"] == max(flux_table["flux"])])
            )
        except TypeError:
            self.logger.info("Mag min extraction failed.")

        try:
            o["jd_min"] = float(flux_table["time"][np.argmax(flux_table["flux"])])
        except TypeError:
            self.logger.info("JD at mag min extraction failed.")

        # Check for non-signficant obs between det and last
        ultab = flux_table[band_mask & ~sig_mask]
        if sum((ultab["time"] >= o["jd_det"]) & (ultab["time"] <= o["jd_last"])) == 0:
            # No non-detection among signifcant bands between start and end
            o["bool_pure"] = True
        else:
            o["bool_pure"] = False

        # We start measuring bandfeatures t_cadence days prior to first detection
        # to allow some nondetection data to be included.
        # (Note: we previously only used significant bands here. Prob wrong. )
        time_mask = (flux_table["time"] > (o["jd_det"] - self.t_cadence)) & (
            flux_table["time"] < (o["jd_last"] + self.t_cadence)
        )
        o.update(self.get_bandfeatures(flux_table[time_mask]))

        # If there is a peak we additionally check whether this is within t_cadence
        # days of detector or last, and call this fastrise and fastfall
        o["bool_fastrise"], o["bool_fastfall"] = None, None
        if o.get("bool_peaked", False):
            if np.abs(o["jd_det"] - o["jd_peak"]) < self.t_cadence:
                o["bool_fastrise"] = True
            else:
                o["bool_fastrise"] = False
            if np.abs(o["jd_last"] - o["jd_peak"]) < self.t_cadence:
                o["bool_fastfall"] = True
            else:
                o["bool_fastfall"] = False
            o["t_rise"] = o["jd_peak"] - o["jd_det"]
            o["t_fall"] = o["jd_last"] - o["jd_peak"]

        # Are there long gaps among the detections?
        jdsorted = np.unique(flux_table["time"])
        if len(jdsorted) > 1:
            if (jdsorted[1:] - jdsorted[0:-1]).max() > self.max_tgap:
                o["bool_hasgaps"] = True
            else:
                o["bool_hasgaps"] = False
        else:
            o["bool_hasgaps"] = None

        # Color
        # Define time subsets at detection, last (significant) and peak (if defined)
        # In each, get the mean flux in each band.
        # We assume that the zeropoint is the same for these fluxes!
        # Also, only exist for positive fluxes (...)

        # Detection colors
        fluxdict, fluxbands = getMeanflux(
            flux_table,
            o["jd_det"] - self.t_cadence / 2,
            o["jd_det"] + self.t_cadence / 2,
        )
        o["det_bands"] = getBandBits(fluxbands)
        for colbands in self.color_list:
            if fluxdict.get(colbands[0], -1) > 0 and fluxdict.get(colbands[1], -1) > 0:
                o[f"{colbands[0]}-{colbands[1]}_det"] = -2.5 * np.log10(
                    fluxdict[colbands[0]] / fluxdict[colbands[1]]
                )
        # Last obs colors
        fluxdict, fluxbands = getMeanflux(
            flux_table,
            o["jd_last"] - self.t_cadence / 2,
            o["jd_last"] + self.t_cadence / 2,
        )
        o["last_bands"] = getBandBits(fluxbands)
        for colbands in self.color_list:
            if fluxdict.get(colbands[0], -1) > 0 and fluxdict.get(colbands[1], -1) > 0:
                o[f"{colbands[0]}-{colbands[1]}_last"] = -2.5 * np.log10(
                    fluxdict[colbands[0]] / fluxdict[colbands[1]]
                )
        # Peak colors, if found
        if o.get("bool_peaked", False):
            fluxdict, fluxbands = getMeanflux(
                flux_table,
                o["jd_peak"] - self.t_cadence / 2,
                o["jd_peak"] + self.t_cadence / 2,
            )
            o["peak_bands"] = getBandBits(fluxbands)
            for colbands in self.color_list:
                if (
                    fluxdict.get(colbands[0], -1) > 0
                    and fluxdict.get(colbands[1], -1) > 0
                ):
                    o[f"{colbands[0]}-{colbands[1]}_peak"] = -2.5 * np.log10(
                        fluxdict[colbands[0]] / fluxdict[colbands[1]]
                    )

        o["success"] = True
        return o


class BaseLightCurveFeatures(LogicalUnit, AmpelBaseModel):
    """
    Calculate various features of the light curve using the light-curve
    package described in https://ui.adsabs.harvard.edu/abs/2021MNRAS.502.5147M%2F/abstract

    Lifted from T2LightCurveFeatures
    Installed as  python3 -mpip install light-curve

    This v2 contains feature selection and choice of flux as unit based on most common features to influence
    a binary tree classifier (compare rank_tabulator_features and previous version of this unit).
    Also decided to require 4 detection in a band for using.
    """

    #: Features to extract from the light curve.
    #: See: https://docs.rs/light-curve-feature/0.2.2/light_curve_feature/features/index.html
    lightcurve_features_flux: dict[str, None | dict[str, Any]] = {
        "Eta": None,  # 2
        "MaximumSlope": None,  # 2
        "Periodogram": {"peaks": 1},  # 4?
        "Skew": None,  # 3
        "StandardDeviation": None,  # 2
        "ExcessVariance": None,  # 2
        "LinearFit": None,  # 3
        "AndersonDarlingNormal": None,  # 4
        "Kurtosis": None,  # 4
        "StetsonK": None,  # 2
    }

    #: Bandpasses to use
    lightcurve_bands: dict[str, Any] = {
        "ztfg": "ztfg",
        "ztfr": "ztfr",
        "ztfi": "ztfi",
        "lsstu": "lsstu",
        "lsstg": "lsstg",
        "lsstr": "lsstr",
        "lssti": "lssti",
        "lsstz": "lsstz",
        "lssty": "lssty",
    }

    def init_lightcurve_extractor(self) -> None:
        self.fluxextractor = light_curve.Extractor(
            *(
                getattr(light_curve, k)(**(v or {}))
                for k, v in self.lightcurve_features_flux.items()
            )
        )

    def post_init(self) -> None:
        self.init_lightcurve_extractor()

    def extract_lightcurve_features(self, flux_table: Table) -> dict[str, float]:
        result = {}
        #
        for band, fid in self.lightcurve_bands.items():
            if (in_band := flux_table[flux_table["band"] == fid]) is None:
                continue

            # Conversion to mag not used if features determined in flux space
            # m, merr, mtime = getMag(in_band, err=True, time=True)

            # We wrap this in a try statement, since there are a few ways these can fail, typically with poor data in one or another aspect
            try:
                # Flux based, requires at least 4 detections
                if len(in_band["flux"]) >= 4:
                    lcout = {
                        f"{k}_{band}_flux": v
                        for k, v in zip(
                            self.fluxextractor.names,
                            self.fluxextractor(
                                in_band["time"], in_band["flux"], in_band["fluxerr"]
                            ),
                            strict=False,
                        )
                    }
                    result.update(lcout)
            except ValueError:
                self.logger.info("lightkurve extract fail")
                pass

        return result


class T2TabulatorRiseDecline(
    AbsStateT2Unit,
    AbsTabulatedT2Unit,
    T2TabulatorRiseDeclineBase,
    BaseLightCurveFeatures,
):
    plot_prob: float = 0.0
    path_testplot: str = "./plots/"

    #    def __init__(self, **kwargs):
    #        super().__init__(**kwargs)
    #    def super().post_init(self):
    #        super().__init__(**kwargs)

    def test_plot(self, name, table, t2result):
        """
        for debugging

        Create panel for each band, showing data + the significant dets.

        Mark times for detection, last det + peak (if set)
        Somehow indicate tilt (from peak pos?)
        Write summary of boolean conclusion.

        Save as stockid + ndet

        """

        import matplotlib.pyplot as plt

        bands = set(table["band"])

        fig, axs = plt.subplots(2, 3)

        for k, band in enumerate(bands):
            bt = table[table["band"] == band]
            ax = axs[int(k / 3), k % 3]

            ax.errorbar(bt["time"], bt["flux"], yerr=bt["fluxerr"], fmt="o", label=band)

            # Detection times
            if "jd_det" in t2result:
                ax.axvline(t2result["jd_det"], color="grey", linewidth=3, alpha=0.5)
            if "jd_last" in t2result:
                ax.axvline(t2result["jd_last"], color="grey", linewidth=3, alpha=0.5)
            if t2result["bool_peaked"]:
                ax.axvline(t2result["jd_peak"], color="red", linewidth=2, alpha=0.8)
            if "jd_peak_" + band in t2result:
                ax.axvline(t2result["jd_peak_" + band], color="green", linewidth=1)

            # We next wich to indicate the slopes we have measured
            if "rise_slope_" + band in t2result:
                xvals = np.array([-t2result["t_lc"] / 4, t2result["t_lc"] / 4])
                yvals = xvals * t2result["rise_slope_" + band] + np.mean(bt["flux"])
                xvals += t2result["jd_det"] + t2result["t_lc"] / 4
                col = "black"
                if np.abs(t2result["rise_slopesig_" + band]) > 3:
                    col = "red"
                ax.plot(xvals, yvals, color=col)
            if "fall_slope_" + band in t2result:
                xvals = np.array([-t2result["t_lc"] / 4, t2result["t_lc"] / 4])
                yvals = xvals * t2result["fall_slope_" + band] + np.mean(bt["flux"])
                xvals += t2result["jd_last"] - t2result["t_lc"] / 4
                col = "black"
                if np.abs(t2result["fall_slopesig_" + band]) > 3:
                    col = "red"
                ax.plot(xvals, yvals, color=col)

            #            ax.legend()
            ax.set_xlabel("MJD")
            ax.set_ylabel(band)

            # Create text string
            title = "ndet: %s " % (t2result["ndet"])

            for boolprop in [
                "peaked",
                "pure",
                "rise",
                "hasgaps",
                "fall",
                "fastrise",
                "fastfall",
            ]:
                if t2result[f"bool_{boolprop}"]:
                    title += f"{boolprop} "

            ax.set_title(title, {"fontsize": 8})

        # Store figure
        os.makedirs(self.path_testplot, exist_ok=True)
        path = os.path.join(
            self.path_testplot, "{}_{}.pdf".format(name, t2result["ndet"])
        )
        plt.tight_layout()
        plt.savefig(path)
        plt.clf()

    def process(
        self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
    ) -> dict[str, Any]:
        """
        Process datapoints belonging to one state of one transient.
        A commong Table is generated which is used as input
        to the feature generator.
        """

        # Convert input datapoints to standardized Astropy Table
        # Using standard tabulators
        flux_table = self.get_flux_table(datapoints)

        # Cut the flux table if requested
        if self.max_ndet > 0 and len(flux_table) > self.max_ndet:
            flux_table = self.cut_flux_table(flux_table)

        # Calculate get_features
        features = self.compute_stats(flux_table)

        # Calculate light_curve features
        lcfeat = self.extract_lightcurve_features(flux_table)
        features.update(lcfeat)
        # features.update(
        #    self.extract_lightcurve_features(flux_table)
        #    )

        # Create averaged values
        avgfeat = self.average_filtervalues(features)
        features.update(avgfeat)

        # if self.do_testplot:
        if features.get("success") and np.random.uniform() < self.plot_prob:
            self.test_plot(compound.get("stock"), flux_table, features)

        return features
