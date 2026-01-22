#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2TabulatorRiseDecline.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                28.12.2018
# Last Modified Date:  03.08.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import os
import random
import re
import string
import sys

# Exponential fit frequently returns overflow warnings
import warnings
from collections.abc import Iterable, Sequence
from contextlib import suppress
from typing import Any

import light_curve
import numpy as np
from astropy.table import Table
from light_curve.light_curve_py import RainbowFit
from scipy.interpolate import make_smoothing_spline
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.base.LogicalUnit import LogicalUnit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document

warnings.filterwarnings("ignore")


# Parameter limits - fits exceedings these are not considered
RISEDEC_LIMITS: dict[str, dict] = {}

# Tau fit
RISEDEC_LIMITS["tau"] = {
    "fall": [0.01, 75],
    "rise": [0.01, 75],
    "exp": [-10, 10],
}
# Spline fit
RISEDEC_LIMITS["spline"] = {
    "offsetratio": [0.0, 10.0],
    "ratioevo": [-99.0, 10.0],
    "splineflux": [-99.0, 100000.0],
    #    'peaksflux':[-99.,10000.], # Not yet implemented as splineflux check comes first and "should" be sufficient
}
# Bandfeatures
RISEDEC_LIMITS["band"] = {
    "rise_slope": [-100.0, 1000.0],
    #    'rise_slopesig':[0.,100.], # Not yet implemented - again assumed slope check sufficient.
    "fall_slope": [-100.0, 1000.0],
}
# Limits from lightkurve package
RISEDEC_LIMITS["lightkurve"] = {
    "linear_fit_slope": [-100.0, 100.0],
    "linear_fit_slope_sigma": [-100.0, 10.0],
    "bazin_fit_baseline": [-100.0, 100.0],
    "bazin_fit_amplitude": [-99.0, 10000.0],
    "bazin_fit_reduced_chi2": [-99.0, 100.0],
    "bazin_fit_rise_time": [-99.0, 1000.0],
    "bazin_fit_fall_time": [-99.0, 1000.0],
    "standard_deviation": [-99.0, 1000.0],
    "maximum_slope": [-99.0, 10000.0],
    "villar_fit_amplitude": [-99.0, 10000.0],
    "villar_fit_baseline": [-99.0, 1000.0],
    "villar_fit_rise_time": [-99.0, 1000.0],
    "villar_fit_fall_time": [-99.0, 1000.0],
    "villar_fit_plateau_duration": [-99.0, 1000.0],
    "villar_fit_reduced_chi2": [-99.0, 100.0],
    "excess_variance": [-99.0, 100.0],
    "kurtosis": [-99.0, 100.0],
    "period_0": [-99.0, 300.0],
    "eta_e": [-99.0, 1000000.0],
}
RISEDEC_LIMITS["rainbow"] = {
    "rainbow_rise_time": [-99.0, 300.0],
    "rainbow_t_color": [-99.0, 300.0],
    "rainbow_Tmin": [-99.0, 100000.0],
    "rainbow_Tmax": [-99.0, 100000.0],
    "rainbow_fall_time": [-99.0, 300.0],
    "rainbow_amplitude": [-99.0, 10000.0],
}


# Define the exponential model function
def exponential_model(x, a, b):
    return a * np.exp(b * x)


# Define the exponential rise and fill model based on Villar et al.
def supernova_villar_model(t, A, t_0, tau_rise, tau_fall):
    return A * np.exp(-(t - t_0) / tau_fall) / (1 + np.exp(-(t - t_0) / tau_rise))


# Function to fit the Villar model to the data
def fit_supernova_villar(t, f, f_err, peak_uncertainty, debugplot=False):
    lower_bounds = [
        0,
        -peak_uncertainty * 10,
        0,
        0,
    ]  # Lower bounds for A, t_0, tau_rise, tau_fall
    upper_bounds = [
        np.inf,
        peak_uncertainty * 10,
        100,
        100,
    ]  # Upper bounds for A, t_0, tau_rise, tau_fall

    # Perform curve fitting
    popt, pcov = curve_fit(
        supernova_villar_model,
        t,
        f,
        sigma=f_err,
        bounds=(lower_bounds, upper_bounds),
        maxfev=10000,
    )

    # popt contains the optimized parameters: A, t_0, tau_rise, and tau_fall
    A, t_0, tau_rise, tau_fall = popt

    # Calculate the fitted values
    z_fit = supernova_villar_model(t, *popt)
    chi2dof = sum((z_fit - f) ** 2 / f_err**2) / len(t)

    if debugplot:
        # Plot debug figure ...
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.figure()
        plt.scatter(t, f, label="Data", color="blue")
        plt.plot(t, z_fit, label="Fitted Villar Model", color="red")
        plt.title(
            f"t_0 = {t_0:.2f}, tau_rise = {tau_rise:.2f}, tau_fall = {tau_fall:.2f}, chi= {chi2dof:.2f}, len= {len(t)}"
        )
        plt.xlabel("Time")
        plt.ylabel("Brightness (z)")
        plt.legend()

        def make_unique(filename):
            return f"{filename.rsplit('.', 1)[0]}_{''.join(random.choices(string.ascii_lowercase + string.digits, k=4))}.{filename.rsplit('.', 1)[1]}"

        plt.savefig(make_unique("expdebug.png"))
        plt.close()

    # Optionally, return the parameters and the fitted curve
    return A, t_0, tau_rise, tau_fall, z_fit, chi2dof, pcov


# Function to fit exponential curve to data
def fit_exponential_rise(x, z):
    # Perform curve fitting
    popt, _ = curve_fit(exponential_model, x, z)

    # popt contains the optimized parameters a and b
    a, b = popt

    # Calculate the fitted values
    z_fit = exponential_model(x, *popt)

    # Optionally, return the parameters and the fitted curve
    return a, b, z_fit


def check_lightkurve(fitmethod: str, fitresults: float) -> bool:
    """
    Investigate outcome from lightkurve package fit.
    Will look whether limits for the {fitmethod} is provided in RISEDEC_LIMITS.
    If so, check whether the fit results are within these limits.
    Return False if limits exists and fit results are outside these limits, otherwise True.
    """

    if fitmethod == "rainbow":
        sys.exit("should not be here")

    if fitmethod not in RISEDEC_LIMITS["lightkurve"]:
        # No limits for this fit method
        return True

    return bool(
        RISEDEC_LIMITS["lightkurve"][fitmethod][0]
        < fitresults
        < RISEDEC_LIMITS["lightkurve"][fitmethod][1]
    )


def check_rainbow(fitresults: dict) -> bool:
    """
    Directly check the Rainbow fit method. Only accept if all limits are passed.
    """
    for prop, meas in fitresults.items():
        if prop not in RISEDEC_LIMITS["rainbow"]:
            continue
        if not (
            RISEDEC_LIMITS["rainbow"][prop][0]
            < meas
            < RISEDEC_LIMITS["rainbow"][prop][1]
        ):
            return False
    return True


def spline_analysis(
    tab: Table,
    band: str,
    spline_lam=0.1,
    do_plot=False,
    int_step=3.0,
    width=2,
    height=None,
    distance=3,
    wlen=None,
    plateau_size=None,
    dt_fluxevo: float = 20.0,
) -> dict:
    """
    Do weighted spline of band lc and extract features.
    tab assumed to be of one band, with time, flux and flux_unc fields.

    Input:
    tab: AstropyTable, containing the light curve data for one band with at least three columns:
        'jd': Julian Date (time)
        'flux': Flux values
        'fluxerr': Uncertainties in the flux
    band: name of band
    spline_lam: Regularization parameter for the spline smoothing (default is 0.1).
    do_plot: Whether or not to plot the spline and peaks (default is False).
    int_step: The interpolation step size (default is 3.).
    width, height, distance, wlen, plateau_size: Parameters for find_peaks, which help control how peaks are detected (these can be adjusted to suit the data).
    dt_fluxevo: Time offset for computing the peak flux evolution.


    Output:

    <band>_splineflux: The total flux from the spline interpolation over the entire light curve.
    <band>_splineerr:  The total error in the flux calculated from the spline of the flux uncertainties.
    <band>_splinethalf:  The duration (in time) between the points where the spline flux is greater than half of its maximum value.
    <band>_dethalf: The number of data points  that lie within the range of times where the spline flux is greater than half of its maximum value.
    <band>_peaksnbr: The number of peaks identified in the smoothed flux.
    <band>_peaksflux (list):  The flux values at the detected peaks.
    <band>_peaksjd (list): The Julian Date values at the detected peaks.
    <band>_offsetpeakflux (list): The flux values at the peaks, shifted by a specified time offset (dt_fluxevo).
    <band>_offsetpeakratio: The ratio of the flux at the shifted peak to the flux at the original peak.
    <band>_dtrise63: The time difference between the first peak and the point where the flux first rises to 63% of the peak flux (rise time).
    <band>_dtfall63: The time difference between the first peak and the point where the flux falls to 63% of the peak flux (fall time).
    <band>_peakstdiff: The time differences between consecutive peaks.
    <band>_valleyflux: The flux values at the valleys between peaks, calculated as the midpoint between consecutive peaks.

    """

    if len(tab) < 6:
        return {}

    spl = make_smoothing_spline(
        tab["time"], tab["flux"], w=1.0 / tab["fluxerr"] ** 2, lam=spline_lam
    )
    dspl = make_smoothing_spline(
        tab["time"], tab["fluxerr"], w=1.0 / tab["fluxerr"] ** 2, lam=spline_lam
    )
    tphase = np.arange(
        tab["time"].min() + int_step, tab["time"].max() - int_step, int_step
    )
    finterp = spl(tphase)

    if do_plot:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.plot(tphase, finterp, "-.")

    prominence = tab["flux"].max() / 5
    peaks, _ = find_peaks(
        finterp,
        distance=distance,
        prominence=prominence,
        width=width,
        wlen=wlen,
        plateau_size=plateau_size,
        height=height,
    )

    if do_plot:
        plt.plot(tphase, finterp)
        _ = [
            plt.axvline(p, color="k", alpha=0.3, linestyle="dashed")
            for p in tphase[peaks]
        ]
        plt.plot(tphase[peaks], finterp[peaks], "X", ms=12, alpha=0.5)

    # Collect features
    odict: dict[str, float | list[float]] = {}
    if len(peaks) > 0:
        # Proceed if integrated flux makes sense
        splineflux = float(sum(finterp * int_step))
        if not (
            RISEDEC_LIMITS["spline"]["splineflux"][0]
            < splineflux
            < RISEDEC_LIMITS["spline"]["splineflux"][1]
        ):
            return {}

        odict[band + "_splineflux"] = splineflux
        odict[band + "_splineerr"] = float(sum(dspl(tphase) * int_step))

        foo = tphase[finterp > (finterp.max() / 2)]
        odict[band + "_splinethalf"] = float(foo[-1] - foo[0])
        odict[band + "_dethalf"] = int(
            sum((tab["time"] > foo[0]) & (tab["time"] < foo[-1]))
        )
        odict[band + "_peaksnbr"] = len(peaks)
        odict[band + "_peaksflux"] = [float(f) for f in finterp[peaks]]
        odict[band + "_peakstime"] = [float(f) for f in tphase[peaks]]

        # Extracting flux at offset position
        try:
            odict[band + "_offsetpeakflux"] = [
                float(f) for f in finterp[peaks + int(dt_fluxevo / int_step)]
            ]
            odict[band + "_offsetpeakratio"] = [
                off / flux
                for off, flux in zip(
                    odict[band + "_offsetpeakflux"],  # type: ignore[arg-type]
                    odict[band + "_peaksflux"],  # type: ignore[arg-type]
                    strict=False,
                )
            ]
        except IndexError:
            # Happens if not enough elements after peak
            pass

        # Looking at rise and deline time to 63% flux (~0.5 mag) for first peak.
        foo = np.where(finterp / 0.63 < odict[band + "_peaksflux"][0])[0]  # type: ignore[index]
        try:
            odict[band + "_dtrise63"] = float(
                tphase[peaks[0]] - tphase[max(foo[foo < peaks[0]])]
            )
            odict[band + "_dtfall63"] = float(
                tphase[min(foo[foo > peaks[0]])] - tphase[peaks[0]]
            )
        except ValueError:
            # In case the array does not extend to these limits
            pass

        if len(peaks) > 1:
            odict[band + "_peakstdiff"] = [
                float(f) for f in tphase[peaks][1:] - tphase[peaks][0:-1]
            ]
            odict[band + "_valleyflux"] = [
                float(f)
                for f in finterp[
                    [int(f) for f in np.floor((peaks[1:] + peaks[0:-1]) / 2)]
                ]
            ]

    return odict


def spline_colevo(
    featuredict: dict, color_list: list[list[str]], dt_fluxevo: float = 20.0
) -> dict:
    """
    Check whether there are multiple bands available, sufficient to calculate color properties.

    Input Parameters:

    featuredict (type: dict):
        A dictionary containing the extracted features from the light curve analysis: peak fluxes (<band>_peaksflux) and offset peak fluxes (<band>_offsetpeakflux) for different bands.
    color_list (type: list of tuples):
        List of tuples for which color values will be calculated. Each tuple contains two band names, like ('g', 'r').
    dt_fluxevo (type: int, default is 20):
        The time step (in days) assumed used when calculating the flux evolution ratio (ratioevo).

    Output:
    cdict (type: dict):
            peakratio_<band1>_<band2>: The ratio of the first peak flux of <band1> to <band2>.
            offsetratio_<band1>_<band2>: The ratio of the offset peak flux of <band1> to <band2>.
            ratioevo_<band1>_<band2>: The ratio evolution, calculated as the difference between the peak ratio and the offset ratio, normalized by dt_fluxevo.

    """
    # Look for color terms
    # Again, we only do this for the first peak - different peak numbers in different bands would need more careful checking
    cdict = {}
    for colpair in color_list:
        if (
            colpair[0] + "_peaksflux" not in featuredict
            or colpair[1] + "_peaksflux" not in featuredict
        ):
            continue

        cdict["peakratio_" + colpair[0] + "_" + colpair[1]] = (
            featuredict[colpair[0] + "_peaksflux"][0]
            / featuredict[colpair[1] + "_peaksflux"][0]
        )
        try:
            if (
                RISEDEC_LIMITS["spline"]["offsetratio"][0]
                < (
                    offsetratio := featuredict[colpair[0] + "_offsetpeakflux"][0]
                    / featuredict[colpair[1] + "_offsetpeakflux"][0]
                )
                < RISEDEC_LIMITS["spline"]["offsetratio"][1]
            ):
                cdict["offsetratio_" + colpair[0] + "_" + colpair[1]] = offsetratio
            if (
                RISEDEC_LIMITS["spline"]["ratioevo"][0]
                < (
                    ratioevo := (
                        cdict["peakratio_" + colpair[0] + "_" + colpair[1]]
                        - offsetratio
                    )
                    / dt_fluxevo
                )
                < RISEDEC_LIMITS["spline"]["ratioevo"][1]
            ):
                cdict["ratioevo_" + colpair[0] + "_" + colpair[1]] = ratioevo
        except KeyError:
            # In case offsetpeak could not be defined
            pass

    return cdict


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


def floaty(value) -> float:
    """
    Make sure that we get a single float that can 
    be stored in a database.
    Reduction rule for non-scalars: mean of all values.
    """

    # Scalar numbers
    if isinstance(value, (int, float, np.number)):
        return float(value)

    # Pandas DataFrame or Series - not used here
    #if isinstance(value, (pd.DataFrame, pd.Series)):
    #    return float(value.to_numpy().mean())

    # NumPy array or array-like
    if isinstance(value, np.ndarray):
        return float(value.mean())

    raise TypeError(
        f"Unsupported type {type(value)}. "
        "Expected float, numpy array, pandas DataFrame or Series."
    )


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


class FitFailed(RuntimeError): ...


class T2TabulatorRiseDeclineBase:
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
    color_list: list[list[str]] = [
        ["lsstu", "lsstg"],
        ["lsstg", "lsstr"],
        ["lsstr", "lssti"],
        ["lssti", "lsstz"],
        ["lsstz", "lssty"],
        ["ztfg", "ztfr"],
        ["ztfr", "ztfi"],
    ]
    max_tgap: int = 30
    min_expfit_det: int = 4  # Fit exponential model if a band has this many detections or more (Warning: possibly time-consuming). Turn off by providing large value.

    # Cut the lightcurve if longer than this limit.
    # Motivated by worse classification for longer (inbalanced training?)
    # max_ndet: int = 20
    # For new training round, change this
    max_ndet: int = 200000

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

        # Get absolute height scale for peak finding
        heightscale = ftable["flux"].max() / 5

        for band in set(ftable["band"]):
            bt = ftable[ftable["band"] == band]

            max_flux = bt["flux"].max()
            max_flux_time = bt[bt["flux"] == max_flux]["time"][0]
            banddata["jd_peak_" + band] = max_flux_time

            # Time above half flux
            halft = bt[bt["flux"] > max_flux / 2]["time"]
            if len(halft) > 2:
                banddata["time_halfpeak_" + band] = halft.max() - halft.min()

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
                    if (
                        RISEDEC_LIMITS["band"]["rise_slope"][0]
                        < (rise := fit[1])
                        < RISEDEC_LIMITS["band"]["rise_slope"][1]
                    ):
                        banddata["rise_slope_" + band] = rise
                        banddata["rise_slopesig_" + band] = rise / np.sqrt(cov[1][1])

                except RuntimeError:
                    pass
            #                    self.logger.debug("Risetime curve fit failed.")
            if len(fallt) > 1:
                try:
                    fit, cov = curve_fit(
                        linearFunc,
                        fallt["time"] - tscale,
                        fallt["flux"],
                        sigma=fallt["fluxerr"],
                        absolute_sigma=True,
                    )
                    if (
                        RISEDEC_LIMITS["band"]["fall_slope"][0]
                        < (fall := fit[1])
                        < RISEDEC_LIMITS["band"]["fall_slope"][1]
                    ):
                        banddata["fall_slope_" + band] = fall
                        banddata["fall_slopesig_" + band] = fall / np.sqrt(cov[1][1])

                except RuntimeError:
                    pass
            #                    self.logger.debug("Falltime curve fit failed.")

            # Possible fit exponential components if sufficient data
            if len(bt) > self.min_expfit_det and len(fallt) > 1 and len(riset) > 1:
                # Fit exponetial rise/fall a la Villar
                try:
                    (
                        _,
                        _,
                        tau_rise,
                        tau_fall,
                        _,
                        chi_dof,
                        _,
                    ) = fit_supernova_villar(
                        bt["time"] - max_flux_time,
                        bt["flux"],
                        bt["fluxerr"],
                        peak_uncertainty=self.t_cadence,
                    )
                    # Stor parameters (not A nor t_0) if fit was successful
                    if (
                        RISEDEC_LIMITS["tau"]["rise"][0]
                        < tau_rise
                        < RISEDEC_LIMITS["tau"]["rise"][1]
                    ) and (
                        RISEDEC_LIMITS["tau"]["fall"][0]
                        < tau_fall
                        < RISEDEC_LIMITS["tau"]["fall"][1]
                    ):
                        banddata["tau_rise_" + band] = tau_rise
                        banddata["tau_fall_" + band] = tau_fall
                        banddata["tau_chidof_" + band] = chi_dof
                        banddata["tau_dof_" + band] = len(bt)

                except ValueError as exc:
                    # triggered in attempt to SVD final jacobian containing inf or NaN
                    if "array must not contain infs or NaNs" in exc.args:
                        pass
                    else:
                        raise
                except RuntimeError:
                    pass
            #                    self.logger.debug("Rise/fall fit failed.")
            elif len(bt) > self.min_expfit_det / 2:
                # Fit exponetial, if less data
                try:
                    _, tau, _ = fit_exponential_rise(
                        bt["time"] - max_flux_time, bt["flux"]
                    )
                    # Stor parameters (not A nor t_0)
                    if (
                        RISEDEC_LIMITS["tau"]["exp"][0]
                        < tau
                        < RISEDEC_LIMITS["tau"]["exp"][1]
                    ):
                        banddata["tau_exp_" + band] = tau
                        banddata["tau_dof_" + band] = len(bt)

                except RuntimeError:
                    pass
            #                    self.logger.debug("Exp fit failed.")

            # Check for flux decline until the time for a possible second bump
            if (
                sum(
                    isecpeak := (
                        np.abs(fallt["time"] - max_flux_time - self.dt_fluxevo)
                        < self.t_cadence / 2
                    )
                )
                > 0
            ):
                banddata["fluxevo_ratio_" + band] = (
                    np.mean(fallt["flux"][isecpeak]) / max_flux
                )

            # Run the spline analysis
            banddata.update(
                spline_analysis(
                    bt,
                    band,
                    int_step=self.t_cadence,
                    height=heightscale,
                    dt_fluxevo=self.dt_fluxevo,
                )
            )

        # Combine spline data into colors if available
        if len(self.color_list) > 0:
            banddata.update(
                spline_colevo(banddata, self.color_list, dt_fluxevo=self.dt_fluxevo)
            )

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

        try:
            det_table = flux_table[band_mask & sig_mask]
        except TypeError:
            # Can happen if no detections are found
            o["success"] = False
            o["cause"] = "Cannot create det table."
            return o

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

        o["mag_det"] = floaty(getMag(flux_table[flux_table["time"] == o["jd_det"]]))
        o["band_det_id"] = getBandBits(
            [flux_table[flux_table["time"] == o["jd_det"]]["band"][0]]
        )

        o["mag_last"] = floaty(getMag(flux_table[flux_table["time"] == o["jd_last"]]))
        o["band_last_id"] = getBandBits(
            [flux_table[flux_table["time"] == o["jd_last"]]["band"][0]]
        )

        # Check fails irregularly - this was due to getMag sometimes returning arrays.
        # Try w floaty
        #with suppress(TypeError):
        o["mag_min"] = floaty(
            getMag(flux_table[flux_table["flux"] == max(flux_table["flux"])])
        )

        #with suppress(TypeError):
        o["jd_min"] = floaty(flux_table["time"][np.argmax(flux_table["flux"])])

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


class BaseLightCurveFeatures(LogicalUnit):
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
    # Minimally requires effective band passes as method parameters:
    # {"band_wave_aa": {"ztgf": 4000, ...}
    rainbow_fit: None | dict[str, Any] = None

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
        # Load single band feature extractors
        self.fluxextractor = light_curve.Extractor(
            *(
                getattr(light_curve, k)(**(v or {}))
                for k, v in self.lightcurve_features_flux.items()
            )
        )
        # Some lightkurve (multiband? experimental?) features require special setup
        if self.rainbow_fit:
            self.rainbowextractor = RainbowFit.from_angstrom(
                self.rainbow_fit["band_wave_aa"], with_baseline=False
            )
        else:
            self.rainbowextractor = None

    def post_init(self) -> None:
        self.init_lightcurve_extractor()

    def extract_lightcurve_features(self, flux_table: Table) -> dict[str, float]:
        result = {}
        # Derive single band features
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
                            #                            self.lightcurve_features_flux.keys(),
                            self.fluxextractor.names,
                            self.fluxextractor(
                                in_band["time"], in_band["flux"], in_band["fluxerr"]
                            ),
                            strict=False,
                        )
                        if check_lightkurve(k, v)
                    }
                    result.update(lcout)
            except ValueError:
                self.logger.debug("lightkurve extract fail")
                pass
        # Multiband features (e.g. Rainbow)
        if self.rainbowextractor and len(flux_table["flux"]) > 8:
            try:
                lcout = self.rainbowextractor(
                    flux_table["time"],
                    flux_table["flux"],
                    sigma=flux_table["fluxerr"],
                    band=flux_table["band"],
                )
                rainbowout = dict(
                    zip(
                        ["rainbow_" + s for s in self.rainbowextractor.names],
                        lcout,
                        strict=False,
                    )
                )
                if check_rainbow(rainbowout):
                    result.update(rainbowout)
            except (OSError, RuntimeError):
                self.logger.debug("lightkurve rainbow extract fail")
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

        import matplotlib.pyplot as plt  # noqa: PLC0415

        bands = set(table["band"])

        _, axs = plt.subplots(2, 3)

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
            title = f"ndet: {t2result['ndet']} "

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
