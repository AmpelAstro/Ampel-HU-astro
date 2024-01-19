#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2PhaseLimit.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 28.04.2021
# Last Modified Date: 21.03.2022
# Last Modified By  : mf@physik.hu-berlin.de

import gc
import os
from collections.abc import Iterable
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


class T2PhaseLimit(AbsStateT2Unit, AbsTabulatedT2Unit):
    """
    Lightcurve analysis tools assume a certain transient life-time (order weeks for SNe).
    These can be confused by spurious early/late detections and/or background
    variability at the position of the SN.
    This unit tries to identify the most likely time-range for a SN explosion and
    evaluate whether sufficient data for continued exploration exist.
    """

    # *Conservative* estimates of how fast the transient rises and declines (in days)
    # These times are derived relative to the *median* detection date,
    # and are thus taken as symmetric around this
    half_time: float

    # Min number of detections remaining in the target range for subsequent analysis
    min_det: int = 3

    # Rejection sigma.
    rej_sigma: float = 5.0  # A lower number will more aggressively reject data
    min_sigma: float = (
        0.0  # Lower limit for including photometry when calculating time-range
    )

    # Spurious light in references typically manifests as negative flux.
    # Will require the fraction of obs with negative to be lower than this fraction, in each filter.
    neg_frac_lim: float = 0.05

    # Max magnitude to consider (constant low flux can correspond to flux in reference)
    max_flux: float

    # Constrain time-range of flux above these ranges
    risedec_fractions: list = []

    # Plot
    plot_suffix: None | str
    plot_dir: None | str

    def _magtoflux(self, mag: float) -> float:
        return 10 ** (-((mag) - 25) / 2.5)

    def __init__(self, *args, **kwargs):
        if "max_flux" not in kwargs:
            if "max_mag" not in kwargs:
                max_mag = 22.5  # previous default value
            else:
                max_mag = kwargs.pop("max_mag")
            kwargs["max_flux"] = self._magtoflux(max_mag)
        super().__init__(*args, **kwargs)

    def _get_risedec(self, bandflux, fracs) -> dict:
        """

        Determine the rise and decline time within a certain flux fraction, relative to a peak time.
        """
        rcinfo = {}

        fmax = max(bandflux["flux"])

        bandtpeak = bandflux["time"][(bandflux["flux"] == fmax)]

        # Rise
        rflux = bandflux[(bandflux["time"] < bandtpeak)]

        for frac in fracs:
            fmask = rflux["flux"] < fmax * (1 - frac)
            if sum(fmask) == 0:
                # No measurement, limit
                rcinfo["rise_" + str(frac) + "_limit"] = float(
                    bandtpeak - min(rflux["time"])
                )
            else:
                rcinfo["rise_" + str(frac)] = float(
                    bandtpeak - max(rflux["time"][fmask])
                )

        # Declline
        dflux = bandflux[(bandflux["time"] > bandtpeak)]

        for frac in fracs:
            fmask = dflux["flux"] < fmax * (1 - frac)
            if sum(fmask) == 0:
                # No measurement, limit
                rcinfo["fall_" + str(frac) + "_limit"] = float(
                    max(dflux["time"]) - bandtpeak
                )
            else:
                rcinfo["fall_" + str(frac)] = float(
                    min(dflux["time"][fmask]) - bandtpeak
                )

        return rcinfo

    def process(
        self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
    ) -> UBson | UnitResult:
        """
        Iteratively reject datapoints with outlying phase estimate.
        Once completed, check whether remaining data is compatible with expectations.
        Note: This unit does not differentiate between filters.
        Note: Upper limits are currently not taken into consideration.
        :returns: dict with estimate of likely peak and median time,
                tentative start and end-date and whether the phase looks like a sn.
        E.g.:
        {
                't_start' : 243101,
                't_end' : 243199,
                't_masked_duration' : 40,
                't_peak' : 243120,
                't_median' : 243124,
                'mag_peak' : 19.,
                'max_neg_frac' : 0.02,
                'n_remaining' : 5,
                't2eval' : 'OK' | 'fail:duration' | 'fail:detections' | 'fail:neg_flux' | 'warning:neg_iband'
        }
        """

        flux_table = self.get_flux_table(datapoints)

        #
        fflux_table = flux_table[
            (flux_table["flux"] < self.max_flux)
            & (flux_table["flux"] < 1e10)
            & (np.abs(flux_table["flux"] / flux_table["fluxerr"]) > self.min_sigma)
        ]

        # Enough data to even start evaluation?
        if len(fflux_table) < self.min_det:
            return {
                "t_start": None,
                "t_end": None,
                "t_masked_duration": None,
                "t_peak": None,
                "t_median": None,
                "flux_peak": None,
                "n_remaining": None,
                "t2eval": "fail:detections",
            }

        jd = fflux_table["time"]
        t_median = np.median(jd)
        mask = jd > 0

        # Iteratively reject datapoints until we reach a minimum (or find no more)
        # Using the median absolute deviation
        while sum(mask) >= self.min_det:
            sig_est = 1.48 * np.median(np.abs(jd[mask] - t_median))
            new_mask = np.abs(jd - t_median) < self.rej_sigma * sig_est
            if sum(new_mask) < sum(mask):
                mask = new_mask
                t_median = np.median(jd[mask])
            else:
                break

        if sum(mask) > 0:
            t_masked_duration = np.max(jd[mask]) - np.min(jd[mask])
        else:
            t_masked_duration = 0

        # Based on the median, define course phase range
        t_start = t_median - 2 * self.half_time
        t_end = t_median + 2 * self.half_time

        # (Re) retrieve data and magnitudes in this range
        # dps = light_curve.get_tuples('jd','magpsf',filters=
        # 	[{'attribute': 'jd', 'operator': '>', 'value': t_start},
        # 	{'attribute': 'jd', 'operator': '<', 'value': t_end},
        # 	{'attribute': 'magpsf', 'operator': '>', 'value': 0},
        # 	] )
        dps = flux_table[
            (flux_table["time"] > t_start)
            & (flux_table["time"] < t_end)
            & (flux_table["flux"] < 1e10)
        ]  # do we need to check for magpsf>0?
        if len(dps) > 0:
            dps.sort("flux")
            dps.reverse()
            t_peak = dps["time"][0]
            flux_peak = dps["flux"][0]
            n_remaining = len(dps)
        else:
            n_remaining = 0
            t_peak = None
            flux_peak = None

        # Tighter phase range based on peak time
        if t_peak is not None:
            t_start = min([t_median, t_peak]) - self.half_time
            t_end = max([t_median, t_peak]) + self.half_time
        else:
            t_start = t_median - self.half_time
            t_end = t_median + self.half_time

        # Investigate rise/fall ranges - use significant flux detections for this
        risedec_dict = {}
        if flux_peak is not None:
            for band in np.unique(fflux_table["band"]):
                # Get flux subset below flux limit
                risedec_dict[f"risedec_{band}"] = self._get_risedec(
                    fflux_table[(fflux_table["band"] == band)], self.risedec_fractions
                )

        # Investigate negative flux
        neg_frac_bands = []
        for filtid in np.unique(flux_table["band"]):
            in_band = flux_table[flux_table["band"] == filtid]
            len_ispos = len(in_band[in_band["flux"] > 0])
            filter_frac = 1 - float(len_ispos) / len(in_band)
            if filter_frac > self.neg_frac_lim:
                neg_frac_bands.append(filtid)

        # Evaluate
        if any(band.endswith("g") for band in neg_frac_bands):
            t2eval = "fail:neg_flux"
        elif any(band.endswith("r") for band in neg_frac_bands):
            t2eval = "fail:neg_flux"
        elif t_masked_duration > 4 * self.half_time:
            t2eval = "fail:duration"
        elif n_remaining < self.min_det:
            t2eval = "fail:detections"
        elif any(band.endswith("i") for band in neg_frac_bands):
            t2eval = "warning:neg_iband"
        else:
            t2eval = "OK"

        # Do plot?
        if self.plot_suffix and self.plot_dir:
            # How to construct name?
            tname = compound.get("stock")
            # Need plotting tools to define id mapper
            # tname = ZTFIdMapper.to_ext_id(light_curve.stock_id)

            fig = plt.figure(figsize=(6, 5))

            all_jd = flux_table["time"]
            bins: Any
            _1, bins, _2 = plt.hist(all_jd, bins=100, label="All alerts")

            plt.hist(jd, bins=bins, label="Clipped det.")
            plt.hist(jd[mask], bins=bins, label="Retained det.")

            plt.axvline(x=t_start, color="k", linestyle="dashed", label="Start/end")
            plt.axvline(x=t_end, color="k", linestyle="dashed")

            if t_peak is not None:
                plt.axvline(x=t_peak, label="t(peak)")
            plt.axvline(x=t_median, label="t(median)")

            plt.title(f"{tname} {t2eval}")
            plt.legend(loc="best")

            plt.tight_layout()
            plt.savefig(
                os.path.join(self.plot_dir, f"t2phaselimit_{tname}.{self.plot_suffix}")
            )

            plt.close("fig")
            plt.close("all")
            del fig
            gc.collect()

        return {
            "t_start": t_start,
            "t_end": t_end,
            "t_masked_duration": t_masked_duration,
            "t_peak": t_peak,
            "t_median": t_median,
            "flux_peak": flux_peak,
            "n_remaining": n_remaining,
            "t2eval": t2eval,
            "t_risedec": risedec_dict,
        }
