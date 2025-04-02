#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2GetLensSNParameters.py
# License           : BSD-3-Clause
# Author            : alice.townsend@physik.hu-berlin.de
# Date              : 22.11.2021
# Last Modified Date: 19.10.2022
# Last Modified By  : alice.townsend@physik.hu-berlin.de


from bisect import bisect_left
from collections.abc import Sequence

import numpy as np
from astropy.table import Table

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2GetLensSNParameters(T2RunSncosmo):
    unc: float

    def _get_fit_metrics(self, sncosmo_result, sncosmo_table, sncosmo_model) -> dict:
        # -------------------------------- Overriden method -----------------------------#
        lc_metrics = super()._get_fit_metrics(
            sncosmo_result, sncosmo_table, sncosmo_model
        )
        # Record the closest detection to peak in each band
        t0 = sncosmo_model.get("t0")  # get time at peak of light curve
        sncosmo_table.add_column(
            "temporary", name="epoch", index=0
        )  # add column with temp values for epoch
        sncosmo_table.sort(["time"])  # sort data chronologically

        def take_closest(myList, myNumber):
            """
            Assumes myList is sorted. Returns closest value to myNumber.
            If two numbers are equally close, return the smallest number.
            """
            pos = bisect_left(myList, myNumber)
            if pos == 0:
                return myList[0]
            if pos == len(myList):
                return myList[-1]
            before = myList[pos - 1]
            after = myList[pos]
            if after - myNumber < myNumber - before:
                return after
            return before

        for band in np.unique(sncosmo_table["band"]):
            table_name = str(band) + "cut"
            globals()[table_name] = Table(
                sncosmo_table[0:0]
            )  # create empty table for each band
            band_table = sncosmo_table[
                (sncosmo_table["band"] == band)
            ]  # split data into different bands
            time_minus7 = take_closest(
                band_table["time"], t0 - 7
            )  # find closest time to t0-7 days in each band
            if (t0 - 7 - self.unc) <= time_minus7 <= (t0 - 7 + self.unc):
                band_table["epoch"][np.where(band_table["time"] == time_minus7)] = (
                    "minus7"
                )
                for row in band_table[
                    np.where(band_table["time"] == time_minus7)
                ]:  # iterate over temporary table that is a subset of data
                    globals()[table_name].add_row(row)
            time_closest = take_closest(
                band_table["time"], t0
            )  # find closest time to peak in each band
            if (t0 - self.unc) <= time_closest <= (t0 + self.unc):
                band_table["epoch"][np.where(band_table["time"] == time_closest)] = "t0"
                for row in band_table[
                    np.where(band_table["time"] == time_closest)
                ]:  # iterate over temporary table that is a subset of data
                    globals()[table_name].add_row(row)
            time_plus7 = take_closest(
                band_table["time"], t0 + 7
            )  # find closest time to t0+7 days in each band
            if (t0 + 7 - self.unc) <= time_plus7 <= (t0 + 7 + self.unc):
                band_table["epoch"][np.where(band_table["time"] == time_plus7)] = (
                    "plus7"
                )
                for row in band_table[
                    np.where(band_table["time"] == time_plus7)
                ]:  # iterate over temporary table that is a subset of data
                    globals()[table_name].add_row(row)

        # Calculate the colour at different epochs
        def calculate_colour_peak(band1, band2, epoch):
            """
            Calculate the colour from two different bands of a user's choosing.
            Bluer colour is band1.
            epoch is the epoch within which we are searching (i.e. t0, minus7, or plus7).
            """
            band1_table = str(band1) + "cut"
            band2_table = str(band2) + "cut"
            if band1 not in sncosmo_table["band"] or band2 not in sncosmo_table["band"]:
                return None, None
            if (
                epoch not in globals()[band1_table]["epoch"]
                or epoch not in globals()[band2_table]["epoch"]
            ):
                return None, None
            colour = -2.5 * np.log10(
                globals()[band1_table][
                    np.where(globals()[band1_table]["epoch"] == epoch)
                ]["flux"]
                / globals()[band2_table][
                    np.where(globals()[band2_table]["epoch"] == epoch)
                ]["flux"]
            )
            unc1 = (
                2.5
                * 0.434
                * (
                    globals()[band1_table][
                        np.where(globals()[band1_table]["epoch"] == epoch)
                    ]["fluxerr"]
                    / globals()[band1_table][
                        np.where(globals()[band1_table]["epoch"] == epoch)
                    ]["flux"]
                )
            )
            unc2 = (
                2.5
                * 0.434
                * (
                    globals()[band2_table][
                        np.where(globals()[band2_table]["epoch"] == epoch)
                    ]["fluxerr"]
                    / globals()[band2_table][
                        np.where(globals()[band2_table]["epoch"] == epoch)
                    ]["flux"]
                )
            )
            unc = (unc1**2 + unc2**2) ** 0.5
            return colour.data[0], unc.data[0]

        (
            lc_metrics["r_i_colour_peak"],
            lc_metrics["r_i_colour_peak_err"],
        ) = calculate_colour_peak("ztfr", "ztfi", "t0")
        (
            lc_metrics["r_i_colour_plus7"],
            lc_metrics["r_i_colour_plus7_err"],
        ) = calculate_colour_peak("ztfr", "ztfi", "plus7")
        (
            lc_metrics["r_i_colour_minus7"],
            lc_metrics["r_i_colour_minus7_err"],
        ) = calculate_colour_peak("ztfr", "ztfi", "minus7")
        (
            lc_metrics["g_r_colour_peak"],
            lc_metrics["g_r_colour_peak_err"],
        ) = calculate_colour_peak("ztfg", "ztfr", "t0")
        (
            lc_metrics["g_r_colour_plus7"],
            lc_metrics["g_r_colour_plus7_err"],
        ) = calculate_colour_peak("ztfg", "ztfr", "plus7")
        (
            lc_metrics["g_r_colour_minus7"],
            lc_metrics["g_r_colour_minus7_err"],
        ) = calculate_colour_peak("ztfg", "ztfr", "minus7")
        (
            lc_metrics["g_i_colour_peak"],
            lc_metrics["g_i_colour_peak_err"],
        ) = calculate_colour_peak("ztfg", "ztfi", "t0")
        (
            lc_metrics["g_i_colour_plus7"],
            lc_metrics["g_i_colour_plus7_err"],
        ) = calculate_colour_peak("ztfg", "ztfi", "plus7")
        (
            lc_metrics["g_i_colour_minus7"],
            lc_metrics["g_i_colour_minus7_err"],
        ) = calculate_colour_peak("ztfg", "ztfi", "minus7")

        # Calculate observed magnitude close to peak
        def calculate_obsmag_peak(band1, epoch):
            """
            Calculates observed magnitude for a particular epoch
            """
            band1_table = str(band1) + "cut"
            if band1 not in sncosmo_table["band"]:
                return None
            if epoch not in globals()[band1_table]["epoch"]:
                return None
            obs_mag = (
                -2.5
                * np.log10(
                    globals()[band1_table][
                        np.where(globals()[band1_table]["epoch"] == epoch)
                    ]["flux"]
                )
                + 25
            )
            return obs_mag.data[0]

        lc_metrics["obsmag_ztfg_peak"] = calculate_obsmag_peak("ztfg", "t0")
        lc_metrics["obsmag_ztfr_peak"] = calculate_obsmag_peak("ztfr", "t0")
        lc_metrics["obsmag_ztfi_peak"] = calculate_obsmag_peak("ztfi", "t0")
        lc_metrics["obsmag_ztfg_plus7"] = calculate_obsmag_peak("ztfg", "plus7")
        lc_metrics["obsmag_ztfr_plus7"] = calculate_obsmag_peak("ztfr", "plus7")
        lc_metrics["obsmag_ztfi_plus7"] = calculate_obsmag_peak("ztfi", "plus7")
        lc_metrics["obsmag_ztfg_minus7"] = calculate_obsmag_peak("ztfg", "minus7")
        lc_metrics["obsmag_ztfr_minus7"] = calculate_obsmag_peak("ztfr", "minus7")
        lc_metrics["obsmag_ztfi_minus7"] = calculate_obsmag_peak("ztfi", "minus7")

        thresholds = {
            "r_i_colour_peak": -0.28,
            "r_i_colour_plus7": -0.23,
            "r_i_colour_minus7": -0.34,
            "g_r_colour_peak": 0.12,
            "g_r_colour_plus7": 0.33,
            "g_r_colour_minus7": -0.08,
            "g_i_colour_peak": 0.06,
            "g_i_colour_plus7": 0.33,
            "g_i_colour_minus7": -0.08,
        }

        if any(
            lc_metrics.get(key) is not None and lc_metrics[key] > value
            for key, value in thresholds.items()
        ):
            lc_metrics["pass_colour_cuts"] = True
        else:
            lc_metrics["pass_colour_cuts"] = False

        return lc_metrics

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        t2_output = super().process(compound, datapoints, t2_views)
        if (
            "sncosmo_result" in t2_output
            and "fit_metrics" in t2_output["sncosmo_result"]
        ):
            colour_pass = t2_output["sncosmo_result"]["fit_metrics"]["pass_colour_cuts"]
            peak_pass = (
                t2_output["sncosmo_result"]["fit_metrics"]["restpeak_model_absmag_B"]
                < -19.5
            )
        else:
            colour_pass = False
            peak_pass = False

        ampel_z_output = self.get_ampelZ(t2_views)
        if "ampel_z" in ampel_z_output:
            redshift = ampel_z_output["ampel_z"]
            redshift_pass = redshift > 0.1
            dist = ampel_z_output["ampel_dist"]
            dist_pass = dist < 3.0
        else:
            redshift_pass = False
            dist_pass = False

        if colour_pass and peak_pass and redshift_pass and dist_pass:
            t2_output["pass_all_cuts"] = True
        else:
            t2_output["pass_all_cuts"] = False

        return t2_output
