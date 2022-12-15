#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2DustEchoEval.py
# License           : BSD-3-Clause
# Author            : Eleni
# Date              : 01.09.2021
# Last Modified Date: 15.12.2022
# Last Modified By  : Simeon

from typing import Dict, Sequence, Tuple, List, Literal
from ampel.view.LightCurve import LightCurve
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.view.T2DocView import T2DocView
from ampel.struct.UnitResult import UnitResult
from ampel.ztf.util.ZTFIdMapper import to_ztf_id
from ampel.types import UBson
from ampel.enum.DocumentCode import DocumentCode

from scipy.signal import argrelextrema
from scipy.signal import find_peaks
import numpy as np
from nltk import flatten
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import uncertainties.umath as umath
import uncertainties.unumpy as unumpy
import os


class T2DustEchoEval(AbsTiedLightCurveT2Unit):
    """ """

    flux: bool = False

    directory: str

    rej_sigma: float = 5.0

    data_type: str = "wise"

    filters: Sequence[str]

    t2_dependency: Sequence[StateT2Dependency[Literal["T2BayesianBlocks"]]]

    def isMonotonic(self, A):
        return all(A[i] <= A[i + 1] for i in range(len(A) - 1)) or all(
            A[i] >= A[i + 1] for i in range(len(A) - 1)
        )

    def isDecreasing(self, A):
        return all(A[i] >= A[i + 1] for i in range(len(A) - 1))

    def isIncreasing(self, A):
        return all(A[i] <= A[i + 1] for i in range(len(A) - 1))

    PlotColor: Sequence[str] = ["red", "blue"]

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(
        self, light_curve: LightCurve, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:

        assert self.data_type in ["ztf_alert", "ztf_fp", "wise"]

        if self.data_type in ["ztf_alert", "ztf_fp"]:
            self.ztfid = to_ztf_id(light_curve.stock_id)
            self.PlotColor = ["green", "red", "orange"]

        else:
            self.ztfid = None

        if not light_curve.photopoints:
            self.logger.error("LightCurve has no photopoint")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        if not t2_views:  # Should not happen actually, T2Processor catches that case
            self.logger.error("Missing tied t2 views")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        excess_region = {
            "excess_jd": ([]),
            "excess_mag": ([]),
            "baseline_jd": (
                []
            ),  # the jd of the last datapoint of the baseline before peak excess region
            "baseline_mag": (
                []
            ),  # the mag of the last datapoint of the baseline before peak excess region
            "start_baseline_jd": (
                []
            ),  # the begining of baseline before the peak exces region
            "max_mag": ([]),
            "max_mag_jd": ([]),
            "significance": ([]),
            "strength_sjoert": ([]),
            "strength": ([]),
            "e_rise": ([]),
            "e_fade": ([]),
        }

        t2_output = {"description": ([]), "values": ([])}

        if self.flux:
            intensity_low_limit = 0.0
            intensity_high_limit = 100000000000000000.0
        else:
            intensity_low_limit = 7.0
            intensity_high_limit = 15.0

        for t2_view in t2_views:
            self.filters_lc = self.filters
            if t2_view.unit == "T2BayesianBlocks":
                self.logger.debug("Parsing t2 results from {}".format(t2_view.unit))
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )

                for key in self.filters:
                    if key not in t2_res.keys():
                        self.filters_lc.remove(key)

                    elif not t2_res[key]:
                        self.filters_lc.remove(key)

                if all(
                    (
                        t2_res[key]["nu_of_excess_regions"][0] == 0
                        and t2_res[key]["nu_of_baseline_blocks"][0] <= 1
                    )
                    for key in self.filters_lc
                ):
                    t2_output["description"].append("Only baseline")
                    for key in self.filters_lc:
                        excess_region["max_mag"].append(t2_res[key]["baseline"])
                else:

                    if t2_res["start_excess"] != None or t2_res["size_excess"] > 1.0:

                        if not t2_output["description"] and all(
                            flatten(min(t2_res[key]["max_mag_excess_region"]))[0]
                            < intensity_low_limit
                            or flatten(min(t2_res[key]["max_mag_excess_region"]))[0]
                            > intensity_high_limit
                            for key in self.filters_lc
                        ):
                            t2_output["description"].append(
                                "Very low or high magnitude"
                            )

                        elif all(
                            max(flatten(t2_res[key]["sigma_from_baseline"][0]))
                            < self.rej_sigma
                            for key in self.filters_lc
                        ):
                            t2_output["description"].append("Very low sigma")

                        elif t2_res["coincide_peak_block"] == -1:
                            t2_output["description"].append(
                                "No coincide excess regions"
                            )

                        elif any(
                            t2_res[key]["strength_sjoert"][0]
                            < t2_res[key]["significance"][0]
                            for key in self.filters_lc
                        ):
                            t2_output["description"].append("Low significance")

                        if not t2_output["description"]:
                            #  Check the excess region
                            #  Recalculate the excess region
                            maybe_interesting = False

                            for key in self.filters_lc:
                                if t2_res[key]["nu_of_excess_regions"][0] != 0:
                                    if self.flux:
                                        idx = np.argmax(
                                            t2_res[key]["max_mag_excess_region"]
                                        )
                                    else:
                                        idx = np.argmin(
                                            t2_res[key]["max_mag_excess_region"]
                                        )
                                    peak = t2_res[key]["jd_excess_regions"][idx]
                                    excess_region["max_mag"].append(
                                        t2_res[key]["max_mag_excess_region"][idx]
                                    )
                                    excess_region["max_mag_jd"].append(
                                        t2_res[key]["max_jd_excess_region"][idx]
                                    )
                                    ######################################################
                                    # Check if we have a stage transition case (time scale > 3years and nu of baye blocks <=1
                                    #                                        if (t2_res[key]['jd_excess_regions'][idx][-1]-t2_res[key]['jd_excess_regions'][idx][0] >= 1095. and t2_res[key]['nu_of_excess_blocks'][idx] <= 1) or t2_res[key]['max_baye_block_timescale'][0] >= 1095.:
                                    if (
                                        t2_res[key]["jd_excess_regions"][idx][-1]
                                        - t2_res[key]["jd_excess_regions"][idx][0]
                                        >= 1095.0
                                        and t2_res[key]["nu_of_excess_blocks"][idx] <= 1
                                    ):
                                        t2_output["description"].append(
                                            "Stage transition"
                                        )
                                    else:
                                        # fist check if there is a baseline before increase
                                        # Check the fluctuation inside the excess region
                                        #########################################################################
                                        # check if there is magnitude flactuanion before the excess region
                                        difference = []
                                        for baseline in np.array(
                                            t2_res[key]["jd_baseline_regions"]
                                        ):
                                            difference.append(peak[0] - baseline[-1])
                                        diff, position = min(
                                            (
                                                (b, nu)
                                                for nu, b in enumerate(
                                                    np.array(difference)
                                                )
                                                if b > 0
                                            ),
                                            default=(None, None),
                                        )

                                        # if diff is None or (diff is not None and len([value for value in difference if value > 0]) == 1):
                                        if diff is None:
                                            t2_output["description"].append(
                                                "Only declination"
                                            )
                                            excess_region["baseline_jd"].append(0)
                                            excess_region["start_baseline_jd"].append(0)
                                            excess_region["baseline_mag"].append(0)
                                        else:
                                            if (
                                                min(
                                                    peak[0] - value[-1]
                                                    for value in np.array(
                                                        t2_res[key][
                                                            "jd_baseline_regions"
                                                        ]
                                                    )
                                                    if peak[0] - value[-1] > 0.0
                                                )
                                                > 500
                                            ):
                                                # Move the gap closer to excess region:
                                                baseline_jd = (
                                                    t2_res[key]["jd_excess_regions"][
                                                        idx
                                                    ][0]
                                                    - 182.0
                                                )
                                                maybe_interesting = True
                                            else:
                                                baseline_jd = t2_res[key][
                                                    "jd_baseline_regions"
                                                ][position][-1]
                                            baseline_mag = t2_res[key][
                                                "mag_edge_baseline"
                                            ][position][0]
                                            t2_output["description"].append(
                                                "Baseline before excess region"
                                            )

                                            excess_region["baseline_jd"].append(
                                                baseline_jd
                                            )
                                            excess_region["start_baseline_jd"].append(
                                                t2_res[key]["jd_baseline_regions"][
                                                    position
                                                ][0]
                                            )
                                            excess_region["baseline_mag"].append(
                                                baseline_mag
                                            )
                                        #########################################################################
                                        #        if any(fluctuation < 3 for fluctuation in t2_res[key]['significance_of_fluctuation'][0]):
                                        #            t2_output['description'].append('Fluctuation before peak')
                                        #        else:
                                        # Check if there is magnitude fluctuations after the excess region
                                        difference = []
                                        for baseline in t2_res[key][
                                            "jd_baseline_regions"
                                        ]:
                                            difference.append(peak[-1] - baseline[0])
                                        diff, position = min(
                                            (
                                                (b, nu)
                                                for nu, b in enumerate(
                                                    np.array(difference)
                                                )
                                                if b < 0
                                            ),
                                            default=(None, None),
                                        )
                                        if (
                                            diff is not None
                                            and t2_res[key]["description"][idx][0] == 0
                                        ):
                                            # Outlier in the middle of the LC, Outlier in the last epoch will not be marked as Outlier
                                            t2_output["description"].append("Outlier")
                                        elif (
                                            diff is None
                                            and t2_res[key]["description"][idx][0] == 0
                                        ):
                                            maybe_interesting = True
                                            t2_output["description"].append(
                                                "Excess region with one data point in the end of LC"
                                            )
                                        elif t2_res[key]["description"][idx][0] != 0:
                                            t2_output["description"].append(
                                                "Excess region exists"
                                            )

                                        excess_jd = t2_res[key]["jd_excess_regions"][
                                            idx
                                        ][-1]
                                        excess_mag = t2_res[key]["mag_edge_excess"][
                                            idx
                                        ][-1]
                                        excess_region["excess_jd"].append(excess_jd)
                                        excess_region["excess_mag"].append(excess_mag)
                                #################################################################################
                                else:
                                    t2_output["description"].append("Only baseline")
                                    excess_region["max_mag"].append(
                                        t2_res[key]["baseline"]
                                    )
                    else:
                        t2_output["description"].append(
                            "Only baseline"
                        )  # no excess region
                        for key in self.filters_lc:
                            excess_region["max_mag"].append(t2_res[key]["baseline"])

        ##################################################################################
        if (
            len(excess_region["excess_jd"]) >= 2
            and len(excess_region["baseline_jd"]) >= 2
            and not all(value == 0 for value in excess_region["baseline_jd"])
            and not (
                any(
                    "Outlier" in description for description in t2_output["description"]
                )
            )
        ):
            t2_output["values"].append(excess_region)

            for fid, passband in enumerate(self.filters_lc, 1):
                if self.flux:
                    if self.data_type in ["ztf_alert", "ztf_fp"]:
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "flux_Jy", "flux_err_Jy"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": passband,
                                },
                            ],
                        )

                        if phot_tuple is None:
                            continue
                    elif self.data_type == "wise":
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "mean_flux", "flux_rms"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": passband,
                                },
                                {
                                    "attribute": "flux_ul",
                                    "operator": "==",
                                    "value": "False",
                                },
                            ],
                        )

                else:
                    phot_tuple = light_curve.get_ntuples(
                        ["jd", "magpsf", "sigmapsf"],
                        [
                            {
                                "attribute": "filter",
                                "operator": "==",
                                "value": passband,
                            },
                            {"attribute": "mag_ul", "operator": "==", "value": "False"},
                        ],
                    )

                df = pd.DataFrame(phot_tuple, columns=["jd", "mag", "mag.err"])
                df = df.sort_values(by=["jd"], ignore_index=True)

                if excess_region["baseline_jd"][fid - 1] != 0:
                    excess = df[
                        (df["jd"] > excess_region["baseline_jd"][fid - 1])
                        & (df["jd"] <= excess_region["excess_jd"][fid - 1])
                    ]
                    excess_after_max = excess[
                        excess["jd"] >= excess_region["max_mag_jd"][fid - 1][0]
                    ]
                    baseline_region = df[
                        (df["jd"] >= excess_region["start_baseline_jd"][fid - 1])
                        & (df["jd"] <= excess_region["baseline_jd"][fid - 1])
                    ]

                    baseline = np.mean(
                        unumpy.uarray(
                            np.array(baseline_region["mag"].values),
                            np.array(baseline_region["mag.err"].values),
                        )
                    ).nominal_value

                    baseline_rms = t2_res[passband]["baseline_rms"][0]
                    baseline_sigma = t2_res[passband]["baseline_sigma"][0]
                    excess_region["significance"].append(
                        t2_res[passband]["significance"][0]
                    )
                    excess_region["strength_sjoert"].append(
                        t2_res[passband]["strength_sjoert"][0]
                    )
                    excess_region["strength"].append(
                        abs(
                            t2_res[passband]["baseline"][0]
                            - excess_region["max_mag"][fid - 1][0]
                        )
                        / baseline_sigma
                    )

                    if self.flux:
                        excess_region["e_rise"].append(
                            (
                                excess_region["max_mag_jd"][fid - 1][0]
                                - excess_region["baseline_jd"][fid - 1]
                            )
                            / np.log(excess_region["max_mag"][fid - 1][0] / baseline)
                        )
                        if excess["mag"][excess["jd"].idxmax()] == max(excess["mag"]):
                            excess_region["e_fade"].append(0)
                        else:
                            excess_region["e_fade"].append(
                                -(
                                    excess_after_max["jd"].values[-1]
                                    - excess_after_max["jd"].values[0]
                                )
                                / np.log(
                                    excess_after_max["mag"].values[-1]
                                    / excess_after_max["mag"].values[0]
                                )
                            )
                    else:
                        excess_region["e_rise"].append(
                            -(
                                excess_region["max_mag_jd"][fid - 1][0]
                                - excess_region["baseline_jd"][fid - 1]
                            )
                            / np.log(excess_region["max_mag"][fid - 1][0] / baseline)
                        )
                        if excess["mag"][excess["jd"].idxmax()] == min(excess["mag"]):
                            excess_region["e_fade"].append(0)
                        else:
                            excess_region["e_fade"].append(
                                (
                                    excess_after_max["jd"].values[-1]
                                    - excess_after_max["jd"].values[0]
                                )
                                / np.log(
                                    excess_after_max["mag"].values[-1]
                                    / excess_after_max["mag"].values[0]
                                )
                            )
                else:
                    baseline_rms = t2_res[passband]["baseline_rms"][0]
                    baseline_sigma = t2_res[passband]["baseline_sigma"][0]
                    excess_region["significance"].append(
                        t2_res[passband]["significance"][0]
                    )
                    excess_region["strength_sjoert"].append(
                        t2_res[passband]["strength_sjoert"][0]
                    )
                    excess_region["strength"].append(
                        abs(
                            t2_res[passband]["baseline"][0]
                            - excess_region["max_mag"][fid - 1][0]
                        )
                        / baseline_sigma
                    )

                    excess_region["e_rise"].append("nan")
                    excess_region["e_fade"].append("nan")

            if all(
                value < 1000.0
                for value in excess_region["e_rise"]
                if value is not "nan"
            ) and all(
                value < 5000.0 and value > 400.0
                for value in excess_region["e_fade"]
                if value is not "nan"
            ):
                if maybe_interesting is False:
                    if any(value == "nan" for value in excess_region["e_rise"]):
                        t2_output["status"] = "3"  #
                    else:
                        t2_output["status"] = "1"
                else:
                    if any(value == "nan" for value in excess_region["e_rise"]):
                        t2_output["status"] = "3_maybe_interesting"  #
                    else:
                        t2_output["status"] = "1_maybe_interesting"
            else:
                if maybe_interesting is False:
                    if any(value == "nan" for value in excess_region["e_rise"]):
                        t2_output["status"] = "4"  #
                    else:
                        t2_output["status"] = "2"  #
                else:
                    if any(value == "nan" for value in excess_region["e_rise"]):
                        t2_output["status"] = "4_maybe_interesting"  #
                    else:
                        t2_output["status"] = "2_maybe_interesting"  #
        else:
            t2_output["status"] = "No further investigation"
            t2_output["values"].append(excess_region)

        #################### Plotting part #######################
        if t2_output["status"] != "No further investigation":
            fig, ax = plt.subplots(2, figsize=(18, 15))

            has_data = {}
            for fid, passband in enumerate(self.filters_lc, 1):
                if self.flux:
                    if self.data_type in ["ztf_alert", "ztf_fp"]:
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "flux_Jy", "flux_err_Jy"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": passband,
                                },
                            ],
                        )

                        if phot_tuple is None:
                            has_data[fid] = False
                            continue
                        else:
                            has_data[fid] = True

                    elif self.data_type == "wise":
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "mean_flux", "flux_rms"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": passband,
                                },
                                {
                                    "attribute": "flux_ul",
                                    "operator": "==",
                                    "value": "False",
                                },
                            ],
                        )
                        has_data[fid] = True
                else:
                    phot_tuple = light_curve.get_ntuples(
                        ["jd", "magpsf", "sigmapsf"],
                        [
                            {
                                "attribute": "filter",
                                "operator": "==",
                                "value": passband,
                            },
                            {"attribute": "mag_ul", "operator": "==", "value": "False"},
                        ],
                    )
                    has_data[fid] = True

                df = pd.DataFrame(phot_tuple, columns=["jd", "mag", "mag.err"])
                df["Outlier"] = False

                if t2_res[passband]["jd_outlier"]:
                    for block in t2_res[passband]["jd_outlier"]:
                        df.loc[
                            (df["jd"] >= block[0]) & (df["jd"] <= block[1]), "Outlier"
                        ] = True
                    outlier_datapoints = df[df["Outlier"] == True]
                    ax[0].errorbar(
                        outlier_datapoints["jd"] - 2400000.5,
                        outlier_datapoints["mag"],
                        yerr=outlier_datapoints["mag.err"],
                        label="Outlier",
                        fmt="s",
                        color="lightgray",
                        ecolor="lightgray",
                        markeredgecolor="black",
                        elinewidth=3,
                        capsize=0,
                    )

                datapoints = df[df["Outlier"] == False]
                ax[0].errorbar(
                    datapoints["jd"] - 2400000.5,
                    datapoints["mag"],
                    yerr=datapoints["mag.err"],
                    label=passband,
                    fmt="s",
                    color=self.PlotColor[fid - 1],
                    ecolor="lightgray",
                    markeredgecolor="black",
                    elinewidth=3,
                    capsize=0,
                )

            for fid, key in enumerate(self.filters_lc, 1):
                if has_data[fid] == False:
                    continue
                if self.flux:
                    idx = np.argmax(t2_res[key]["max_mag_excess_region"])
                else:
                    idx = np.argmin(t2_res[key]["max_mag_excess_region"])

                ax[0].axvline(
                    x=(t2_res[key]["jd_excess_regions"][idx][0] - 40) - 2400000.5,
                    label="T2 guess " + str(key),
                    color=self.PlotColor[fid - 1],
                )
                ax[0].axvline(
                    x=(t2_res[key]["jd_excess_regions"][idx][-1] + 40) - 2400000.5,
                    color=self.PlotColor[fid - 1],
                )
                ax[0].axvspan(
                    (t2_res[key]["jd_excess_regions"][idx][0] - 40) - 2400000.5,
                    (t2_res[key]["jd_excess_regions"][idx][-1] + 40) - 2400000.5,
                    alpha=0.05,
                    color=self.PlotColor[fid - 1],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"],
                    label="Baseline " + str(key),
                    ls="-",
                    color=self.PlotColor[fid - 1],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"][0] - t2_res[key]["baseline_sigma"][0],
                    label="Baseline sigma " + str(key),
                    ls="--",
                    color=self.PlotColor[fid - 1],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"][0] + t2_res[key]["baseline_sigma"][0],
                    label="Baseline sigma " + str(key),
                    ls="--",
                    color=self.PlotColor[fid - 1],
                )

            if self.flux:
                ax[0].set_ylabel("Flux", fontsize=20)
            else:
                ax[0].invert_yaxis()
                ax[0].set_ylabel("Difference magnitude", fontsize=20)
            ax[0].spines["left"].set_linewidth(3)
            ax[0].spines["bottom"].set_linewidth(3)
            ax[0].spines["top"].set_linewidth(3)
            ax[0].spines["right"].set_linewidth(3)
            ax[0].set_xlabel("MJD", fontsize=20)
            ax[0].tick_params(
                axis="both", which="major", labelsize=18, pad=10, length=10, width=3
            )

            handles, labels = ax[0].get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax[0].legend(
                by_label.values(),
                by_label.keys(),
                fontsize="x-large",
                loc=2,
                bbox_to_anchor=(1.05, 1),
                borderaxespad=0.0,
            )

            ax[1].patch.set_facecolor("white")
            ax[1].axes.xaxis.set_visible(False)
            ax[1].axes.yaxis.set_visible(False)

            position = ([0.1, 0.8], [0.25, 0.8], [0.4, 0.8])
            ax[1].text(
                0.1, 0.8, "Final class: " + str(t2_output["status"]), fontsize=25
            )
            ax[1].text(
                position[0][0],
                position[0][-1] - 0.18,
                "Max sigma from baseline",
                fontsize=21,
            )
            #            ax[1].hlines(position[0][-1]-0.22, position[0][0], position[0][0]+0.4, color='red')
            ax[1].text(
                position[0][0], position[0][-1] - 0.26, "Max magnitude", fontsize=21
            )
            ax[1].text(position[0][0], position[0][-1] - 0.34, "Baseline", fontsize=21)
            ax[1].text(
                position[0][0],
                position[0][-1] - 0.42,
                "Time scale of excess region",
                fontsize=21,
            )
            ax[1].text(
                position[0][0], position[0][-1] - 0.50, "e-folding rise", fontsize=21
            )
            ax[1].text(
                position[0][0], position[0][-1] - 0.58, "e-folding fade", fontsize=21
            )

            for fid, key in enumerate(self.filters_lc, 1):

                if has_data[fid] == False:
                    continue

                if self.flux:
                    idx = np.argmax(t2_res[key]["max_mag_excess_region"])
                else:
                    idx = np.argmin(t2_res[key]["max_mag_excess_region"])
                timescale = (
                    t2_res[key]["jd_excess_regions"][idx][-1]
                    - t2_res[key]["jd_excess_regions"][idx][0]
                )
                ax[1].text(
                    position[fid - 1][0] + 0.4,
                    position[fid - 1][-1] - 0.2,
                    str(key) + "\n\n",
                    color=self.PlotColor[fid - 1],
                    fontsize=25,
                )
                ax[1].text(
                    position[fid - 1][0] + 0.4,
                    position[0][-1] - 0.18,
                    format(max(flatten(t2_res[key]["sigma_from_baseline"][0])), ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    position[fid - 1][0] + 0.4,
                    position[0][-1] - 0.26,
                    format(t2_res[key]["max_mag_excess_region"][idx][0], ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    position[fid - 1][0] + 0.4,
                    position[0][-1] - 0.34,
                    format(t2_res[key]["baseline"][0], ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    position[fid - 1][0] + 0.4,
                    position[0][-1] - 0.42,
                    format(timescale, ".2f"),
                    fontsize=21,
                )
                if excess_region["e_rise"][fid - 1] != "nan":
                    ax[1].text(
                        position[fid - 1][0] + 0.4,
                        position[0][-1] - 0.50,
                        format(excess_region["e_rise"][fid - 1], ".2f"),
                        fontsize=21,
                    )
                    ax[1].text(
                        position[fid - 1][0] + 0.4,
                        position[0][-1] - 0.58,
                        format(excess_region["e_fade"][fid - 1], ".2f"),
                        fontsize=21,
                    )
                else:
                    ax[1].text(
                        position[fid - 1][0] + 0.4,
                        position[0][-1] - 0.50,
                        "NaN",
                        fontsize=21,
                    )
                    ax[1].text(
                        position[fid - 1][0] + 0.4,
                        position[0][-1] - 0.58,
                        "NaN",
                        fontsize=21,
                    )

            output_dir = os.path.join(
                self.directory, "Stage_" + str(t2_output["status"])
            )
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

            plt.savefig(
                output_dir + "/" + str(light_curve.stock_id) + ".pdf",
                bbox_inches="tight",
            )
            plt.close()

        t2_output: dict[str, UBson] = t2_output
        return t2_output
