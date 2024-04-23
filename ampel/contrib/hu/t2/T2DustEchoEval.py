#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2DustEchoEval.py
# License           : BSD-3-Clause
# Author            : Eleni
# Date              : 01.09.2021
# Last Modified Date: 10.01.2023
# Last Modified By  : Jannis <jannis.necker@desy.de>

import os
from collections.abc import Sequence
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # type: ignore
from uncertainties import unumpy  # type: ignore

from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.contrib.hu.util.flatten import flatten
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id


class T2DustEchoEval(AbsTiedLightCurveT2Unit):
    """ """

    flux: bool = False

    directory: str

    rej_sigma: float = 5.0

    data_type: str = "wise"

    filters: list[str]

    t2_dependency: Sequence[StateT2Dependency[Literal["T2BayesianBlocks"]]]

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #

    def process(
        self, light_curve: LightCurve, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:
        assert self.data_type in ["ztf_alert", "ztf_fp", "wise"]

        self.ztfid: None | str = None

        if self.data_type in ["ztf_alert", "ztf_fp"] and isinstance(
            light_curve.stock_id, int
        ):
            self.ztfid = to_ztf_id(light_curve.stock_id)

        self.plot_colors: dict[str, str] = {
            "Wise_W1": "red",
            "Wise_W2": "blue",
            "ZTF_g": "green",
            "ZTF_r": "red",
            "ZTF_i": "orange",
            "ztfg": "green",
            "ztfr": "red",
            "ztfi": "orange",
        }

        for fil in self.filters:
            assert fil in list(self.plot_colors.keys())

        if not light_curve.photopoints:
            self.logger.error("LightCurve has no photopoint")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        if not t2_views:  # Should not happen actually, T2Processor catches that case
            self.logger.error("Missing tied t2 views")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        excess_region: dict = {
            "baseline_jd": [],  # the jd of the last datapoint of the baseline before peak excess region
            "start_baseline_jd": [],  # the begining of the baseline before the peak excess region
            "baseline_mag": [],  # the mag of the last datapoint of the baseline before peak excess region
            "excess_jd": [],
            "excess_mag": [],
            "max_mag": [],
            "max_mag_jd": [],
            "significance": [],
            "strength_sjoert": [],
            "strength": [],
            "e_rise": [],
            "e_fade": [],
        }

        t2_output: dict[str, Any] = {"description": [], "values": []}

        if self.flux:
            intensity_low_limit = 0.0
            intensity_high_limit = 100000000000000000.0
        else:
            intensity_low_limit = 7.0
            intensity_high_limit = 15.0

        for t2_view in t2_views:
            self.filters_lc = self.filters
            if t2_view.unit == "T2BayesianBlocks":
                self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )

                for key in self.filters:
                    if key not in t2_res or t2_res.get(key) is None:
                        self.filters_lc.remove(key)

                if all(
                    (
                        t2_res[key]["nu_of_excess_regions"] == 0
                        and t2_res[key]["nu_of_baseline_regions"] <= 1
                    )
                    for key in self.filters_lc
                ):
                    t2_output["description"].append("Only baseline")
                    for key in self.filters_lc:
                        excess_region["max_mag"] = t2_res[key]["baseline"]
                elif t2_res["start_excess"] is not None or t2_res["size_excess"] > 1.0:
                    if t2_res["coincide_peak_block"] == -1:
                        t2_output["description"].append("No coincide excess regions")

                    elif not t2_output.get("description") and all(
                        min(t2_res[key]["max_mag_excess_region"]) < intensity_low_limit
                        or min(t2_res[key]["max_mag_excess_region"])
                        > intensity_high_limit
                        for key in self.filters_lc
                    ):
                        t2_output["description"].append("Very low or high magnitude")

                    elif all(
                        max(flatten(t2_res[key]["sigma_from_baseline"]))
                        < self.rej_sigma
                        for key in self.filters_lc
                    ):
                        t2_output["description"].append("Very low sigma")

                    elif any(
                        t2_res[key]["strength_sjoert"] < t2_res[key]["significance"]
                        for key in self.filters_lc
                    ):
                        t2_output["description"].append("Low significance")

                    if not t2_output.get("description"):
                        #  Check the excess region
                        #  Recalculate the excess region
                        maybe_interesting = False

                        for key in self.filters_lc:
                            if t2_res[key]["nu_of_excess_regions"] not in (0, None):
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

                                if (
                                    t2_res[key]["jd_excess_regions"][idx][-1]
                                    - t2_res[key]["jd_excess_regions"][idx][0]
                                    >= 1095.0
                                    and t2_res[key]["nu_of_excess_blocks"][idx] <= 1
                                ):
                                    t2_output["description"].append("Stage transition")

                                else:
                                    # fist check if there is a baseline before increase
                                    # Check the fluctuation inside the excess region
                                    #########################################################################
                                    # check if there is magnitude flactuanion before the excess region
                                    difference = [
                                        peak[0] - baseline[-1]
                                        for baseline in np.array(
                                            t2_res[key]["jd_baseline_regions"]
                                        )
                                    ]

                                    diff, position = min(
                                        (
                                            (b, nu)
                                            for nu, b in enumerate(np.array(difference))
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
                                                    t2_res[key]["jd_baseline_regions"]
                                                )
                                                if peak[0] - value[-1] > 0.0
                                            )
                                            > 500
                                        ):
                                            # Move the gap closer to excess region:
                                            baseline_jd = (
                                                t2_res[key]["jd_excess_regions"][idx][0]
                                                - 182.0
                                            )
                                            maybe_interesting = True
                                        else:
                                            baseline_jd = t2_res[key][
                                                "jd_baseline_regions"
                                            ][position][-1]

                                        baseline_mag = t2_res[key]["mag_edge_baseline"][
                                            position
                                        ]
                                        t2_output["description"].append(
                                            "Baseline before excess region"
                                        )

                                        excess_region["baseline_jd"].append(baseline_jd)
                                        start_baseline_jd = t2_res[key][
                                            "jd_baseline_regions"
                                        ][position]
                                        excess_region["start_baseline_jd"].append(
                                            start_baseline_jd[0]
                                        )

                                        excess_region["baseline_mag"].append(
                                            baseline_mag
                                        )
                                    #########################################################################
                                    # Check if there is magnitude fluctuations after the excess region
                                    difference = []
                                    for baseline in t2_res[key]["jd_baseline_regions"]:
                                        difference.append(peak[-1] - baseline[0])
                                    diff, position = min(
                                        (
                                            (b, nu)
                                            for nu, b in enumerate(np.array(difference))
                                            if b < 0
                                        ),
                                        default=(None, None),
                                    )
                                    if (
                                        diff is not None
                                        and t2_res[key]["description"][idx] == 0
                                    ):
                                        # Outlier in the middle of the LC, Outlier in the last epoch will not be marked as Outlier
                                        t2_output["description"].append("Outlier")
                                    elif (
                                        diff is None
                                        and t2_res[key]["description"][idx] == 0
                                    ):
                                        maybe_interesting = True
                                        t2_output["description"].append(
                                            "Excess region with one data point in the end of LC"
                                        )

                                    elif t2_res[key]["description"][idx] != 0:
                                        t2_output["description"].append(
                                            "Excess region exists"
                                        )

                                    excess_jd = t2_res[key]["jd_excess_regions"][idx][
                                        -1
                                    ]
                                    excess_mag = t2_res[key]["mag_edge_excess"][idx]
                                    excess_region["excess_jd"].append(excess_jd)
                                    excess_region["excess_mag"].append(excess_mag)
                            #################################################################################
                            else:
                                t2_output["description"].append("Only baseline")
                                excess_region["max_mag"].append(t2_res[key]["baseline"])
                    else:
                        t2_output["description"].append(
                            "Only baseline"
                        )  # no excess region

                        for key in self.filters_lc:
                            excess_region["max_mag"].append(t2_res[key]["baseline"])

        ##################################################################################
        if (
            len(excess_region.get("excess_jd", [])) >= 2
            and len(excess_region["baseline_jd"]) >= 2
            and not all(value == 0 for value in excess_region.get("baseline_jd", 0))
            and not (
                any(
                    "Outlier" in description for description in t2_output["description"]
                )
            )
        ):
            t2_output["values"] = excess_region

            for fid, key in enumerate(self.filters_lc, 1):
                if self.flux:
                    if self.data_type in ["ztf_alert", "ztf_fp"]:
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "flux_Jy", "flux_err_Jy"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": key,
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
                                    "value": key,
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
                                "value": key,
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
                        excess["jd"] >= excess_region["max_mag_jd"][fid - 1]
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

                    excess_region["significance"].append(t2_res[key]["significance"])
                    excess_region["strength_sjoert"].append(
                        t2_res[key]["strength_sjoert"]
                    )
                    excess_region["strength"].append(
                        t2_res[key]["max_sigma_excess_region"]
                    )

                    if self.flux:
                        excess_region["e_rise"].append(
                            (
                                excess_region["max_mag_jd"][fid - 1]
                                - excess_region["baseline_jd"][fid - 1]
                            )
                            / np.log(excess_region["max_mag"][fid - 1] / baseline)
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
                    excess_region["significance"].append(t2_res[key]["significance"])
                    excess_region["strength_sjoert"].append(
                        t2_res[key]["strength_sjoert"]
                    )
                    excess_region["strength"].append(
                        t2_res[key]["max_sigma_excess_region"]
                    )

                    excess_region["e_rise"].append("nan")
                    excess_region["e_fade"].append("nan")

            if all(
                value < 1000.0 for value in excess_region["e_rise"] if value != "nan"
            ) and all(
                value < 5000.0 and value > 400.0
                for value in excess_region["e_fade"]
                if value != "nan"
            ):
                if maybe_interesting is False:
                    if any(value == "nan" for value in excess_region["e_rise"]):
                        t2_output["status"] = "3"  #
                    else:
                        t2_output["status"] = "1"
                elif any(value == "nan" for value in excess_region["e_rise"]):
                    t2_output["status"] = "3_maybe_interesting"  #
                else:
                    t2_output["status"] = "1_maybe_interesting"
            elif maybe_interesting is False:
                if any(value == "nan" for value in excess_region["e_rise"]):
                    t2_output["status"] = "4"  #
                else:
                    t2_output["status"] = "2"  #
            elif any(value == "nan" for value in excess_region["e_rise"]):
                t2_output["status"] = "4_maybe_interesting"  #
            else:
                t2_output["status"] = "2_maybe_interesting"  #
        else:
            t2_output["status"] = "No further investigation"
            t2_output["values"] = excess_region

        #################### Plotting part #######################
        if t2_output["status"] != "No further investigation":
            fig, ax = plt.subplots(2, figsize=(18, 15))

            has_data = {}
            for fid, key in enumerate(self.filters_lc, 1):
                if self.flux:
                    if self.data_type in ["ztf_alert", "ztf_fp"]:
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "flux_Jy", "flux_err_Jy"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": key,
                                },
                            ],
                        )

                        if phot_tuple is None:
                            has_data[fid] = False
                            continue
                        has_data[fid] = True

                    elif self.data_type == "wise":
                        phot_tuple = light_curve.get_ntuples(
                            ["jd", "mean_flux", "flux_rms"],
                            [
                                {
                                    "attribute": "filter",
                                    "operator": "==",
                                    "value": key,
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
                                "value": key,
                            },
                            {"attribute": "mag_ul", "operator": "==", "value": "False"},
                        ],
                    )
                    has_data[fid] = True

                df = pd.DataFrame(phot_tuple, columns=["jd", "mag", "mag.err"])
                df["Outlier"] = False

                if t2_res[key]["jd_outlier"]:
                    for block in t2_res[key]["jd_outlier"]:
                        df.loc[
                            (df["jd"] >= block[0]) & (df["jd"] <= block[1]), "Outlier"
                        ] = True
                    outlier_datapoints = df[df["Outlier"]]
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

                datapoints = df[~df["Outlier"]]
                ax[0].errorbar(
                    datapoints["jd"] - 2400000.5,
                    datapoints["mag"],
                    yerr=datapoints["mag.err"],
                    label=key,
                    fmt="s",
                    color=self.plot_colors[key],
                    ecolor="lightgray",
                    markeredgecolor="black",
                    elinewidth=3,
                    capsize=0,
                )

            for fid, key in enumerate(self.filters_lc, 1):
                if has_data[fid] is False:
                    continue
                if self.flux:
                    idx = np.argmax(t2_res[key]["max_mag_excess_region"])
                else:
                    idx = np.argmin(t2_res[key]["max_mag_excess_region"])

                ax[0].axvline(
                    x=(t2_res[key]["jd_excess_regions"][idx][0] - 40) - 2400000.5,
                    color=self.plot_colors[key],
                )
                ax[0].axvline(
                    x=(t2_res[key]["jd_excess_regions"][idx][-1] + 40) - 2400000.5,
                    color=self.plot_colors[key],
                )
                ax[0].axvspan(
                    (t2_res[key]["jd_excess_regions"][idx][0] - 40) - 2400000.5,
                    (t2_res[key]["jd_excess_regions"][idx][-1] + 40) - 2400000.5,
                    label="T2 guess " + key,
                    alpha=0.05,
                    color=self.plot_colors[key],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"],
                    label="Baseline " + key,
                    ls="-",
                    color=self.plot_colors[key],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"] - t2_res[key]["baseline_sigma"],
                    label="Baseline sigma " + key,
                    ls="--",
                    color=self.plot_colors[key],
                )
                ax[0].axhline(
                    y=t2_res[key]["baseline"] + t2_res[key]["baseline_sigma"],
                    label="Baseline sigma " + key,
                    ls="--",
                    color=self.plot_colors[key],
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
            by_label = dict(zip(labels, handles, strict=False))
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

            pos = ([0.1, 0.8], [0.25, 0.8], [0.4, 0.8])
            ax[1].text(
                0.1, 0.8, "Final class: " + str(t2_output["status"]), fontsize=25
            )
            ax[1].text(
                pos[0][0],
                pos[0][-1] - 0.18,
                "Max sigma from baseline",
                fontsize=21,
            )

            ax[1].text(pos[0][0], pos[0][-1] - 0.26, "Max magnitude", fontsize=21)
            ax[1].text(pos[0][0], pos[0][-1] - 0.34, "Baseline", fontsize=21)
            ax[1].text(
                pos[0][0],
                pos[0][-1] - 0.42,
                "Time scale of excess region",
                fontsize=21,
            )
            ax[1].text(pos[0][0], pos[0][-1] - 0.50, "e-folding rise", fontsize=21)
            ax[1].text(pos[0][0], pos[0][-1] - 0.58, "e-folding fade", fontsize=21)

            for fid, key in enumerate(self.filters_lc, 1):
                if has_data[fid] is False:
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
                    pos[fid - 1][0] + 0.4,
                    pos[fid - 1][-1] - 0.2,
                    key + "\n\n",
                    color=self.plot_colors[key],
                    fontsize=25,
                )
                ax[1].text(
                    pos[fid - 1][0] + 0.4,
                    pos[0][-1] - 0.18,
                    format(max(flatten(t2_res[key]["sigma_from_baseline"])), ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    pos[fid - 1][0] + 0.4,
                    pos[0][-1] - 0.26,
                    format(t2_res[key]["max_mag_excess_region"][idx], ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    pos[fid - 1][0] + 0.4,
                    pos[0][-1] - 0.34,
                    format(t2_res[key]["baseline"], ".2f"),
                    fontsize=21,
                )
                ax[1].text(
                    pos[fid - 1][0] + 0.4,
                    pos[0][-1] - 0.42,
                    format(timescale, ".2f"),
                    fontsize=21,
                )
                if excess_region["e_rise"][fid - 1] != "nan":
                    ax[1].text(
                        pos[fid - 1][0] + 0.4,
                        pos[0][-1] - 0.50,
                        format(excess_region["e_rise"][fid - 1], ".2f"),
                        fontsize=21,
                    )
                    ax[1].text(
                        pos[fid - 1][0] + 0.4,
                        pos[0][-1] - 0.58,
                        format(excess_region["e_fade"][fid - 1], ".2f"),
                        fontsize=21,
                    )
                else:
                    ax[1].text(
                        pos[fid - 1][0] + 0.4,
                        pos[0][-1] - 0.50,
                        "NaN",
                        fontsize=21,
                    )
                    ax[1].text(
                        pos[fid - 1][0] + 0.4,
                        pos[0][-1] - 0.58,
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

        return t2_output
