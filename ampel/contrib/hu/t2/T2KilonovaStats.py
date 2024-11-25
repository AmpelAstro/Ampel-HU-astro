import json
import os
from collections.abc import Mapping, Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy import stats

from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.util import get_payload
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


def jsonify(string: str):
    return (
        string.replace("'", '"')
        .replace("None", "null")
        .replace("True", "true")
        .replace("False", "false")
        .replace("inf", "Infinity")
    )


def dejsonify(arg):
    # print("got:",arg)
    try:
        return float(arg)
    except ValueError:
        c = {"-Infinity": -np.inf, "Infinity": np.inf, "NaN": np.nan}
        return c[arg]


class T2KilonovaStats(AbsTiedStateT2Unit):
    # Which units should this be chained to
    t2_dependency: Sequence[StateT2Dependency[Literal["T2KilonovaEval"]]]

    data_dir = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/kilonovaness"
    binned_gaus_dict: dict = {}

    def post_init(self):
        # print("T2KILONOVASTATS POST INIT")
        with open(
            os.path.join(self.data_dir, "binned_gaus_params.json"), encoding="utf-8"
        ) as file:
            self.binned_gaus_dict = json.load(file, parse_constant=dejsonify)
            file.close()

        print(self.binned_gaus_dict)

    def get_keys(self, map_dist=np.inf, ndet=np.inf):
        ndet_key = np.inf
        map_key = np.inf
        if ndet:
            for key in self.binned_gaus_dict:
                fkey = dejsonify(key)
                # print("NDET key:: ", key, " is made ", fkey, " ", type(fkey), " ndet = ", ndet, " ", type(ndet))
                if ndet <= fkey:
                    ndet_key = key
                    break

        if map_dist:
            for key in self.binned_gaus_dict[ndet_key]:
                fkey = dejsonify(key)
                # print(fkey)
                if map_dist <= fkey:
                    map_key = key
                    break

        return map_key, ndet_key

    def gaus_cdf(self, xval, loc: float = 0, scale: float = 1):
        return stats.norm.cdf(xval, loc=loc, scale=scale)

    def get_kn_stats(self, map_dist=np.inf, ndet=np.inf):
        map_key, ndet_key = self.get_keys(map_dist=map_dist, ndet=ndet)

        # print(ndet_key, map_key)
        return self.binned_gaus_dict[ndet_key][map_key]

    def get_kn_cumprob(self, kilonovaness, map_dist=np.inf, ndet=np.inf):
        gaus_stats = self.get_kn_stats(map_dist=map_dist, ndet=ndet)
        # print("GAUS STATS ", gaus_stats, " ", type(gaus_stats))

        return self.gaus_cdf(kilonovaness, loc=gaus_stats[0], scale=gaus_stats[1])

    def get_kn_densitiy_stats(self, kilonovaness, map_area, map_dist=np.inf):
        map_key, _ = self.get_keys(map_dist=map_dist, ndet=np.inf)
        map_key = str(dejsonify(map_key)).replace(".0", "")

        data_dir = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/kilonovaness"
        file_base = os.path.join(data_dir, "densities")

        density_files = os.listdir(file_base)
        print()
        print("T2KILONVOANESSSTATS:: ", map_key, type(map_key), density_files)
        print()
        files = [dfile for dfile in density_files if "-" + str(map_key) in dfile]
        if files:
            file = files[0]
        else:
            print("T2KilonovaStats: No density file found for ", map_key, " Mpc.")
            return {}

        dist_range = file[: file.find("_")]

        tmp_df = pd.read_csv(os.path.join(file_base, file), sep=";")
        tmp_df["kilonovaness"] = tmp_df["kilonovaness"].astype(int)
        tmp_df.index = tmp_df["kilonovaness"]

        exp_kn = map_area / 1000 * tmp_df.iloc[kilonovaness]["rate-1000"]
        exp_kn_pls = map_area / 1000 * tmp_df.iloc[kilonovaness]["plus"]
        exp_kn_min = exp_kn - map_area / 1000 * tmp_df.iloc[kilonovaness]["minus"]

        # statement = f"For map of area = {map_area} squ.deg. at {dist_range} Mpc, expect {exp_kn:.1f} (+{exp_kn_pls:.1f} -{exp_kn_min:.1f}) candidates with kn = {kilonovaness}."
        # print(statement)

        return {
            "exp_kn": exp_kn,
            "exp_kn_pls": exp_kn_pls,
            "exp_kn_min": exp_kn_min,
            "dist_range": dist_range,
        }

    def get_kn_total_stats(self, kilonovaness, map_area, map_dist=np.inf, ndet=np.inf):
        gaus_stats = (
            1 - self.get_kn_cumprob(kilonovaness, map_dist=map_dist, ndet=ndet)
        ) * 100
        # print(
        #    f"Map at this distance with this ndet, {gaus_stats:.5f}% of events have kn >= {kilonovaness}"
        # )

        res = self.get_kn_densitiy_stats(kilonovaness, map_area, map_dist=map_dist)
        res["gaus_percent"] = gaus_stats
        return res

    def kilonovaness_statistics(
        self, t2res: Mapping[str, Any]
    ) -> None | dict[str, Any]:
        # print("KEYS:: ", t2res.keys())
        return self.get_kn_total_stats(
            kilonovaness=t2res["kilonovaness"],
            map_area=t2res["map_area"],
            map_dist=t2res["map_dist"],
        )

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        print("T2KILONOVASTATS:: process")
        for t2_view in t2_views:
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            t2_res = get_payload(t2_view)

            if t2_view.unit == "T2KilonovaEval":
                return self.kilonovaness_statistics(t2_res)

        return {}
