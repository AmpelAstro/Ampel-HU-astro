#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/PlotLightcurveSample.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 28.09.2021
# Last Modified Date: 28.09.2021
# Last Modified By  : jnordin@physik.hu-berlin.de

import os
import re
from collections.abc import Generator
from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.view.TransientView import TransientView


class PlotLightcurveSample(AbsPhotoT3Unit):
    """

    Unit plots results from lightcurve fitters (RunSncosmo, RunParsnip).
    Potentially differentiate between BTS classes, include peak magnitude and correlation strength.

    Filtering is done in the T3 selecter - could potentially include multiple channels.
    However, assumes each SN only appears once.


    Todo:
    Construct some  Hubble diagram (if "alpha/beta" are included).

    {
        'lc_model' : 'sugar',
        'param' : ['q0','q1','q2','q3','Av'],
        'bts_classes' : {
            'SN Ia' : ['SN Ia', 'SN Ia-91T'],
            'SN II' : ['SN II', 'SN IIP', 'SN IIb', 'SN IIb', 'SN IIn'],
            'SN Ibc' : ['SN Ic-BL', 'SN Ib/c', 'SN Ic', 'SN Ib', 'SN Ibn'],
            'SLSNL' : ['SLSN-I', 'SLSN-II'],
            'TDE' : ['TDE'],
        },
        'z_source' : 'AMPELz'
    }
    """

    # Which model config should be used? Will be matched to either
    # model (T2RunParsnip) or model_name (T2RunSncosmo)
    unit_name: str = "T2RunSncosmo"
    lc_model: str
    model_config: Optional[int]

    # Which parameters to plot (correlate)
    param: List[str]
    param_bounds: Optional[Dict[str, List[float]]]
    include_absmag: bool = False

    # Which BTS classes should be highlighted (with key as title)?
    bts_classes: Dict[str, List[str]] = {
        "SN Ia": ["SN Ia", "SN Ia-91T"],
        "SN II": ["SN II", "SN IIP", "SN IIb", "SN IIb", "SN IIn"],
        "SN Ibc": ["SN Ic-BL", "SN Ib/c", "SN Ic", "SN Ib", "SN Ibn"],
        "SLSNL": ["SLSN-I", "SLSN-II"],
        "TDE": ["TDE"],
    }
    name_filter: dict[str, str] = {"ZTF name": "ZTF", "TNS ID": "TNS"}

    ## Methods for filtering data
    # Limit results to one of the two redshift sources (BTS, AMPELz)
    z_source: Optional[str]
    # If an ampel_z_group is provided, only use values at or lower than this value
    max_ampelz_group: Optional[int]
    # Cut on a minimal number of df
    min_dof: Optional[int]
    # Cut on an allowed range of chisq / dof
    chidof_range: Optional[List[float]]

    # Output parametrs
    plot_dir: str  # Plots will be stored here. Can be made Optional w db save
    save_table: bool = False
    plot_suffix: str = "pdf"  # ex pdf / png

    def _get_t2results(self, tran_view, t2unit, config=None):
        """
        Stable method for getting the latest completed t2result for a unit.
        """

        t2results = []
        if t2docs := tran_view.get_t2_docs(unit_id=t2unit):
            for t2doc in list(t2docs):
                if not t2doc["status"] >= 0:
                    continue
                if config is not None:
                    if not t2doc["config"] == config:
                        continue
                try:
                    t2results.append(t2doc["body"][-1]["result"])
                except KeyError:
                    pass
        return t2results

    def _get_parsnip_results(self, tran_view):
        """
        Look for results from RunParsnip corresponding to the config, if present.
        """

        t2res = self._get_t2results(tran_view, "T2RunParsnip", self.model_config)
        t2res = [t2 for t2 in t2res if self.lc_model == t2["model"]]
        t2res = [t2 for t2 in t2res if "prediction" in t2.keys()]
        if self.z_source:
            t2res = [
                t2
                for t2 in t2res
                if t2["z_source"] and re.match(self.z_source, t2["z_source"])
            ]
        if not len(t2res) > 0:
            return None
        if len(t2res) > 1:
            self.logger.info(
                "Multiple T2 results, grabbing one",
                extra={"stock": tran_view.id, "t2res": t2res},
            )
        # Check fit quality
        if self.min_dof and self.min_dof > t2res[0]["prediction"]["model_dof"]:
            return None
        if self.chidof_range:
            chidof = (
                t2res[0]["prediction"]["model_chisq"] / self.min_dof
                > t2res[0]["prediction"]["model_dof"]
            )
            if chidof < self.chidof_range[0] or chidof > self.chidof_range[1]:
                return None
        # Now get the results
        params = {}
        for pname in self.param:
            params[pname] = t2res[0]["prediction"][pname]
            try:
                params[pname + "_err"] = t2res[0]["prediction"][pname + "_error"]
            except KeyError:
                pass
        if self.include_absmag:
            params["absmag"] = t2res[0]["prediction"]["luminosity"]
            params["absmag_err"] = t2res[0]["prediction"]["luminosity_error"]

        return params

    def _get_sncosmo_results(self, tran_view):
        """
        Look for results from RunSncosmo corresponding to the config, if present.
        """

        if t2res := tran_view.get_t2_views(unit=self.unit_name):
            t2res = [
                t2 for t2 in t2res if t2.body and "sncosmo_result" in t2.body[-1].keys()
            ]
            t2res = [t2 for t2 in t2res if t2.body[-1]["sncosmo_result"]["success"]]
            t2res = [t2 for t2 in t2res if self.lc_model == t2.body[-1]["model_name"]]
            if self.z_source:
                t2res = [
                    t2
                    for t2 in t2res
                    if re.match(self.z_source, t2.body[-1]["z_source"])
                ]
            if self.model_config:
                t2res = [t2 for t2 in t2res if self.model_config == t2.config]
            if not len(t2res) > 0:
                return None
            if len(t2res) > 1:
                # self.logger.info('Multiple T2 results, grabbing last', extra={'stock':tran_view.id, 't2res': t2res} )
                print("why multiple?")
                print(t2res)
            t2body = t2res[-1].body[-1]

            # Check fit quality
            if self.min_dof and self.min_dof > t2body["sncosmo_result"]["ndof"]:
                return None
            if self.chidof_range:
                chidof = (
                    t2body["sncosmo_result"]["chisq"] / t2body["sncosmo_result"]["ndof"]
                )
                if chidof < self.chidof_range[0] or chidof > self.chidof_range[1]:
                    return None

            # Now get the results
            params = {}
            for pname in self.param:
                ix = t2body["sncosmo_result"]["param_names"].index(pname)
                params[pname] = t2body["sncosmo_result"]["parameters"][
                    t2body["sncosmo_result"]["param_names"].index(pname)
                ]
                params[pname + "_err"] = t2body["sncosmo_result"]["errors"][pname]

            # For some reason not always showing up?
            if self.include_absmag:
                params["absmag"] = t2body["fit_metrics"]["restpeak_model_absmag_B"]
                # params['absmag_err'] = -2.5 / np.log(10) * df['x0_err'] / df['x0']
            return params
        else:
            return None  # Explicit, nothing found.

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """
        Loop through provided TransientViews and retrieve data to plot.
        """

        table_rows: list[dict[str, Any]] = []
        for k, tran_view in enumerate(gen, 1):
            sninfo: dict[str, UBson] = {}
            # Collect base information applying to all states
            # If here, add stock info (name, channel etcs)
            if names := (tran_view.stock or {}).get("name", []):
                for label, name_str in self.name_filter.items():
                    r = re.compile(name_str)
                    # While names will mostly be unique, it might not always be the case.
                    sninfo[label] = list(filter(r.match, names))  # type: ignore[arg-type]
                    # Avoid list when possible
                    if (
                        isinstance((item := sninfo[label]), (list, tuple))
                        and len(item) == 1
                    ):
                        sninfo[label] = item[0]
            sninfo["stock"] = tran_view.id
            assert tran_view.stock
            sninfo["channels"] = tran_view.stock.get("channel")  # type: ignore[assignment]

            # Get AmpelZ info
            ampelz, ampelz_group = None, 99
            if isinstance(
                t2res := tran_view.get_latest_t2_body(unit="T2DigestRedshifts"), dict
            ):
                if "ampel_z" in t2res.keys():
                    # In case of multiple entries, no way to distinguish these (would have to add config selection to get_t2results)
                    sninfo["ampel_z_group"] = t2res["group_z_nbr"]
                    sninfo["ampel_z"] = t2res["ampel_z"]
            if self.max_ampelz_group:
                if (
                    isinstance(ampel_z_group := sninfo.get("ampel_z_group"), int)
                    and self.max_ampelz_group < ampel_z_group
                ):
                    # Skip this one, z too uncertain
                    continue

            # Get BTS class if available
            if isinstance(
                t2res := tran_view.get_latest_t2_body(unit="T2MatchBTS"), dict
            ):
                # Get plot class from BTS matching (if any)
                bts_class = "Unknown"
                if "bts_type" in t2res.keys():
                    for bts_name, bts_subclasses in self.bts_classes.items():
                        if t2res["bts_type"] in bts_subclasses:
                            bts_class = bts_name
                sninfo["class_from_bts"] = bts_class

            # Finally, we collect fit information from either T2RunSncosmo or T2RunParsnip
            lcinfo = self._get_sncosmo_results(tran_view)
            # Have to reawaken parsnip before this can be tested
            # if not lcinfo:
            #    lcinfo = self._get_parsnip_results(tran_view)
            if lcinfo:
                # No fit info fulfilling requirements
                # continue
                sninfo.update(lcinfo)

            table_rows.append(sninfo)

        # Convert
        df = pd.DataFrame.from_dict(table_rows)

        # Cut down based on parameter bounds
        if self.param_bounds:
            for pname, pbound in self.param_bounds.items():
                df = df[(df[pname] >= pbound[0]) & (df[pname] <= pbound[1])]

        # Create main pairgrid plot
        plt.figure(1, figsize=(20, 20))
        g = sns.PairGrid(df, hue="class_from_bts", vars=self.param)
        g.map_lower(plt.scatter)
        g.map_upper(sns.kdeplot)
        g.map_diag(sns.kdeplot, lw=3, legend=True)
        g.add_legend()
        ppath = os.path.join(
            self.plot_dir, "t3plotlightcurvesample." + self.plot_suffix
        )
        plt.savefig(ppath, dpi=300)
        plt.clf()
        plt.close()

        # Compare properties with absolute mag, if chosen
        if self.include_absmag:
            plt.figure(2, figsize=(20, 10))
            g = sns.PairGrid(
                df,
                hue="class_from_bts",
                x_vars=["absmag"] + self.param,
                y_vars="absmag",
            )
            g.map_diag(sns.histplot, color=".3")
            g.map_offdiag(sns.scatterplot)
            g.add_legend()
            ppath = os.path.join(
                self.plot_dir, "t3plotlightcurvesample_absmag." + self.plot_suffix
            )
            plt.savefig(ppath, dpi=300)
            plt.clf()
            plt.close()

        if self.save_table:
            df.to_csv(os.path.join(self.plot_dir, "plot_lightcurvesample.csv"))

        return None
