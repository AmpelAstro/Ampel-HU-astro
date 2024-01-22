#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2RunPossis.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                11.12.2021
# Last Modified Date:  17.03.2022
# Last Modified By:    mf@physik.hu-berlin.de


import copy
from collections.abc import Sequence
from typing import Literal
from urllib.parse import urljoin
from urllib.request import urlopen

import numpy as np
import sncosmo  # type: ignore[import]
from sfdmap2.sfdmap import SFDMap  # type: ignore[import]

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2RunPossis(T2RunSncosmo):
    """
    Load a POSSIS kilnova model and fit to a LightCurve object as process is called.

    Load one of the POSSIS models and create an sncosmo_model
    model for fit by T2RunSncosmo.
    :param possis_base_url: str, path to (github) possis repository
    :param model_gen: str, name of model generation (subfolder)
    :mej_dyn: float, Possis parameter
    :mej_wind: float, Possis parameter
    :phi: int, Possis parameter
    :cos_theta: float, Possis parameter

    Dynamically fix model explosion time
    :param explosion_time_jd: Union[None, float, Literal['StockTriggerTime']]


    """

    # Parameters determining which POSSIS model will be read
    possis_base_url: str = "https://raw.githubusercontent.com/mbulla/kilonova_models/f810e0ec7e7a6ae32624738329e73d561b081372"
    model_gen: str = "bns_m3_3comp"
    mej_dyn: float = 0.01
    mej_wind: float = 0.09
    phi: int = 45
    cos_theta: float = 0.3  # Typically 0., 0.1, ...1.0

    sncosmo_model_name: str = "_".join(
        map(str, [model_gen, mej_dyn, mej_wind, phi, cos_theta])
    )

    possis_models: dict = {
        "model_ind": {
            "model_gen": model_gen,
            "mej_dyn": mej_dyn,
            "mej_wind": mej_wind,
            "phi": phi,
            "cos_theta": cos_theta,
        }
    }

    sncosmo_data: dict = {}

    # Fix time to specific explosion timestamp
    # StockTriggerTime assumes the value is updated during runtime
    explosion_time_jd: None | float | Literal["TriggerTime"]

    # Which units should this be changed to
    t2_dependency: Sequence[  # type: ignore[assignment]
        StateT2Dependency[
            Literal[
                "T2DigestRedshifts",
                "T2MatchBTS",
                "T2PhaseLimit",
                "T2PropagateStockInfo",
                "T2HealpixProb",
            ]
        ]
    ]

    def post_init(self) -> None:
        """
        Retrieve POSSIS models.

        Note that this could be done once at instance init.
        """

        for model_ind in self.possis_models:
            model_gen = self.possis_models[model_ind]["model_gen"]
            self.possis_models[model_ind]["sncosmo_model_name"] = "_".join(
                map(
                    str,
                    [
                        model_gen,
                        self.possis_models[model_ind]["mej_dyn"],
                        self.possis_models[model_ind]["mej_wind"],
                        self.possis_models[model_ind]["phi"],
                        self.possis_models[model_ind]["cos_theta"],
                    ],
                )
            )

        # print("T2RUNPOSSIS::", self.possis_models)

        for model_ind, model_dict in self.possis_models.items():
            model_gen = model_dict["model_gen"]
            mej_dyn = (
                model_dict["mej_dyn"] if model_dict.get("mej_dyn") else self.mej_dyn
            )
            mej_wind = (
                model_dict["mej_wind"] if model_dict.get("mej_wind") else self.mej_wind
            )
            phi = model_dict["phi"] if model_dict.get("phi") else self.phi
            cos_theta = (
                model_dict["cos_theta"]
                if model_dict.get("cos_theta")
                else self.cos_theta
            )
            sncosmo_model_name = (
                model_dict["sncosmo_model_name"]
                if model_dict.get("sncosmo_model_name")
                else self.sncosmo_model_name
            )

            apply_mwcorrection = (
                model_dict["apply_mwcorrection"]
                if model_dict.get("apply_mwcorrection")
                else self.apply_mwcorrection
            )

            model_url = urljoin(
                self.possis_base_url + "/",
                f"{model_gen}/nph1.0e+06_mejdyn{mej_dyn:05.3f}_mejwind{mej_wind:05.3f}_phi{phi}.txt",
            )

            # print(model_url)

            # Find file
            with urlopen(model_url) as fh:
                # Read model Parameters from first three lines
                lines = [next(fh).decode() for _ in range(3)]
                nobs = int(lines[0])
                nwave = int(lines[1])
                line3 = lines[2].split(" ")
                ntime = int(line3[0])
                t_i = float(line3[1])
                t_f = float(line3[2])
                model_cos_theta = np.linspace(0, 1, nobs)  # 11 viewing angles

                # print(cos_theta, model_cos_theta)
                phase = np.linspace(t_i, t_f, ntime)  # epochs

                # Limit to one angle
                # Note: U. Feindt developed model where angle was fit, left out for now
                if phi == 90:
                    theta_mask = np.array([True])
                else:
                    theta_mask = np.isclose(cos_theta, model_cos_theta)
                if not sum(theta_mask) == 1:
                    raise ValueError(f"Model cos_theta {model_cos_theta} not defined")

                # Read model data
                mdata = np.genfromtxt(fh)
            wave = mdata[0 : int(nwave), 0]
            flux = np.array(
                [
                    mdata[i * int(nwave) : i * int(nwave) + int(nwave), 1:]
                    for i in range(int(nobs))
                ]
            ).T

            # Reduce to one angle
            flux_1angle = flux[:, :, theta_mask].squeeze()
            # Create model
            source = sncosmo.TimeSeriesSource(
                phase, wave, flux_1angle, name=sncosmo_model_name
            )

            tmp_dustmap = None

            # Setup model, with or without MW correction
            if apply_mwcorrection:
                dust = sncosmo.models.CCM89Dust()
                tmp_sncosmo_model = sncosmo.Model(
                    source=source,
                    effects=[dust],
                    effect_names=["mw"],
                    effect_frames=["obs"],
                )
                tmp_dustmap = SFDMap()
                tmp_fit_params = copy.deepcopy(tmp_sncosmo_model.param_names)
                tmp_fit_params.remove("mwebv")
            else:
                tmp_sncosmo_model = sncosmo.Model(source=source)
                tmp_fit_params = copy.deepcopy(tmp_sncosmo_model.param_names)

            # If redshift _should_ be provided we remove this from fit parameters
            if self.redshift_kind is not None or self.backup_z is not None:
                tmp_fit_params.remove("z")

            # If explosion time should be fixed, do so
            # If explosion time should be fixed, do so
            if isinstance(self.explosion_time_jd, float):
                tmp_sncosmo_model.set(t0=self.explosion_time_jd)
                tmp_fit_params.remove("t0")

            self.sncosmo_data[model_ind] = {
                "sncosmo_model": tmp_sncosmo_model,
                "fit_params": tmp_fit_params,
                "dustmap": tmp_dustmap if tmp_dustmap else None,
                "apply_mwcorrection": apply_mwcorrection,
            }

        # self.default_param_vals = self.sncosmo_model.parameters

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """
        Load model and fit the loaded model to the data provided as a LightCurve.
        If requested, retrieve redshift and explosion time from t2_views.

        Return dict of dicts from model fitting.
        """

        result_dict = {}

        for model_ind in self.possis_models:
            # print("T2RUNPOSSIS evaluating::", model_ind)

            self.sncosmo_model = self.sncosmo_data[model_ind]["sncosmo_model"]
            self.fit_params = self.sncosmo_data[model_ind]["fit_params"]
            self.dustmap = self.sncosmo_data[model_ind]["dustmap"]
            self.apply_mwcorrection = self.sncosmo_data[model_ind]["apply_mwcorrection"]

            self.default_param_vals = self.sncosmo_model.parameters

            # Check if model explosion time should be fixed from t2
            if isinstance(self.explosion_time_jd, str):
                for t2_view in t2_views:
                    if t2_view.unit not in ["T2PropagateStockInfo", "T2HealpixProb"]:
                        continue
                    self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                    t2_res = (
                        res[-1]
                        if isinstance(res := t2_view.get_payload(), list)
                        else res
                    )
                    if "trigger_time" not in t2_res:
                        self.logger.info("No explosion time", extra={"t2res": t2_res})
                        return UnitResult(code=DocumentCode.T2_MISSING_INFO)
                    self.explosion_time_jd = float(t2_res["trigger_time"])
                    # Reset model
                    self.logger.debug(
                        "reset explosion time",
                        extra={"explosion_time": self.explosion_time_jd},
                    )
                    self.sncosmo_model.set(t0=self.explosion_time_jd)
                    self.fit_params.remove("t0")

            model_results = super().process(compound, datapoints, t2_views)

            result_dict[model_ind] = model_results

        # print("T2RUNPOSSIS::", result_dict)

        # Restart sncosmo processing
        # return super().process(compound, datapoints, t2_views)

        return result_dict
