#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2PolynomialFit.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 29.03.2024
# Last Modified Date: 29.03.2024
# Last Modified By  : jnordin@physik.hu-berlin.de

from collections.abc import Iterable
from typing import Any

import numpy as np
from astropy.table.column import Column

from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


class T2PolynomialFit(AbsStateT2Unit, AbsTabulatedT2Unit):
    """
    Fit polynomial to each band.
    Return coefficients and chi2/dof
    """

    order: int = 1

    def eval_polyfit(
        self, time: Column, flux: Column, flux_err: Column
    ) -> dict[str, Any]:
        """
        Fit and evaluate a polynomical fit to input data.
        """
        pfit = np.polynomial.Polynomial.fit(time, flux, self.order, w=1 / flux_err)
        dof = len(time) - self.order - 1

        return {
            "chi2dof": sum((flux - pfit(time)) ** 2 / flux_err**2) / dof,
            "dof": dof,
            "coef": list(pfit.convert().coef),
        }

    def process(
        self,
        compound: T1Document,
        datapoints: Iterable[DataPoint],
    ) -> UBson | UnitResult:
        """
        Executed for each transient state.
        """

        # Obtain photometric table
        flux_table = self.get_flux_table(datapoints)

        # Fit polynmial to each band
        fitinfo: dict[str, Any] = {}
        for band in set(flux_table["band"]):
            i = flux_table["band"] == band
            if sum(i) <= (self.order + 1):  # Not enough data in band
                continue
            fitinfo[band] = self.eval_polyfit(
                flux_table[i]["time"], flux_table[i]["flux"], flux_table[i]["fluxerr"]
            )

        # Average over bands
        if len(fitinfo) > 0:
            for o in range(self.order + 1):
                fitinfo[f"p{o}"] = np.mean(
                    [v["coef"][o] for k, v in fitinfo.items() if isinstance(v, dict)]
                )
            fitinfo["chi2dof"] = np.mean(
                [v["chi2dof"] for k, v in fitinfo.items() if isinstance(v, dict)]
            )
            fitinfo["order"] = self.order
            return fitinfo

        return {"order": self.order, "chi2dof": None}
