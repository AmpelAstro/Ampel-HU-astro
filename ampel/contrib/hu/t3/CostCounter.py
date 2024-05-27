#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t3/CostCounter.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                8.1.2024
# Last Modified Date:  8.1.2024
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from collections.abc import Generator
from typing import TypedDict

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.struct.T3Store import T3Store
from ampel.types import T3Send, UBson
from ampel.view.SnapView import SnapView


class Costs(TypedDict):
    stocks: int
    states: int
    dps: int
    t2docs: dict[str, int]
    t2duration: dict[str, float]


class CostCounter(AbsT3Unit):
    """

    Derive metrics for the total cost, as parsed by the provided documents:
    - Number of states
    - dps associated with each state (and thus DB operations)
    - t2documents associated with each state
    - total compute time for each t2unit

    Potentially: number of alerts probed in session (T3SessionAlertsNumber)
    a proxy for number of alerts processed?

    To be complemented with total number of alerts that would have been processed
    if no initial DB query was made (only e.g. a test time range).

    """

    def process(
        self, gen: Generator[SnapView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson:
        """
        Loop through transients and add "costs".
        """

        costs: Costs = {
            "stocks": 0,
            "states": 0,
            "dps": 0,
            "t2docs": {},
            "t2duration": {},
        }

        # We will here loop through transients and react individually
        for view in gen:
            costs["stocks"] += 1
            costs["dps"] += len(view.t0) if view.t0 else 0
            costs["states"] += len(view.t1) if view.t1 else 0

            for t2v in view.t2 or []:
                unit = str(t2v.unit)
                if unit not in costs["t2docs"]:
                    costs["t2docs"][unit] = 0
                    costs["t2duration"][unit] = 0
                costs["t2docs"][unit] += 1

                # Add duration if any exists
                for metadict in t2v.meta:
                    costs["t2duration"][unit] += metadict.get("duration", 0)

        return costs  # type: ignore[return-value]
