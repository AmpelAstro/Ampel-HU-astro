#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t3/AbsScoreCalculator.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                8.1.2024
# Last Modified Date:  8.1.2024
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from collections.abc import Generator
from typing import Any, Literal

from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit
from ampel.base.decorator import abstractmethod
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.view.SnapView import SnapView


class AbsScoreCalculator(AbsT3ReviewUnit, abstract=True):
    """

    Abstract class for units which evaluated the final performance of a pipepline.
    Each requested t2 output will be fed into an _evaluate_ class, which has to be
    implemented by the child.

    The intended use case is for workflow optimization, where variants of the same setup can be rerun
    with varying parameters and the the score values compared. This can, in turn, be compared with the
    resource required for each runthrough.

    The summed score for all transients are returned, individully for each unit as well as
    for all of them together. Each specific application implements _evaluate_ to allow nonlinear behaviour, combined analysis of
    multiple t2 outputs, and cases where a specific object is looked for.
    """

    # List of T2 unit names which should be fed into the evaluate method
    t2scores_from: list[str] = []

    # In case multiple results exists for the same t2 unit, should the sum, max or min score be used?
    scores_parsing: Literal["sum", "max", "min"] = "max"

    def post_init(self) -> None:
        self.t2scores = {t2name: 0.0 for t2name in self.t2scores_from}

    def process(
        self, gen: Generator[SnapView, T3Send, None], t3s: None | T3Store = None
    ) -> UBson | UnitResult:
        """
        Loop through transients and calculate scores.
        """

        # We will here loop through transients and react individually
        for view in gen:
            for t2unit in self.t2scores_from:
                # Need to take into account possibly multiple t2 results for each kind of unit.
                scores = [
                    self.evaluate(t2unit, t2_result.body[-1])
                    for t2_result in view.get_t2_views(unit=t2unit)
                    if t2_result.body and isinstance(t2_result.body[-1], dict)
                ]
                if len(scores) == 0:
                    continue

                if self.scores_parsing == "sum":
                    self.t2scores[t2unit] += sum(scores)
                elif self.scores_parsing == "max":
                    self.t2scores[t2unit] += max(scores)
                elif self.scores_parsing == "min":
                    self.t2scores[t2unit] += min(scores)

        return {"scores_unit": self.t2scores, "scores_sum": sum(self.t2scores.values())}

    @abstractmethod
    def evaluate(self, t2unit: str, t2_result: dict[str, Any]) -> float:
        """
        Replace with evaluate method for workflow test
        """
        raise NotImplementedError
