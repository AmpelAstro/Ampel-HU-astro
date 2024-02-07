#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ScoreSingleObject
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                8.1.2024
# Last Modified Date:  8.1.2024
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from typing import Any

from ampel.contrib.hu.t3.AbsScoreCalculator import AbsScoreCalculator


class ScoreSingleObject(AbsScoreCalculator):
    """

    Calculate score based on how early a specific SN is detected.
    Based on data from T2InfantEval.

    """

    # Identifying the SN, done through coordinate matching
    ra: float = 210.910674637
    dec: float = 54.3116510708
    maxdist: float = 2.0 / 60 / 60  # max squared dist, in deg.

    # Calculating the score
    tzero: float = 2460083.6538773
    powerscale: float = -2.0  # timediff to tzero to this power gives scaling.

    t2scores_from: list[str] = ["T2InfantCatalogEval"]

    def evaluate(self, t2unit: str, t2_result: dict[str, Any]) -> float:
        """
        Return score if the matched transient correspond to a specific transient, with
        value given by time differents to provided tzero.
        """

        print(t2unit)
        print(t2_result)

        if t2unit != "T2InfantCatalogEval":
            return 0

        # Check if this is the desired object
        cdist = (t2_result["ra"] - self.ra) ** 2 + (t2_result["dec"] - self.dec) ** 2
        print(cdist)
        if cdist > self.maxdist**2:
            return 0

        print("score1")

        # Return score
        return (t2_result["t_max"] - self.tzero) ** self.powerscale
