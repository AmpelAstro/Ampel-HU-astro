#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t0/RandFilter.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                14.12.2017
# Last Modified Date:  24.11.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from random import uniform
from ampel.abstract.AbsAlertFilter import AbsAlertFilter


class RandFilter(AbsAlertFilter):
    passing_rate: float

    def post_init(self):
        self.logger.info(
            f"RandFilter initialized with passing rate {self.passing_rate}"
        )

    def process(self, alert) -> bool:
        return uniform(0, 1) < self.passing_rate
