#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/RandFilter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 14.12.2017
# Last Modified Date: 08.03.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from random import uniform
from typing import Any
from ampel.abstract.AbsAlertFilter import AbsAlertFilter


class RandFilter(AbsAlertFilter[Any]):

	version: float = 1.0
	passing_rate: float


	def post_init(self):
		self.logger.info(f"RandFilter initialized with passing rate {self.passing_rate}")


	def apply(self, alert: Any) -> bool:

		if uniform(0, 1) < self.passing_rate:
			return True

		return False
