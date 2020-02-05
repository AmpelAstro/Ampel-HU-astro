#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/RandFilter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 14.12.2017
# Last Modified Date: 03.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from random import uniform
from typing import Optional, Sequence, ClassVar
from ampel.alert.PhotoAlert import PhotoAlert
from ampel.abstract.AbsPhotoAlertFilter import AbsPhotoAlertFilter


class RandFilter(AbsPhotoAlertFilter):
	""" """

	version: ClassVar[float] = 1.0
	passing_rate: float


	def post_init(self):
		""" """
		self.logger.info(f"RandFilter with passing rate {self.passing_rate}")


	def apply(self, alert: PhotoAlert) -> Optional[Sequence[str]]:
		""" """

		if uniform(0, 1) < self.passing_rate:
			return self.on_match_t2_units

		return None
