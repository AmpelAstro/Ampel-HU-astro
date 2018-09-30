#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/RandFilter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 14.12.2017
# Last Modified Date: 04.07.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>


from ampel.base.abstract.AbsAlertFilter import AbsAlertFilter
from random import uniform

class RandFilter(AbsAlertFilter):


	version = 0.2


	def __init__(self, on_match_t2_units, base_config=None, run_config=None, logger=None):
		"""
		"""
		self.on_match_default_t2_units = on_match_t2_units

		if run_config is None:
			raise ValueError("run config required (threshold defined there)")

		self.passing_rate = run_config['passingRate']

		if logger is not None:
			logger.info("RandFilter with passing rate {}".format(self.passing_rate))
			self.logger = logger


	def apply(self, ampel_alert):
		"""
		"""

		rv = uniform(0,1)
		if rv < self.passing_rate:
			return self.on_match_default_t2_units
		else:
			self.logger.info("NO", extra={'rv': rv})
			return None
