#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/NoFilter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 14.12.2017
# Last Modified Date: 08.03.2019
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from ampel.abstract.AbsT0AlertFilter import AbsT0AlertFilter

class NoFilter(AbsT0AlertFilter):
	
	def __init__(self, on_match_t2_units, base_config=None, run_config=None, logger=None):
		self.on_match_default_t2_units = on_match_t2_units
		self.logger = logger

	def apply(self, ampel_alert):
		return self.on_match_default_t2_units
