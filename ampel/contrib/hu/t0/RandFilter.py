#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/RandFilter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 14.12.2017
# Last Modified Date: 14.11.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import logging
from random import uniform
from typing import Any, Dict
from pydantic import BaseModel
from ampel.abstract.AbsPhotoAlertFilter import AbsPhotoAlertFilter
#from ampel.model.t0.T0FilterModel import T0FilterModel


class RandFilter(AbsPhotoAlertFilter):

	class InitConfig(AbsPhotoAlertFilter.InitConfig):
	#class InitConfig(T0FilterModel):
		""" """
		passing_rate: float


	# pylint: disable=super-init-not-called
	def __init__(
		self, logger: logging.Logger, init_config: BaseModel = None, 
		resources: Dict[str, Any] = None
	):
		"""
		"""
		self.on_match_t2_units = init_config.on_match_t2_units

		if init_config is None:
			raise ValueError("Init config is required (threshold defined there)")

		self.passing_rate = init_config.passing_rate

		if logger is not None:
			logger.info(f"RandFilter with passing rate {self.passing_rate}")
			self.logger = logger


	def apply(self, alert):
		"""
		"""

		rv = uniform(0,1)
		if rv < self.passing_rate:
			return self.on_match_t2_units
		else:
			return None
