#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TransientInfoPrinter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 11.06.2018
# Last Modified Date: 11.06.2018
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from ampel.base.TransientView import TransientView
from ampel.abstract.AbsT3Unit import AbsT3Unit


class TransientInfoPrinter(AbsT3Unit):
	"""
	"""

	version = 0.1


	def __init__(self, logger, base_config=None):
		"""
		"""
		self.logger = logger


	def run(self, run_config, transients=None):
		"""
		"""

		self.logger.info("run(...) called")

		for tran_view in transients:
			TransientView._print_info(tran_view, self.logger)

		self.logger.info("end of run(...)")
