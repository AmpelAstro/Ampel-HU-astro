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


	def __init__(self, logger, run_config=None, base_config=None):
		"""
		"""
		self.logger = logger
		self.count = 0


	def add(self, transients):
		"""
		"""

		if transients is not None and len(transients) > 0:
			self.logger.info("Adding new transient(s)")
			for tran_view in transients:
				TransientView._print_info(tran_view, self.logger)
				self.count += 1
		else:
			self.logger.info("No transient provided")


	def run(self):
		"""
		"""
		self.logger.info("Total number of transient printed: %i" % self.count)


