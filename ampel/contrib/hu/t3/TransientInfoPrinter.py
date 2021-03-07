#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/TransientInfoPrinter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 11.06.2018
# Last Modified Date: 10.06.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import logging
from typing import Optional, Tuple, Any, Dict, Union
from ampel.type import T3AddResult
from ampel.view.TransientView import TransientView
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit

class TransientInfoPrinter(AbsPhotoT3Unit):

	logfile: Optional[str]	# logging (INFO and onwards) goes to this file.

	def post_init(self) -> None:

		self.count = 0

		# log to file if so configured
		if self.logfile:
			fh = logging.FileHandler(self.logfile)
			fh.setLevel(logging.INFO)
			self.logger.addHandler(fh)
			self.logger.info("Added logging handle to: {logfile}")

		if self.context:
			self.logger.info(f"Context: {self.context}")


	def add(self, transients: Tuple[TransientView, ...]) -> T3AddResult:

		if not transients:
			self.logger.info("No transient provided")
			return None

		self.logger.info(f"Printing info of {len(transients)} transient(s)")
		self.count += len(transients)

		for tran_view in transients:
			TransientInfoPrinter._print_info(tran_view, self.logger)

		return None


	def done(self) -> Optional[Union[bool, Dict[str, Any]]]:
		self.logger.info(f"Total number of transient printed: {self.count}")
		return None


	@staticmethod
	def _print_info(tran: TransientView, logger: Any):

		logger.info(f"Ampel ID: {tran.id}") # type: ignore
		if tran.extra and 'name' in tran.extra:
			logger.info(f"Transient names: {tran.extra['name']}")

		if tran.stock:
			logger.info(f"Channel: {tran.stock['channel']}")

		created = tran.get_time_created()
		modified = tran.get_time_modified()

		logger.info("Created: %s" % (created if created is not None else 'Not available'))
		logger.info("Modified: %s" % (modified if modified is not None else 'Not available'))

		if tran.stock:
			logger.info(f"Tags: {tran.stock['tag']}")

		logger.info("Content: %s" % TransientView.content_summary(tran))

		if tran.extra:
			logger.info(f"Extra: {tran.extra}")

		if tran.stock:
			# should be timely sorted
			logger.info("Journal:")
			for el in tran.stock['journal']:
				logger.info(str(el))
