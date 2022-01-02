#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/TransientInfoPrinter.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                11.06.2018
# Last Modified Date:  30.07.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

import logging
from typing import Any
from collections.abc import Generator
from ampel.types import UBson, T3Send
from ampel.struct.UnitResult import UnitResult
from ampel.view.T3Store import T3Store
from ampel.view.TransientView import TransientView
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.util.pretty import prettyjson


class TransientInfoPrinter(AbsPhotoT3Unit):

	logfile: None | str	# logging (INFO and onwards) goes to this file.

	def post_init(self) -> None:

		self.count = 0

		# log to file if so configured
		if self.logfile:
			fh = logging.FileHandler(self.logfile)
			fh.setLevel(logging.INFO)
			self.logger.addHandler(fh) # type: ignore
			self.logger.info("Added logging handle to: {logfile}")


	def process(self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None) -> UBson | UnitResult:

		self.logger.info("Printing transients info")
		self.logger.info("=" * 80)
		count = 0

		for tran_view in gen:
			count += 1
			TransientInfoPrinter._print_info(tran_view, self.logger)
			self.logger.info("=" * 80)

		self.logger.info(f"Printed info for {count} transients")

		return None


	@staticmethod
	def _print_info(tran: TransientView, logger: Any) -> None:

		logger.info(f"Ampel ID: {tran.id}") # type: ignore
		if tran.extra and 'name' in tran.extra:
			logger.info(f"Transient names: {tran.extra['name']}")

		if tran.stock:
			logger.info(f"Channel: {tran.stock['channel']}")

		created = tran.get_time_created()
		updated = tran.get_time_updated()

		logger.info("Created: %s" % (created if created is not None else 'Not available'))
		logger.info("Modified: %s" % (updated if updated is not None else 'Not available'))

		if tran.stock:
			logger.info(f"Tags: {tran.stock['tag']}")

		logger.info("Content: %s" % tran.content_summary())

		if tran.extra:
			logger.info(f"Extra: {tran.extra}")

		if tran.stock:
			# should be timely sorted
			logger.info("Journal:")
			for el in tran.stock['journal']:
				logger.info(prettyjson(el))

		if tran.logs:
			logger.info("Logs:")
			for le in tran.logs:
				logger.info(" " + str(le))
