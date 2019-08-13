#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TransientInfoPrinter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 11.06.2018
# Last Modified Date: 09.03.2019
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import logging
from pydantic import BaseModel
from ampel.base.TransientView import TransientView
from ampel.base.abstract.AbsT3Unit import AbsT3Unit

class TransientInfoPrinter(AbsT3Unit):
	"""
	"""
	
	class RunConfig(BaseModel):
		"""
		"""
		logfile		: str = None	# logging (INFO and onwards) goes to this file.
	
	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		self.logger = logger
		self.count = 0
		
		# add logging to file if so configured
		logfile = None if run_config is None else run_config.logfile
		if not logfile is None:
			fh = logging.FileHandler(logfile)
			fh.setLevel(logging.INFO)
			self.logger.addHandler(fh)
			self.logger.info("added logging handle to: %s"%logfile)
			
		if global_info is not None:
			self.logger.info(
				"Provided global info: %s" % 
				global_info
			)

	def add(self, transients):
		"""
		"""
		if transients is None:
			self.logger.info("No transient provided")
			return

		self.logger.info(
			"Printing info of %i transient(s)" % 
			len(transients)
		)
		print("")

		self.count += len(transients)

		ret=[]
		for tran_view in transients:
			TransientInfoPrinter._print_info(tran_view, self.logger)
			print("")

		return ret


	def done(self):
		"""
		"""
		self.logger.info(
			"Total number of transient printed: %i" % 
			self.count
		)


	@staticmethod
	def _print_info(tran, logger):
		""" 
		"""
		logger.info("Ampel ID: %i" % tran.tran_id)

		logger.info("Transient names: %s" % tran.tran_names)

		if tran.channel is not None:
			logger.info("Channel: %s" % str(tran.channel))

		created = tran.get_time_created(True)
		modified = tran.get_time_modified(True)

		logger.info("Created: %s" % (created if created is not None else 'Not available'))
		logger.info("Modified: %s" % (modified if modified is not None else 'Not available'))
		logger.info("Tags: %s" % (tran.tags if tran.tags is not None else "Not set"))
		logger.info("Latest state: %s" % (
			tran.latest_state.hex() if tran.latest_state is not None else "Not set"
		))
		logger.info("Content: %s" % TransientView.content_summary(tran))

		# should be timely sorted
		logger.info("Journal:")
		for el in tran.journal: 
			logger.info(str(el))
