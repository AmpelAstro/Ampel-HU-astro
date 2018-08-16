#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TransientInfoDumper.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 15.08.2018
# Last Modified Date: 15.08.2018
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

from ampel.base.TransientView import TransientView
from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.utils.json import AmpelEncoder, object_hook
import json
import requests
import uuid
from urllib.parse import urlencode
from io import StringIO

class TransientViewDumper(AbsT3Unit):
	"""
	"""

	version = 0.1
	resources = ('desycloud.default',)

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		self.logger = logger
		self.count = 0
		self.outfile = StringIO()
		self.filename = str(uuid.uuid1()) + '.json'
		self.dest = base_config['desycloud.default'] + '/AMPEL/' + self.filename
		self.public = "https://desycloud.desy.de/index.php/s/W6J4sm7H6aAk7aq/download?{}".format(urlencode({'files': self.filename}))
		# don't bother preserving immutable types
		self.encoder = AmpelEncoder(lossy=True)

	def add(self, transients):
		"""
		"""
		if transients is not None:

			batch_count = len(transients)
			self.count += batch_count

			for tran_view in transients:
				self.outfile.write(self.encoder.encode(tran_view))
				self.outfile.write("\n")


	def done(self):
		"""
		"""
		mb = len(self.outfile.getvalue().encode()) / 2.0 ** 20
		self.logger.info("{:.1f} MB of JSONy goodness".format(mb))
		self.logger.info("Total number of transient printed: %i" % self.count)
		requests.put(self.dest, data=self.outfile.getvalue()).raise_for_status()
		self.logger.info(self.public)
