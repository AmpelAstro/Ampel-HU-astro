#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/ChannelSummaryPublisher.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 13.11.2018
# Last Modified Date: 13.11.2018
# Last Modified By  : m. giomi <matteo.giomi@desy.de>

import json
import requests
import uuid
from urllib.parse import urlencode
from io import StringIO

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.utils.json import AmpelEncoder
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.pipeline.common.ZTFUtils import ZTFUtils


class ChannelSummaryPublisher(AbsT3Unit):
	"""
		Create a json file with summary statistics for the channel. For the transients
		detected in the last N days, this json file contains, i.e. coords, RB score,
		first detection, latest detection, and the total number of transients detected
		by the channel.
	"""

	version = 0.1
	resources = ('desycloud.default',)

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		self.logger = AmpelLogger.get_logger() if logger is None else logger
		self.count = 0
		self.outfile = StringIO()
		self.filename = "channel_summay_"str(uuid.uuid1()) + '.json'
		
		self.dest = base_config['desycloud.default'] + '/AMPEL/' + self.filename
		self.public = "https://desycloud.desy.de/index.php/s/W6J4sm7H6aAk7aq/download?{}".format(urlencode({'files': self.filename}))
		# don't bother preserving immutable types
		self.encoder = AmpelEncoder(lossy=True)
		
		# pick up which alert keys you want to be part of the summary
		self.default_alert_metrics = ['ztf_name', 'ra', 'dec', 'rb', 'first_detection', 'last_detection']
		if run_config is None:
			self.alert_metrics = self.default_alert_metrics
		else:
			self.alert_metrics = run_config.get('alert_metrics')
		self.logger.info("Channel summary will include following metrics: %s"%repr(self.alert_metrics))
		

	def extract_from_transient_view(tran_view):
		"""
			given transient view object return a dictionary 
			with the desired metrics
		"""
		
		out = {}
		
		if 'ztf_name' in self.alert_metric:
			out['ztf_name'] = ZTFUtils.to_ztf_id(tran_view.tran_id)
		
		# sort photopoints
		pps = sorted([pp.content for pp in tran_view.photopoints], key=lambda x: x['jd'])
		
		# some metric should only be computed for the latest pp
		for key in self.alert_metrics:
			if key in ['ra', 'dec', 'rb']:
				out[key] = pps[-1].get(key)
		
		# look at full det history for start and end of detection
		if 'first_detection' in self.alert_metrics:
			out['first_detection'] = pps[0]['jd']
		if 'first_detection' in self.alert_metrics:
			out['first_detection'] = pps[-1]['jd']
		
		return out


	def add(self, transients):
		"""
			load the stats from the alerts
		"""
		if transients is not None:
			for tran_view in transients:
				info_dict = self.extract_from_transient_view(tran_view)
				self.outfile.write(self.encoder.encode(info_dict))
				self.outfile.write("\n")
				


	def done(self):
		"""
		"""
		mb = len(self.outfile.getvalue().encode()) / 2.0 ** 20
		self.logger.info("{:.1f} MB of JSONy goodness".format(mb))
		self.logger.info("Total number of transient printed: %i" % self.count)
		requests.put(self.dest, data=self.outfile.getvalue()).raise_for_status()
		self.logger.info(self.public)
