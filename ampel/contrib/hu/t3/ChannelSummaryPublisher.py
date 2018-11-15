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
import time, datetime
from pytz import timezone
from urllib.parse import urlencode
from io import StringIO

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.utils.json import AmpelEncoder
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils


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
		self.summary = {}
		
		self.dest = base_config['desycloud.default'] + '/AMPEL/ZTF/'
		# pick up which alert keys you want to be part of the summary
		self.default_alert_metrics = ['ztf_name', 'ra', 'dec', 'rb', 'detections', 'first_detection', 'last_detection']
		if run_config is None:
			self.alert_metrics = self.default_alert_metrics
		else:
			self.alert_metrics = run_config.get('alert_metrics')
		self.logger.info("Channel summary will include following metrics: %s"%repr(self.alert_metrics))
		self._channels = set()
		self.session = requests.Session()

	def extract_from_transient_view(self, tran_view):
		"""
			given transient view object return a dictionary 
			with the desired metrics
		"""
		
		out = {}
		
		out['ztf_name'] = ZTFUtils.to_ztf_id(tran_view.tran_id)
		
		# sort photopoints
		if tran_view.photopoints is None or len(tran_view.photopoints) == 0:
			return out
		pps = sorted([pp.content for pp in tran_view.photopoints], key=lambda x: x['jd'])

		# some metric should only be computed for the latest pp
		for key in self.alert_metrics:
			if key in ['ra', 'dec', 'rb']:
				out[key] = pps[-1].get(key)
		
		# look at full det history for start and end of detection
		if 'first_detection' in self.alert_metrics:
			out['first_detection'] = pps[0]['jd']
		if 'last_detection' in self.alert_metrics:
			out['last_detection'] = pps[-1]['jd']
		if 'detections' in self.alert_metrics:
			out['detections'] = len(pps)
		
		return out


	def add(self, transients):
		"""
			load the stats from the alerts
		"""
		if transients is not None:
			for tran_view in transients:
				if isinstance(tran_view.channel, list):
					raise ValueError("Only single-channel views are supported")
				self._channels.add(tran_view.channel)
				info_dict = self.extract_from_transient_view(tran_view)
				key = info_dict.pop("ztf_name")
				self.summary[key] = info_dict


	def done(self):
		"""
		"""
		if len(self._channels) == 0:
			return
		elif len(self._channels) > 1:
			raise ValueError("Got multiple channels ({}) in summary".format(list(self._channels)))

		# The latest ZTF night is the yesterday, Pacific time
		timestamp = datetime.datetime.fromtimestamp(time.time(), timezone('US/Pacific')) - datetime.timedelta(days=1)
		filename = timestamp.strftime("channel-summary-%Y%m%d.json")

		channel = list(self._channels)[0]
		basedir = '{}/{}'.format(self.dest, channel)
		rep = self.session.head(basedir)
		if not rep.ok:
			self.session.request('MKCOL', basedir).raise_for_status()

		outfile = StringIO()
		outfile.write(AmpelEncoder(lossy=True).encode(self.summary))
		outfile.write('\n')
		mb = len(outfile.getvalue().encode()) / 2.0 ** 20
		self.logger.info("{:.1f} MB of JSONy goodness".format(mb))
		self.session.put('{}/{}'.format(basedir, filename), data=outfile.getvalue()).raise_for_status()
