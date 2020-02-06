#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/ChannelSummaryPublisher.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 13.11.2018
# Last Modified Date: 06.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import backoff, requests, datetime

from typing import List
from pytz import timezone
from astropy.time import Time
from io import StringIO, BytesIO
from typing import Optional, Dict, Any

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.utils.json import AmpelEncoder, load


class ChannelSummaryPublisher(AbsT3Unit):
	"""
	Create a json file with summary statistics for the channel. For the transients
	detected in the last N days, this json file contains, i.e. coords, RB score,
	first detection, latest detection, and the total number of transients detected
	by the channel.
	"""

	require = ('desycloud.default', )

	dryRun: bool = False
	baseDir: str = '/AMPEL/ZTF'
	alertMetrics: List[str] = [
		'ztf_name', 'ra', 'dec', 'rb', 'detections',
		'first_detection', 'last_detection'
	]


	def post_init(self, context: Optional[Dict[str, Any]]) -> None:

		self.summary = {}
		self._jd_range = [float('inf'), -float('inf')]
		self.dest = self.resource['desycloud.default'] + self.run_config.baseDir
		self.logger.info(
			f"Channel summary will include following metrics: {repr(self.alertMetrics)}"
		)
		self._channels = set()
		self.session = requests.Session()


	def extract_from_transient_view(self, tran_view):
		"""
		given transient view object return a dictionary
		with the desired metrics
		"""
		metrics = set(self.run_config.alertMetrics)
		out = {}
		out['ztf_name'] = tran_view.tran_names[0]
		out['tns_names'] = tran_view.tran_names[1:]

		# skip transient if it has no associated photopoints.
		if tran_view.photopoints is None or len(tran_view.photopoints) == 0:
			return
		# sort photopoints
		pps = sorted([pp.content for pp in tran_view.photopoints], key=lambda x: x['jd'])

		# some metric should only be computed for the latest pp
		for key in metrics:
			if key in ['ra', 'dec', 'rb']:
				out[key] = pps[-1].get(key)

		# look at full det history for start and end of detection
		if 'first_detection' in metrics:
			out['first_detection'] = pps[0]['jd']
		if 'last_detection' in metrics:
			out['last_detection'] = pps[-1]['jd']
		if 'detections' in metrics:
			out['detections'] = len(pps)
		if pps[-1]['jd'] < self._jd_range[0]:
			self._jd_range[0] = pps[-1]['jd']
		if pps[-1]['jd'] > self._jd_range[1]:
			self._jd_range[1] = pps[-1]['jd']

		return out


	def add(self, transients):
		"""
		load the stats from the alerts
		"""
		if transients is not None:
			for tran_view in transients:
				if isinstance(tran_view.channel, list):
					raise ValueError("Only single-channel views are supported")
				info_dict = self.extract_from_transient_view(tran_view)
				if not info_dict:
					continue
				key = info_dict.pop("ztf_name")
				self.summary[key] = info_dict
				self._channels.add(tran_view.channel)


	@backoff.on_exception(backoff.expo, (TimeoutError, requests.exceptions.HTTPError),
		giveup=lambda exc: isinstance(exc, requests.exceptions.HTTPError) and \
		exc.response.status_code not in {400, 403, 405, 423, 500}
	)
	def done(self):
		"""
		"""
		if len(self._channels) == 0:
			return
		elif len(self._channels) > 1:
			raise ValueError(f"Got multiple channels ({list(self._channels)}) in summary")

		# Find the date of the most recent observation, in Pacific time
		timestamp = Time(self._jd_range[-1], format='jd').to_datetime(timezone('US/Pacific'))
		# If before noon, it belongs to the night that started yesterday
		if timestamp.hour < 12:
			timestamp -= datetime.timedelta(days=1)
		filename = timestamp.strftime("channel-summary-%Y%m%d.json")

		channel = list(self._channels)[0]
		basedir = '{}/{}'.format(self.dest, channel)
		rep = self.session.head(basedir)
		if not (rep.ok or self.run_config.dryRun):
			self.session.request('MKCOL', basedir).raise_for_status()
		try:
			rep = self.session.get('{}/{}'.format(basedir, filename))
			rep.raise_for_status()
			partial_summary = next(load(BytesIO(rep.content)))
			partial_summary.update(self.summary)
			self.summary = partial_summary
		except (requests.exceptions.HTTPError, StopIteration):
			pass

		outfile = StringIO()
		outfile.write(AmpelEncoder(lossy=True).encode(self.summary))
		outfile.write('\n')
		mb = len(outfile.getvalue().encode()) / 2.0 ** 20
		self.logger.info("{}: {} transients {:.1f} MB".format(filename, len(self.summary), mb))
		if not self.run_config.dryRun:
			self.session.put('{}/{}'.format(basedir, filename), data=outfile.getvalue()).raise_for_status()
