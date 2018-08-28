#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/XShooterFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 28.08.2018
# Last Modified Date: 28.08.2018
# Last Modified By  : m. giomi <matteo.giomi@desy.de>


import logging
from numpy import array, logical_and
from urllib.parse import urlparse
from astropy.time import Time

from ampel.contrib.hu.t0.DecentFilter import DecentFilter



class XShooterFilter(DecentFilter):
	"""
		Filter derived from the DecentFilter, in addition selecting very new
		transients which are visible from the South. In particular, the transient
		are accepted if they are detected during the last 6h, at least one non-detection
		during the last 5 days (and no detection during this time).
	"""

	# Static version info
	version = 1.0
	resources = ('catsHTM.default',)
	
	
	def __init__(self, on_match_t2_units, base_config=None, run_config=None, logger=None):
		"""
		"""
		if run_config is None or len(run_config) == 0:
			raise ValueError("Please check you run configuration")

		self.on_match_t2_units = on_match_t2_units
		self.logger = logger if logger is not None else logging.getLogger()
		
		# init the parent DecentFilter
		DecentFilter.__init__(self, self.on_match_t2_units, base_config, run_config, logger=self.logger)
		
		# now add the parameters which are relevant for this
		# new filter. All the others are passed to the DecentFilter
		config_params = (
			'MAX_DEC',					# maximum allowed value for the declination
			'DET_WITHIN',				# the transient must have been detected within the last 'DET_WITHIN' days
			'UL_WITHIN'					# the transient must AT LEAST one ulim within the last 'UL_WITHIN' days
			)
		for el in config_params:
			if el not in run_config:
				raise ValueError("Parameter %s missing, please check your channel config" % el)
			if run_config[el] is None:
				raise ValueError("Parameter %s is None, please check your channel config" % el)
			self.logger.info("Using %s=%s" % (el, run_config[el]))
		
		# remember the pars
		self.max_dec 					= run_config['MAX_DEC']
		self.det_within					= run_config['DET_WITHIN']
		self.ul_within					= run_config['UL_WITHIN']
		self.keys_to_check += ('jd',)

	def apply(self, alert):
		"""
		run the decent filter on the alert
		"""
		
		# cut on declination
		latest = alert.pps[0]
		if latest['dec'] > self.max_dec:
			self.logger.debug("rejected: declination %.2f deg is above maximum allowed of %.2f deg"% 
				(latest['dec'], self.max_dec))
			return None
		
		# --------------------------------------------------------------------- #
		#					CUT ON THE HISTORY OF THE ALERT						#
		# --------------------------------------------------------------------- #
		
		now_jd = Time.now().jd
		self.logger.debug("Setting 'now' to JD %.4f to cut on alert history"%now_jd)
		
		# check on history 1: detected in the last 6h
		detection_jds = array(alert.get_values('jd', upper_limits=False))
		recent_detections = detection_jds > (now_jd - self.det_within)
		if not any(recent_detections):
			self.logger.debug("rejected: no detection within the last %.3f days (latest one %.3f days ago)"%
				(self.det_within, (now_jd - max(detection_jds))))
			return None
		
		# check on the history 2: at least one upper limit in the last 5 days
		ulim_jds = alert.get_values('jd', upper_limits=True)
		if ulim_jds is None:
			self.logger.debug("rejected: this alert has no upper limits")
			return None
		ulim_jds =array(ulim_jds)
		if not any(ulim_jds > (now_jd - self.ul_within)):
			self.logger.debug("rejected: no upper limit in the last %.3f days (latest one %.3f days ago)"%
			(self.ul_within, (now_jd - max(ulim_jds))))
			return None
		
		# check on history 3: no detection within the last 5 days
		not_so_recent_detections = detection_jds[~recent_detections]
		if any(not_so_recent_detections > (now_jd - self.ul_within)):
			self.logger.debug("rejected: at least one detection within the last %.3f days (but not within %.3f days)."%
				(self.ul_within, self.det_within))
			return None
		
		# now apply the DecentFilter
		return DecentFilter.apply(self, alert)
