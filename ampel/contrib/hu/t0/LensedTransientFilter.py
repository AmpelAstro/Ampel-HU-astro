#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/LensedTransientFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 04.27.2018
# Last Modified Date: 24.08.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

from functools import partial

from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.alert.PhotoAlert import PhotoAlert
from ampel.contrib.hu.base.ExtcatsUnit import ExtcatsUnit


class LensedTransientFilter(ExtcatsUnit, AbsAlertFilter[PhotoAlert]):

	min_ndet: int
	ClusListSearchRadius: float
	MasterlensSearchRadius: float
	CaslteQSOSearchRadius: float


	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		searches = {
			'cluslist': (self.ClusListSearchRadius, 'ra_deg', 'dec_deg'),
			'masterlens': (self.MasterlensSearchRadius, 'ra_coord', 'dec_coord'),
			'castleqso': (self.CaslteQSOSearchRadius, 'ra', 'dec'),
		}
		self.queries = {
			k: partial(
				self.get_extcats_query(
					k,
					ra_key=ra_key,
					dec_key=dec_key
				).binaryserach, # sic
				rs_arcsec=radius
			)
			for k, (radius, ra_key, dec_key) in searches.items()
		}


	def apply(self, alert):
		"""
		Mandatory implementation.
		To exclude the alert, return *None*
		To accept it, either
			* return self.on_match_default_flags
			* return a custom combination of T2 unit names
		"""

		# cut on the number of previous detections
		if len(alert.pps) < self.min_ndet:
			return None

		# now consider the last photopoint
		latest = alert.pps[0]

		# check if it a positive subtraction
		if not (
			latest['isdiffpos'] and
			(latest['isdiffpos'] == 't' or latest['isdiffpos'] == '1')
		):
			self.logger.debug("rejected: 'isdiffpos' is %s", latest['isdiffpos'])
			return None

		# and match with the catalogs using position of latest photopoint
		for cat, query in self.queries.items():
			if query(latest["ra"], latest["dec"]):
				return True

		return None
