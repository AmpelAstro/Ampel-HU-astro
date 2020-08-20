#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/LensedTransientFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 04.27.2018
# Last Modified Date: 05.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from extcats import CatalogQuery
from pymongo import MongoClient
from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.alert.PhotoAlert import PhotoAlert


class LensedTransientFilter(AbsAlertFilter[PhotoAlert]):

	require = ('extcats.reader', )

	min_ndet: int
	ClusListSearchRadius: float
	MasterlensSearchRadius: float
	CaslteQSOSearchRadius: float


	def post_init(self):
		"""
		This filter reject candidates if they have less than a certain number
		of detection or if they are not positive subtractions (reference lower than sci),
		or if they do not match with the position of kwokn lenses.
		"""

		self.search_radiuses = {
			'cluslist': self.ClusListSearchRadius,
			'masterlens': self.MasterlensSearchRadius,
			'castleqso': self.CaslteQSOSearchRadius
		}

		# init the catalog query objects for the different lens catalogs
		catq_kwargs = {
			'logger': self.logger,
			'dbclient': MongoClient(self.resource['extcats.reader'])
		}
		cluslist_query = CatalogQuery.CatalogQuery(
			"cluslist", ra_key = 'ra_deg', dec_key = 'dec_deg', **catq_kwargs
		)
		mlens_query = CatalogQuery.CatalogQuery(
			"masterlens", ra_key = 'ra_coord', dec_key = 'dec_coord', **catq_kwargs
		)
		castle_query = CatalogQuery.CatalogQuery(
			"castleqso", ra_key = 'ra', dec_key = 'dec', **catq_kwargs
		)

		# group the catalogs together
		self.catqueries = {
			'cluslist': cluslist_query,
			'masterlens': mlens_query,
			'castleqso': castle_query
		}

		# Feedback
		for cat, rs in self.search_radiuses.items():
			self.logger.info("Catalog: %s --> Search radius: %.2e arcsec" % (cat, rs))


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
		for cat, catquery in self.catqueries.items():
			rs = self.search_radiuses[cat]
			if catquery.binaryserach(latest["ra"], latest["dec"], rs):
				self.logger.debug("searching matches in %s within %.2f arcsec" % (cat, rs))
				return True

		self.logger.debug("rejected: alert position did not match any lens in the catalogs.")

		return None
