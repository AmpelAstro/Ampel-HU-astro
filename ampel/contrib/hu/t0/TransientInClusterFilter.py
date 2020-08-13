#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/TransientInClusterFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 28.08.2018
# Last Modified Date: 05.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from pymongo import MongoClient
from extcats.CatalogQuery import CatalogQuery
from extcats.catquery_utils import get_distances
from ampel.contrib.hu.t0.DecentFilter import DecentFilter


class TransientInClusterFilter(DecentFilter):
	"""
	Filter derived from the DecentFilter, in addition selecting candidates
	with position compatible with that of nearby galaxy clusters.
	"""

	require = 'extcats.reader', 'catsHTM.default'

	big_search_radius_arcmin: float	# conservative search radius around cluster position. Max in RASSEBCS is 16.a arcmin
	cluserter_rad_multiplier: float # if we want to enlarge the search region around each cluster.

	def post_init(self):

		super().post_init()

		# feedback
		for k in self.__annotations__:
			self.logger.info(f"Using {k}={getattr(self, k)}")

		# convert the 'big search radius' from arcmin to arcsecs.
		# Take into account the multiplier as well
		self.big_rs_arcsec = self.big_search_radius_arcmin * 60 * self.cluserter_rad_multiplier
		self.logger.info(
			f"Big search radius to match with RASSEBCS clusters: {self.big_rs_arcsec:.2f} arcsec"
		)

		# init the catalog query objects
		self.rassebcs_query = CatalogQuery(
			"RASSEBCS", ra_key = 'RA', dec_key = 'DEC', logger = self.logger,
			dbclient = MongoClient(self.resource['extcats.reader'])
		)


	def apply(self, alert):
		"""
		run the filter on the alert. First we run the decent filter, then we match
		with the cluster catalog.
		"""

		# if the candidate has passed the decent filter, check if it is compatible
		# with the position of some nearby galaxy cluster
		latest = alert.pps[0]
		alert_ra, alert_dec = latest['ra'], latest['dec']

		# A) find position of all the nearby clusters. If none is found, reject alert.
		nearby_clusters = self.rassebcs_query.findwithin(
			alert_ra, alert_dec, rs_arcsec = self.big_rs_arcsec
		)

		if nearby_clusters is None:
			self.logger.debug(
				f"rejected: no cluster from RASSEBCS within {self.big_rs_arcsec:.2f} arcsec from alert position"
			)
			return None

		# B) for all the nearby clusters, compute their distances to the alert
		# (extcats works in arcsec but ANGULAR_RADIUS is in arcmin)
		clust_dists_2_alert = get_distances(alert_ra, alert_dec, nearby_clusters, "RA", "DEC") / 60.

		# C) if, for any of the nearby clusters, the distance to the alert is smaller
		# than the cluster size, count this as a match
		if not any(nearby_clusters['ANGULAR_RADIUS'] > clust_dists_2_alert):

			self.logger.debug(
				"rejected: distance to alert is larger than the cluster size for all the nearby clusters"
			)

			for ii in range(len(nearby_clusters)):
				self.logger.debug(
					f"Angular radius: {nearby_clusters['ANGULAR_RADIUS'][ii]:.2f}, "
					f"distance to alert: {clust_dists_2_alert[ii]:.2f} (both arcmin)"
				)

			return None

		for ii in range(len(nearby_clusters)):
			self.logger.debug(
				f"Angular radius: {nearby_clusters['ANGULAR_RADIUS'][ii]:.2f}, "
				f"distance to alert: {clust_dists_2_alert[ii]:.2f} (both arcmin)"
			)

		# now run the decent filter (faster to do it afterwards ;-)
		return super().apply(alert)
