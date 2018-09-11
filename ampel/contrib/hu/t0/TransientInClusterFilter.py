#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/XShooterFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 28.08.2018
# Last Modified Date: 28.08.2018
# Last Modified By  : m. giomi <matteo.giomi@desy.de>


import logging
from pymongo import MongoClient
from urllib.parse import urlparse

from extcats import CatalogQuery
from extcats.catquery_utils import get_distances

from ampel.contrib.hu.t0.DecentFilter import DecentFilter



class TransientInClusterFilter(DecentFilter):
	"""
		Filter derived from the DecentFilter, in addition selecting candidates
		with position compatible with that of nearby galaxy clusters.
	"""

	# Static version info
	version = 1.0
	resources = ('extcats.reader',)
	
	
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
			'BIG_SEARCH_RADIUS_ARCMIN',					# conservative search radius around cluster position. Max in RASSEBCS is 16.a arcmin
			'CLUSTER_RADIUS_MULTIPLIER'					# if we want to enlarge the search region around each cluster.
			)
		for el in config_params:
			if el not in run_config:
				raise ValueError("Parameter %s missing, please check your channel config" % el)
			if run_config[el] is None:
				raise ValueError("Parameter %s is None, please check your channel config" % el)
			self.logger.info("Using %s=%s" % (el, run_config[el]))
		
		self.big_search_radius_arcmin			= run_config['BIG_SEARCH_RADIUS_ARCMIN']
		self.cl_rad_multiply					= run_config['CLUSTER_RADIUS_MULTIPLIER']
		
		# convert the 'big search radius' from arcmin to arcsecs. 
		# Take into account the multiplier as well
		self.big_rs_arcsec = self.big_search_radius_arcmin*60*self.cl_rad_multiply
		self.logger.info("Big search radius to match with RASSEBCS clusters: %.2f arcsec"%
			self.big_rs_arcsec)
		
		# init the catalog query objects
		db_client = MongoClient(base_config['extcats.reader'])
		self.rassebcs_query = CatalogQuery.CatalogQuery(
													"RASSEBCS",
													ra_key='RA',
													dec_key='DEC',
													logger=self.logger,
													dbclient=db_client)
	
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
														alert_ra,
														alert_dec, 
														rs_arcsec = self.big_rs_arcsec)
		if nearby_clusters is None:
			self.logger.debug(
				"rejected: no cluster from RASSEBCS within %.2f arcsec from alert position"% 
				(self.big_rs_arcsec))
			return None
		
		
		# B) for all the nearby clusters, compute their distances to the alert
		clust_dists_2_alert = get_distances(alert_ra, alert_dec, nearby_clusters, "RA", "DEC")
		clust_dists_2_alert = clust_dists_2_alert/60.	# extcats works in arcsec but ANGULAR_RADIUS is in arcmin
		
		# C) if, for any of the nearby clusters, the distance to the alert is smaller
		# than the cluster size, count this as a match
		if not any(nearby_clusters['ANGULAR_RADIUS']>clust_dists_2_alert):
			self.logger.debug(
				"rejected: distance to alert is larger than the cluster size for all the nearby clusters")
			for ii in range(len(nearby_clusters)):
				self.logger.debug("Angular radius: %.2f, distance to alert: %.2f (both arcmin)"%
					(nearby_clusters['ANGULAR_RADIUS'][ii], clust_dists_2_alert[ii]))
			return None
		
		for ii in range(len(nearby_clusters)):
			self.logger.debug("Angular radius: %.2f, distance to alert: %.2f (both arcmin)"%
				(nearby_clusters['ANGULAR_RADIUS'][ii], clust_dists_2_alert[ii]))
		
		# now run the decent filter (faster to do it afterwards ;-)
		return DecentFilter.apply(self, alert)

