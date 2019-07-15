#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File				: ampel/contrib/hu/t3/rapidBase
# License			: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date				: 15.07.2019
# Last Modified Date: 15.07.2019
# Last Modified By	        : jnordin@physik.hu-berlin.de

import re
import logging
import numpy as np
from pydantic import BaseModel, BaseConfig
from typing import Dict, List

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.base.dataclass.JournalUpdate import JournalUpdate
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils



# get the science records for the catalog match
def get_catalogmatch_srecs(tran_view, logger):
	cat_res = tran_view.get_science_records(t2_unit_id="CATALOGMATCH")
	if len(cat_res) == 0 or cat_res is None or cat_res[-1].get_results() is None:
		logger.info("NO CATALOG MATCH FOR THIS TRANSIENT")
		return {}
	return cat_res[-1].get_results()[-1]['output']


class RapidBase(AbsT3Unit):
	"""
		Select transients for rapid reactions. Intended as base class where the react method can be
                implemented as wished and a testreact method posts test reactions to Slack
	"""
	
	version = 0.1
	
	class RunConfig(BaseModel):
		"""
 		Necessary class to validate configuration.
		"""
		class Config(BaseConfig):
			"""
			Raise validation errors if extra fields are present
			"""
			allow_extra = False
			ignore_extra = False
		
		# Ampel config
		ext_journal			: bool	= True		# weather journal will go to separate collection.
		
		# React config
		do_react			: bool			# Unless set, no full reaction will be triggered

		# Test react config
		do_testreact			: bool			# If set, will post trigger to slack
		slack_token		: str = "xoxb-297846339667-549790069252-FLwKXkra0NL3FNnrvu9XYm4a"
		slack_channel 		: str = "#ampel_tns_test"
		slack_username		: str = "Ampel_TNS_test"


		# Cuts based on T2 catalog redshifts
		require_catalogmatch : bool = True   # Require a redshift max from a T2 output
		redshift_catalogs	: List[str] = [] # List of catalog-like output to search for redshift
		max_redshift	: float	= 0.05	# maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_redshift	: float	= 0.001	# minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_dist	: float = 1.2	# arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
		max_dist	: float = 50 	# arcsec, maximum distance 

		# Cut on alert properties
		min_ndet	: int	= 2		# A candidate need to have at least this many detections
		min_ndet_postul	: int	= 2		# and if it has this minimum nr of detection after the last significant (max_maglim) UL.
		max_age		: float = 2.5		# days, If a detection has an age older than this, skip (stars,age).
		min_age		: float = 0		# Min age of detection history
		min_peak_mag	: float	= 19.5	# range of peak magnitudes for submission
		max_peak_mag	: float = 16	#
		min_n_filters	: int	= 1		# Reported detections in at least this many filters
		min_gal_lat	: float = 14	# Minimal galactic latitide
		ssdistnr_max	: float = 1		# reject alert if ssdistnr smaller than this value for any pp
		ps1_sgveto_rad	: float = 1		# reject alert if PS1 star for any pp
		ps1_sgveto_sgth	: float = 0.8	# 
		rb_minmed	: float = 0.3	# Minimal median RB.
		drb_minmed	: float = 0.5	# Minimal median RB.

		maglim_min	: float = 19.5	# Limiting magnitude to consider upper limits as 'significant'
		maglim_maxago	: float = 2.5	# A limiting magnitude max this time ago
		
		
		# Cut to apply to all the photopoints in the light curve.
		# This will affect all operations, i.e. evaluating the position, 
		# computing number of detections ecc.
		lc_filters		: List[Dict]= [
										{
										'attribute': 'sharpnr',
										'operator': '>=', 
										'value': -10.15
										}, 
										{
										'attribute': 'magfromlim',
										'operator': '>',
										'value': 0
										}
									]
		

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		self.run_config = self.RunConfig() if run_config is None else run_config
		self.name = "RapidBase"
		self.logger.info("Initialized T3 RapidBase instance %s"%self.name)
		self.logger.info("base_config: %s"%self.base_config)
		self.logger.info("run_config: %s"%self.run_config)
		

	def react(self, tran_view, info):
		"""
			Replace with react method adopted to particular facility or output			
		"""		

		success = False		
		description = 'Tried to trigger X'

		# Replace with specific attempt, chance success & descrption
		raise NotImplementedError('No real reaction implemented in RapidBase')


		# Document what we did
		jcontent = {'t3unit': self.name, 'reaction': descrption, 'success':success}
		jup = JournalUpdate(tran_id=tran_view.tran_id, ext=self.run_config.ext_journal, content=jcontent)
	

		return success, jup



	def test_react(self, tran_view, info):
		"""
			Trigger a test slack report			
		"""		

		success = False		


		import datetime, io, json
		from slackclient import SlackClient
		from slackclient.exceptions import SlackClientError
		
		
		sc = SlackClient(self.run_config.slack_token)
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean", filters=self.run_config.lc_filters)
		msg = ("Pancha says: Look up %s at RA %s DEC %s. Added info %s"%
				(ztf_name, ra, dec, info) )
		api = sc.api_call(
			"chat.postMessage",
			channel = self.run_config.slack_channel,
			text = msg,
			username = self.run_config.slack_username,
			as_user = False
			)
		if not api['ok']:
			raise SlackClientError(api['error'])
		else:
			success = True

		description = 'Sent SLACK msg'
		self.logger.info(description, extra={'channel': self.run_config.slack_channel})


		# Document what we did
		jcontent = {'t3unit': self.name, 'reaction': description, 'success':success}
		jup = JournalUpdate(tran_id=tran_view.tran_id, ext=self.run_config.ext_journal, content=jcontent)

		return success, jup




	def accept_tview(self, tran_view):
		"""
			decide weather or not this transient is worth reacting to.
			
			NOTE that even if many of these cuts could defined passed directly to
			the task/job config, some of them still require relatively non trivial
			computation (e.g. 'age' of the transient). This makes this selection method 
			necessary.
		"""
				
		# We are lazy and create an info dict that can be included with the printout 
		# should properly not be part of accept method
		info = {}

		# get the latest light curve
		lc = tran_view.get_latest_lightcurve()
		
		# apply cut on history: consider photophoints which are sharp enough
		pps = lc.get_photopoints(filters=self.run_config.lc_filters)
		self.logger.info("%d photop. passed filter %s"%(len(pps), self.run_config.lc_filters))
		
		# cut on number of detection
		if len(pps) < self.run_config.min_ndet:
			self.logger.info("not enough detections: got %d, required %d"%
				(len(pps), self.run_config.min_ndet))
			return False
		info['detections'] = len(pps)

		# cut on age
		jds = [pp.get_value('jd') for pp in pps]
		most_recent_detection, first_detection = max(jds), min(jds)
		age = most_recent_detection - first_detection
		if age > self.run_config.max_age or age < self.run_config.min_age:
			self.logger.info("age of %.2f days outside of range [%.2f, %.2f]"%
				(age, self.run_config.min_age, self.run_config.max_age))
			return False
		info['age'] = age

		
		# cut on number of detection after last SIGNIFICANT UL
		ulims = lc.get_upperlimits(filters={'attribute':'diffmaglim', 'operator':'>=', 'value':self.run_config.maglim_min})
		if len(ulims) > 0:
			last_ulim_jd = sorted(ulims, key=lambda x: x.get_value('jd'))[-1].get_value('jd')
			pps_after_ndet = lc.get_photopoints( 
				filters = self.run_config.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_ulim_jd}])
			# Check if there are enough positive detection after the last significant UL
			if len(pps_after_ndet) < self.run_config.min_ndet_postul:
				self.logger.info("not enough consecutive detections after last significant UL.",
					extra={'NDet': len(pps), 'lastUlimJD': last_ulim_jd})
				return False
			# Check that there is a recent ul
			if (most_recent_detection-last_ulim_jd)>self.run_config.maglim_maxago:
				self.logger.info("No recent UL.",
					extra={'lastDet': most_recent_detection, 'lastUlimJD': last_ulim_jd})
				return False
			info['last_UL'] = most_recent_detection-last_ulim_jd
		else:
			self.logger.info("no UL")
			return False
		
		
		# cut on number of filters
		used_filters = set([pp.get_value('fid') for pp in pps])
		if len(used_filters) < self.run_config.min_n_filters:
			self.logger.info("requested detections in more than %d bands, got: %d"%
				(self.run_config.min_n_filters, len(used_filters)))
			return False
		
		# cut on range of peak magnitude
		mags = [pp.get_value('magpsf') for pp in pps]
		peak_mag = min(mags)
		if peak_mag > self.run_config.min_peak_mag or peak_mag < self.run_config.max_peak_mag:
			self.logger.info("peak magnitude of %.2f outside of range [%.2f, %.2f]"%
				(peak_mag, self.run_config.min_peak_mag, self.run_config.max_peak_mag))
			return False
		info['peakmag'] = peak_mag
		
		
		# cut on galactic coordinates
		ra, dec = lc.get_pos(ret="mean", filters=self.run_config.lc_filters)
		coordinates = SkyCoord(ra, dec, unit='deg')
		b = coordinates.galactic.b.deg
		if abs(b) < self.run_config.min_gal_lat:
			self.logger.info("transient at b=%.2f too close to galactic plane (cut at %.2f)"%
				(b, self.run_config.min_gal_lat))
			return False
		
		# cut on distance to closest solar system object
		# TODO: how to make this check: ('0.0' in list(phot["ssdistnr"])
		ssdist = np.array([pp.get_value('ssdistnr') for pp in pps])
		ssdist[ssdist==None] = -999
		#print (ssdist)

		close_to_sso = np.logical_and(ssdist < self.run_config.ssdistnr_max, ssdist > 0)
		if np.any(close_to_sso):
			self.logger.info("transient too close to solar system object", extra={'ssdistnr': ssdist.tolist()})
			return False
		
		# check PS1 sg for the full alert history
		distpsnr1, sgscore1 = zip(*lc.get_tuples('distpsnr1', 'sgscore1', filters=self.run_config.lc_filters))
		is_ps1_star = np.logical_and(
									np.array(distpsnr1) < self.run_config.ps1_sgveto_rad,
									np.array(sgscore1) > self.run_config.ps1_sgveto_sgth
								)
		if np.any(is_ps1_star):
			self.logger.info(
				"transient below PS1 SG cut for at least one pp.",
				extra={'distpsnr1': distpsnr1, 'sgscore1': sgscore1}
				)
			return False
		
		# cut on median RB and DRB score
		rbs = [pp.get_value('rb') for pp in pps]
		if np.median(rbs) < self.run_config.rb_minmed:
			self.logger.info(
				"Median RB %below limit.",
				extra={'median_rd': np.median(rbs), 'rb_minmed': self.run_config.rb_minmed}
				)
			return False
		drbs = [pp.get_value('drb') for pp in pps]
		if np.median(drbs) < self.run_config.drb_minmed:
			self.logger.info(
				"Median DRB %below limit.",
				extra={'median_drd': np.median(drbs), 'drb_minmed': self.run_config.drb_minmed}
				)
			return False
		


		# ----------------------------------------------------------------------#
		# 			  CUTS ON T2 RECORDS				#
		# ----------------------------------------------------------------------#


		#####      T2 Catalog matching
		cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)
		
		# check that we got any catalogmatching results (that it was run)
		if self.run_config.require_catalogmatch:

			if len(cat_res)==0:
				self.logger.info("no T2CATALOGMATCH results")
				return False

			# Loop through listed catalogs for match
			zmatch = False
			for catname in self.run_config.redshift_catalogs:
				catinfo = cat_res.get(catname, False)
				if ( catinfo and (self.run_config.min_redshift <  catinfo['z'] < self.run_config.max_redshift) 
					and (self.run_config.min_dist <  catinfo['dist2transient'] < self.run_config.max_dist) ) :
						self.logger.info("z matched.", extra= {'catalog': catname,'z':catinfo['z'],'dist':catinfo['dist2transient']})
						zmatch = True
						info[catname+'_z'] = catinfo['z']
						info[catname+'_dist2transient'] = catinfo['dist2transient']

			if not zmatch:
				self.logger.info('No z match.')
				return False
		
		# tag AGNs
		milliquas = cat_res.get('milliquas', False)
		sdss_spec = cat_res.get('SDSS_spec', False)
		if (milliquas and milliquas['redshift'] > 0):
			info['milliAGN'] = True
		if (sdss_spec and sdss_spec['bptclass'] in [4, 5]):
			info ['sdssAGN'] = True
		
	
		### Potentially other checks on T2 results, eg photoz and lightcurve

		
		# congratulation TransientView, you made it!
		return info

		
		



	def add(self, transients):
		"""
			Loop through transients and check for TNS names and/or candidates to submit
		"""
		
		if transients is None:
			self.logger.info("no transients for this task execution")
			return []

		journal_updates = []		
		# We will here loop through transients and react individually
		for tv in transients:
			matchinfo = self.accept_tview(tv)

			# Check sumission criteria
			if not matchinfo:
				continue

			self.logger.info("Passed reaction threshold", extra={"tranId":tran_view.tran_id})

			# Ok, so we have a transient to react to
			if self.run_config.do_react:
				success, jup = self.react(tv, matchinfo)
				if not jup is None:
					journal_updates.append(jup)


			# Otherwise, test
			if self.run_config.do_testreact:
				test_success, jup = self.test_react(tv,matchinfo)
				if not jup is None:
					journal_updates.append(jup)


		return journal_updates


	def done(self):
		"""
		"""

		# Should possibly do some accounting or verification
		self.logger.info("done running T3")
		
		

