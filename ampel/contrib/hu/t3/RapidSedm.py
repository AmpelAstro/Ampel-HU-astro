#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File				: ampel/contrib/hu/t3/rapidBase
# License			: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date				: 05.08.2019
# Last Modified Date		: 05.08.2019
# Last Modified By	        : jnordin@physik.hu-berlin.de



import re
import logging
import numpy as np
import requests, datetime
from pydantic import BaseModel, BaseConfig
from typing import Dict, List

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.dataclass.JournalUpdate import JournalUpdate
from ampel.logging.AmpelLogger import AmpelLogger
from ampel.ztf.utils.ZTFUtils import ZTFUtils
from ampel.contrib.hu.t3.RapidBase import RapidBase





class RapidSedm(RapidBase):
	"""
		Select transients for rapid reactions. Intended as base class where the react method can be
                implemented as wished and a testreact method posts test reactions to Slack.
		
		This version reacts by setting a target for SEDM observatoins
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
		slack_token		: str = "***REMOVED***"
		slack_channel 		: str = "#ztf_auto"
		slack_username		: str = "AMPEL"


		# Base SEDM trigger info
		sedm_url : str =  'http://pharos.caltech.edu/request'
		sedm_payload : dict = {
					'obj_ra': None,
					'obj_dec': None,
					'obj_epoch': 2000,
					'obj_mag': None,
					'obj_name': None,
					'allocation': 20180319205302725,
					'status': "PENDING",
					'inidate': None,
					'enddate': None,
					'priority': 5.,
					'ifu': True,
					'ab': "n",
					'ifu_use_mag': "y",
					'rc': False,
					'rc_use_mag': "y",
					'do_r': "y",
					'r_exptime': 0,
					'r_repeats': 1,
					'do_g': "y",
					'g_exptime': 0,
					'g_repeats': 1,
					'do_i': "y",
					'i_exptime': 0,
					'i_repeats': 1,
					'do_u': "y",
					'u_exptime': 0,
					'u_repeats': 1,
					'maxairmass': 2.5,
					'min_moon_dist': 30,
					'user_id': 284
					}
		sedm_username : str =  'jnordin'
		sedm_password : str =  '***REMOVED***'


		# Cuts based on T2 catalog redshifts
		require_catalogmatch : bool = True   # Require a redshift max from a T2 output
		redshift_catalogs	: List[str] = [] # List of catalog-like output to search for redshift
		max_redshift	: float	= 0.05	# maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_redshift	: float	= 0.001	# minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_dist	: float = 1.2	# arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
		max_dist	: float = 30 	# arcsec, maximum distance 
		max_absmag	: float	= -13	# max abs mag through peak mag and redshift from catalog mach (require both)
		min_absmag	: float	= -17	# min abs mag through peak mag and redshift from catalog mach (require both)

		# Cut on alert properties
		min_ndet	: int	= 2		# A candidate need to have at least this many detections
		min_ndet_postul	: int	= 2		# and if it has this minimum nr of detection after the last significant (max_maglim) UL.
		max_age		: float = 1.5		# days, If a detection has an age older than this, skip (stars,age).
		min_age		: float = 0		# Min age of detection history
		min_peak_mag	: float	= 19.25	# range of peak magnitudes for submission
		max_peak_mag	: float = 16	#
		min_n_filters	: int	= 1		# Reported detections in at least this many filters
		min_gal_lat	: float = 14	# Minimal galactic latitide
		ssdistnr_max	: float = 1		# reject alert if ssdistnr smaller than this value for any pp
		ps1_sgveto_rad	: float = 1		# reject alert if PS1 star for any pp
		ps1_sgveto_sgth	: float = 0.8	# 
		rb_minmed	: float = 0.3	# Minimal median RB.
		drb_minmed	: float = 0.995	# Minimal median RB.
		min_magrise     : float = -20   # NOT IMPLEMENTED

		maglim_min	: float = 19.25	# Limiting magnitude to consider upper limits as 'significant'
		maglim_maxago	: float = 2.5	# A limiting magnitude max this time ago
		
		
		# Cut to apply to all the photopoints in the light curve.
		# This will affect most operations, i.e. evaluating the position, 
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
		self.name = "RapidSedm"
		self.logger.info("Initialized T3 RapidSedm instance %s"%self.name)
		self.logger.info("base_config: %s"%self.base_config)
		self.logger.info("run_config: %s"%self.run_config)
		

	def react(self, tran_view, info):
		"""
			Send a trigger to the SEDM. Note that we have no good way of investigating the queue at this time			
		"""		

		
		# Assemble required information. These *should* already be present in the default info provided by the info dict returned for a sucessfull accept_tview
		react_dict = {}
		react_dict.update(self.run_config.sedm_payload)
		react_dict['obj_ra'] = info['ra']
		react_dict['obj_dec'] = info['dec']
		react_dict['obj_mag'] = info['latest_mag']    # Assuming that the object is not declining?
		react_dict['obj_name'] = ZTFUtils.to_ztf_id(tran_view.tran_id)
		react_dict['inidate'] = datetime.datetime.utcnow()
		react_dict['enddate'] = datetime.datetime.utcnow()+ datetime.timedelta(days=2)

		# We are still in debug stage, turn down priority
		#react_dict['priority'] = 1

		self.logger.debug('SEDM trigger for %s w dict %s'%(ZTFUtils.to_ztf_id(tran_view.tran_id), react_dict))

		# Make the post
		response = requests.post(self.run_config.sedm_url, data=react_dict,  auth=(self.run_config.sedm_username, self.run_config.sedm_password))
		# Check result
		if not response.status_code==200:
			success = False
		else:
			success = True


		# Document what we did
		jcontent = {'t3unit': self.name, 'reactDict': react_dict, 'success':success}
		jup = JournalUpdate(tran_id=tran_view.tran_id, ext=self.run_config.ext_journal, content=jcontent)
	

		return success, jup



		

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

			self.logger.info("Passed reaction threshold", extra={"tranId":tv.tran_id})

			# Ok, so we have a transient to react to
			if self.run_config.do_react:
				success, jup = self.react(tv, matchinfo)
				if not jup is None:
					journal_updates.append(jup)
				if success:
					self.logger.info("React success", extra={"tranId":tv.tran_id,"success":success})
				else:
					self.logger.info("React failure", extra={"tranId":tv.tran_id,"success":success})
			else:
				success = False
				jup = None


			# Otherwise, test
			matchinfo['SEDM trigger success'] = success 
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
		
		

