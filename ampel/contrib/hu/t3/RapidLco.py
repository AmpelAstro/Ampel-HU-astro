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

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.base.dataclass.JournalUpdate import JournalUpdate
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.pipeline.config.AmpelConfig import AmpelConfig
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils
from ampel.contrib.hu.t3.RapidBase import RapidBase




class RapidLco(RapidBase):
	"""
		Submit LCO triggers for candidates passing criteria.		
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


		# LCO trigger info
		lco_api : str = 'bf81e4a699d2f65529753cb65f8cce176b47d897'
		# A dict of LCO API triggers to be sent for each SN that fulfills all criteria.
		# Assumed to have the following key content 'trigger_name':{'start_delay':X (days), 'end_delay':Y (days), 'api_form':Z}
		# Where the start and end delays define the allowed LCO time range and the api_form provides the request to be submitted.
		# The following keys of the api_form will be changed: name, target:name, target:ra, target:dec, windows:end, windows:start
		lco_payload : dict = { 'lco_u_rapid' : {'start_delay':0, 'end_delay':1,'api_form':								{"group_id":"ZTF_rapid_sample","proposal":"TOM2020A-009","ipp_value":1.05,"operator":"SINGLE","observation_type":"RAPID_RESPONSE","requests":[{"acceptability_threshold":90,"configurations":[{"type":"EXPOSE","instrument_type":"1M0-SCICAM-SINISTRO","instrument_configs":[{"bin_x":1,"bin_y":1,"exposure_count":"2","exposure_time":"750","mode":"full_frame","rotator_mode":"","extra_params":{"defocus":0},"optical_elements":{"filter":"up"}},{"bin_x":1,"bin_y":1,"exposure_count":"1","exposure_time":"60","mode":"full_frame","rotator_mode":"","extra_params":{"defocus":0},"optical_elements":{"filter":"gp"}}],"acquisition_config":{"mode":"OFF","extra_params":{}},"guiding_config":{"mode":"ON","optional":True,"extra_params":{}},"target":{"name":"ZTF_rapid_sample","type":"ICRS","ra":"x","dec":"y","proper_motion_ra":0,"proper_motion_dec":0,"epoch":2000,"parallax":0},"constraints":{"max_airmass":1.6,"min_lunar_distance":30}}],"windows":[{"end":"x"}],"location":{"telescope_class":"1m0"}}]} },
					'lco_u_queue' : {'start_delay':1, 'end_delay':3,'api_form':								{"group_id":"ZTF_rapid_follow","proposal":"TOM2020A-009","ipp_value":1.05,"operator":"SINGLE","observation_type":"NORMAL","requests":[{"acceptability_threshold":90,"configurations":[{"type":"EXPOSE","instrument_type":"1M0-SCICAM-SINISTRO","instrument_configs":[{"bin_x":1,"bin_y":1,"exposure_count":"2","exposure_time":"500","mode":"full_frame","rotator_mode":"","extra_params":{"defocus":0},"optical_elements":{"filter":"up"}},{"bin_x":1,"bin_y":1,"exposure_count":"1","exposure_time":"30","mode":"full_frame","rotator_mode":"","extra_params":{"defocus":0},"optical_elements":{"filter":"gp"}}],"acquisition_config":{"mode":"OFF","extra_params":{}},"guiding_config":{"mode":"ON","optional":True,"extra_params":{}},"target":{"name":"ZTF_rapid_follow","type":"ICRS","ra":"x","dec":"x","proper_motion_ra":0,"proper_motion_dec":0,"epoch":2000,"parallax":0},"constraints":{"max_airmass":1.6,"min_lunar_distance":30}}],"windows":[{"end":"x","start":"x"}],"location":{"telescope_class":"1m0"}}]} }
					}




		# Cuts based on T2 catalog redshifts
		require_catalogmatch : bool = True   # Require a redshift max from a T2 output
		redshift_catalogs	: List[str] = [] # List of catalog-like output to search for redshift
		max_redshift	: float	= 0.05	# maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_redshift	: float	= 0.001	# minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_dist	: float = 1.5	# arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
		max_dist	: float = 30 	# arcsec, maximum distance 
		max_kpc_dist	: float = 999 	# kpc, maximum distance
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
		self.name = "RapidLco"
		self.logger.info("Initialized T3 RapidLco instance %s"%self.name, extra={'base_config':self.base_config,'run_config':self.run_config})
		

	def react(self, tran_view, info):
		"""
			Send a trigger to the LCO. 			
		"""		

		transient_name = ZTFUtils.to_ztf_id(tran_view.tran_id)

		# Step through all LCO submit forms
		success = True # Will be set to false if any submit fails
		submitted = []
		responses = []		

		for submit_name, submit_info in self.run_config.lco_payload.items():
	
			# Create submit dictionary	
			react_dict = AmpelConfig.recursive_unfreeze(submit_info['api_form'])
									
			# Update with information
			react_dict['name'] = submit_name+'_'+transient_name
			react_dict['requests'][0]['configurations'][0]['target']['name'] = transient_name
			react_dict['requests'][0]['configurations'][0]['target']['ra'] = str(info['ra'])
			react_dict['requests'][0]['configurations'][0]['target']['dec'] = str(info['dec'])
		
			# Some keys are not necessarily there
			timenow = datetime.datetime.utcnow() 
			if 'start' in react_dict['requests'][0]['windows'][0].keys():
				dtime = datetime.timedelta(days=submit_info['start_delay'])
				react_dict['requests'][0]['windows'][0]['start'] = '%s'%( (timenow+dtime) )
			if 'end' in react_dict['requests'][0]['windows'][0].keys():
				dtime = datetime.timedelta(days=submit_info['end_delay'])
				react_dict['requests'][0]['windows'][0]['end'] = '%s'%( (timenow+dtime) )


			self.logger.debug('Starting LCO trigger',extra={'target':transient_name, 'react_dict':react_dict})
			# Make a test to validate
			testreply = requests.post(
				'https://observe.lco.global/api/requestgroups/validate/',
				headers={'Authorization': 'Token {}'.format(self.run_config.lco_api)},
				json= react_dict
			)

			# Abort if we have errors
			if len(testreply.json()['errors'])>0:
				self.logger.info('Validating LCO trigger fails for for %s'%(transient_name),extra={'target':transient_name, 'react_dict':react_dict})
				success = False
				continue

			# Submit full trigger
			response = requests.post(
				'https://observe.lco.global/api/requestgroups/',
				headers={'Authorization': 'Token {}'.format(self.run_config.lco_api)},
				json= react_dict
			)

			# Check whether this was successful
			try:
				response.raise_for_status()
			except requests.exceptions.HTTPError as exc:
				self.logger.info('Submit LCO fails for %s'%(transient_name),extra={'target':transient_name, 'react_dict':react_dict,'response':response.content})
				success = False

			# Look at the response
			self.logger.info('Submit LCO succeeds for %s'%(transient_name),extra={'target':transient_name, 'react_dict':react_dict,'response':response.json()})

			submitted.append( react_dict )
			responses.append( response.json() )



		# Document what we did
		jcontent = {'t3unit': self.name, 'reactDicts': submitted, 'success':success, 'lcoResponses':responses}
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

			# Add some more info for display
			matchinfo['LCO paths'] = []


			self.logger.info("Passed reaction threshold", extra={"tranId":tv.tran_id})

			# Ok, so we have a transient to react to
			if self.run_config.do_react:
				success, jup = self.react(tv, matchinfo)
				if not jup is None:
					journal_updates.append(jup)
					for response in jup.content['lcoResponses']:
						matchinfo['LCO paths'].append( 'https://observe.lco.global/userrequests/{}/'.format(response['id']) )
				if success:
					self.logger.info("React success", extra={"tranId":tv.tran_id,"success":success})
				else:
					self.logger.info("React failure", extra={"tranId":tv.tran_id,"success":success})
			else:
				success = False
				jup = None


			# Otherwise, test
			matchinfo['LCO trigger success'] = success
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
		
		

