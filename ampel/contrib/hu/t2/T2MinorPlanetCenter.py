#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2LCQuality.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 10.01.2019
# Last Modified Date: 10.01.2019
# Last Modified By  : jnordin@physik.hu-berlin.de


import logging
logging.basicConfig()

#from astropy.table import Table, unique
import numpy as np


from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils
from ampel.base.abstract.AbsT2Unit import AbsT2Unit
from ampel.core.flags.T2RunStates import T2RunStates
from pydantic import BaseModel, BaseConfig
from astropy.time import Time

class T2MinorPlanetCenter(AbsT2Unit):
	"""
		Check if the *latest* detection of a transient corresponds
		matches something known by the MinorPlanetCenter.	
		
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
		only_latest		: bool	= True		# Will only match the latest photopoint.

		filters			: dict	= None		# Potential filter for photopoint selection



	def __init__(self, logger, base_config):
		"""
		"""
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
		### Step 1. Base determinations based on combined detections
		self.logger.debug('Initiated T2MinorPlanetCenter ')


	def run(self, light_curve, run_config):
		""" 
			Parameters
			-----------
				light_curve: `ampel.base.LightCurve` instance. 
				 	See the LightCurve docstring for more info.
			
				run_config: `dict` or None

			Returns
			-------
				dict with entries as in class doc string.
					
					{
						'ndet' : 3,
						...
					}
		"""
		
		self.logger.debug('Checking %s'%(light_curve.id))

		run_config = self.RunConfig() if run_config is None else run_config
		pps = list( light_curve.get_photopoints(filters=run_config.filters) )

		# Check whether we are running for all or only latest 
		if run_config.only_latest:
			pps.sort(key=lambda x: x.get_value('obs_date'))
			pps = [pps[-1]]
			self.logger.debug('Restricting to latest PP at %s'%(pps[0].get_value('obs_date')))
	
		# Loop through remaining pps and check with MPC
		mpc_checks = {}
		for pp in pps:
			print('%s %s %s'%(pp.get_value('ra'),pp.get_value('dec'),pp.get_value('obs_date')))
			# Convert date to UT
			t = Time(pp.get_value('obs_date'), format='jd')
			print(t)
			t.format = 'iso'
			print(t)
#			print(t.to_value('ymdhms', scale='ut1') )

	

		return mpc_checks

