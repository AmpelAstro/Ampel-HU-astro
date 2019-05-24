#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2Observability.py
# License           : BSD-3-Clause
# Author            : matteo.giomi@desy.de
# Date              : 19.09.2018
# Last Modified Date: 21.09.2018
# Last Modified By  : ulrich.feindt@fysik.su.se

import logging
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from ampel.base import Observatory
from ampel.base.abstract.AbsT2Unit import AbsT2Unit
from ampel.core.flags.T2RunStates import T2RunStates


class T2Observability(AbsT2Unit):
	"""
		cross match the position of a transient to those of sources in a set
		of catalogs and attach the required information to the transient.
	"""
	
	version = 0.1

	def __init__(self, logger, base_config):
		"""
		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
		# empty dict of named EarthLocation objetcs.
		self.observatories = {}
		
		# default parameters for LightCurve.get_pos method
		self.lc_get_pos_defaults = {'ret': "brightest", 'filters': None}
		
		# default parameters for visibility
		self.visibility_defaults = {'airmass_th': 2, 'sun_alt_th': -12, 'min_moon_dist': 30}
		

	def init_observatory(self, name, position):
		"""
			Return the earth location object corresponding to the desired observatory.
			Repeated requests to the same instance will not cause new duplicates.
			
			Returns:
			--------
				
				Observatory instance.
		"""
		
		# check if the observatory has already been init
		obs = self.observatories.get(name)
		if obs is None:
			self.logger.debug("Observatory %s not previously instantiated. Doing it now."%name)
			
			#init your obs depending on the position
			obs = Observatory(name, *position, logger=self.logger)
			self.observatories[name] = obs
			return obs
		else:
			self.logger.debug("Observatory %s already exists."%name)
			return obs

	def run(self, light_curve, run_config):
		""" 
			Parameters
			-----------
				light_curve: "ampel.base.LightCurve" instance. 
				 See the LightCurve docstring for more info.
			
				run_parameters: `dict`
						configuration parameter for this job. 
						
						Eg:
						
						run_config = 
							{
							'get_lc_pos_kwargs': None, # optional see ampel.base.LightCurve doc
							'observatories':
								{
								'SEDm': 
									{	
										pos: {
											'lat': 33.3483717		#from google maps, and for p48
											'long': -116.85972959
											'alt': 1680
										},
										constraints: {
											'airmass_th': 2,
											'sun_alt_th': -12,
											'min_moon_dist': 30
										}
									},
								'SNIFS': 
									{
										..
									},
								...
								}
							}
			
			Returns
			-------
				dict with the keys to append to each transient. 
				
				{
					obs1: {
							night1: {
								start,
								stop,
								},
							night2: {
								...
							},
							night3: {
							}
					},
					obs2:
						{
						...
					},
					..
				}
		"""
		
		# get ra and dec from lightcurve object
		lc_get_pos_kwargs = run_config.get('lc_get_pos_kwargs')
		if lc_get_pos_kwargs is None:
			lc_get_pos_kwargs = self.lc_get_pos_defaults
		self.logger.debug("getting transient position from lightcurve using args: %s", lc_get_pos_kwargs)
		transient_ra, transient_dec = light_curve.get_pos(**lc_get_pos_kwargs)
		self.logger.debug("Transient position (ra, dec): %.4f, %.4f deg"%(transient_ra, transient_dec))
		
		# initialize the catalog quer(ies). Use instance variable to aviod duplicates
		out_dict = {}
		observatories = run_config.get('observatories')
		if observatories is None:
			raise KeyError("run_config missing list of observatories.")
		
		for name, observatory in observatories.items():
			
			my_obs = self.init_observatory(name, observatory['pos'])

			for k in range(3):
					trange = [(Time.now() + k*u.day).iso[:10], (Time.now() + (k+1)*u.day).iso[:10]]
			ret = my_obs.compute_observability(transient_ra, transient_dec, trange, **observatory['constraints'])

			out_dict[name]['night%i'%(k+1)] = {'start': ret[0].iso, 'end': ret[-1].iso}
			
		# return the info as dictionary
		return out_dict
