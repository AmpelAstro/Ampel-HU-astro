#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File			  : ampel/contrib/hu/t2/T2LCQuality.py
# License		   : BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date			  : 10.01.2019
# Last Modified Date: 13.02.2019
# Last Modified By  : jnordin@physik.hu-berlin.de


import logging
logging.basicConfig()

import requests
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import time
from bs4 import BeautifulSoup

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
		searchradius	: float = 1			# Searchradius passed to MPC (arcminutes!!!)
		maglim		: float	= 22		# V-band magnitude limit passed to MPC
		filters			: dict	= None		# Potential filter for photopoint selection



	def __init__(self, logger, base_config):
		"""
		"""
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
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
				dict with entries as in class doc string. Each observation date checked gets an entry in the dict which contains the number of individual matches 'ndet', the angular distances of the matches in degree and the magnitudes of the matches.
					
					{
						obs_date : 
						{
							'ndet': number of detections (float), 
							'ang_distances_deg': angular distances in degree (list of floats or None), 
							'mags': mangitudes (list of floats or None)
						},
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
		angular_separation_deg = []
		mag_vband = []

		NEO_URL = "https://cgi.minorplanetcenter.net/cgi-bin/mpcheck.cgi"

		
		for pp in pps:
			ra = pp.get_value('ra')
			dec = pp.get_value('dec')

			self.logger.debug('Checking MPC',extra={'ra':ra, 'dec':dec,  'obs_date':{pp.get_value('obs_date')} } )

			# Convert date for HTTP request
			t = time.Time(pp.get_value("obs_date"), format="jd", scale="utc")
			year = t.strftime("%Y")
			month = t.strftime("%m")
			day = t.strftime("%d")
			daydecimal = t.mjd - np.fix(t.mjd)
			daydecimal = str(daydecimal).split(".")[1]
			day = day + "." + daydecimal
			day = np.around(np.float(day), decimals=2)

			# Convert coordinates for HTTP request
			radec_skycoord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
			ra = radec_skycoord.ra.to_string(u.hour, sep=" ",pad=True)
			dec = radec_skycoord.dec.to_string(u.deg, sep=" ",pad=True)


			request_data = {"year": f"{year}", "month": f"{month}", "day": f"{day}", "which": "pos",
				"ra": f"{ra}", "decl": f"{dec}", "TextArea": "", "radius": f"{run_config.searchradius}", 
				"limit": f"{run_config.maglim}", "oc": "500", "sort": "d", "mot": "h", 
				"tmot": "s", "pdes": "u", "needed": "f", "ps": "n", "type": "p"}

			# Post the request
			response = requests.post(url=NEO_URL, data=request_data, timeout=30)

			# Parse the result
			soup = BeautifulSoup(response.text, 'html5lib')

			try:
				pre = soup.find_all('pre')[-1]
				results = pre.text.lstrip(" ").split("\n")[3:]
				separations = []
				mags = []
				for result in results:
					if len(result) > 10:
						radec = result[25:46]
						mag = float(result[47:51])
						skycoord = SkyCoord(radec, unit=(u.hourangle, u.deg))
						sep = skycoord.separation(radec_skycoord)
						separations.append(sep.deg)
						mags.append(mag)
			except IndexError:
				separations = []
				mags = []
				pass

			if len(separations) == 0:
				result = {t.jd: {"ndet": 0, "ang_distances_deg": None, "mags": None}}
			else:
				result = {t.jd: {"ndet": len(separations), "ang_distances_deg": separations, "mags": mags}}

			mpc_checks.update(result)

		return mpc_checks
