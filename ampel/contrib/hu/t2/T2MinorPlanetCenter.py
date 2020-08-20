#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File			  	: ampel/contrib/hu/t2/T2MinorPlanetCenter.py
# License		   	: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de, simeon.reusch@desy.de
# Date			  	: 10.01.2019
# Last Modified Date: 13.02.2019
# Last Modified By  : simeon.reusch@desy.de


from typing import Optional, Dict, List, Union, Any
from pydantic import BaseModel, BaseConfig
from astropy.time import Time
from ampel.view.LightCurve import LightCurve
from ampel.type import T2UnitResult
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit


class T2MinorPlanetCenter(AbsLightCurveT2Unit):
	"""
	Check if the *latest* detection of a transient corresponds
	matches something known by the MinorPlanetCenter.
	"""


	class RunConfig(BaseModel):
		""" Necessary class to validate configuration.  """
		class Config(BaseConfig):
			""" Raise validation errors if extra fields are present """
			allow_extra = False
			ignore_extra = False
		
		# Will only match the latest photopoint.
		only_latest		: bool	= True		
		# Search radius passed to MPC (arcminutes!!!)
		searchradius	: float = 1			
		# V-band magnitude limit passed to MPC
		maglim		: float	= 22		
		# Potential filter for photopoint selection
		filters: Optional[Union[Dict, List[Dict]]] = None

		# Potential filter for photopoint selection
		filters: Optional[Union[Dict, List[Dict]]] = None


	def post_init(self):
		"""
		"""		
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
				dict with entries as in class doc string. Each observation date checked gets an entry 
				in the dict which contains the number of individual matches 'ndet', the angular distances 
				of the matches in degree and the magnitudes of the matches.
					
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

		self.logger.debug(f'Checking {light_curve.id}')

		run_config = self.RunConfig() if run_config is None else run_config
		pps = list(
			light_curve.get_photopoints(filters=run_config.filters)
		)

		# Check whether we are running for all or only latest
		if run_config.only_latest:
			pps.sort(key=lambda x: x.body.get('obs_date'))
			pps = [pps[-1]]
			self.logger.debug(
				f"Restricting to latest PP at {pps[0].body.get('obs_date')}"
			)

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
