#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File			  	: ampel/contrib/hu/t2/T2MinorPlanetCenter.py
# License		   	: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de, simeon.reusch@desy.de
# Date			  	: 10.01.2019
# Last Modified Date: 13.02.2019
# Last Modified By  : simeon.reusch@desy.de


from typing import Optional, Dict, List, Union, ClassVar
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import requests
from bs4 import BeautifulSoup # type: ignore[import]

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.T2IngestOptions import T2IngestOptions


class T2MinorPlanetCenter(AbsPointT2Unit):
	"""
	Check if the *latest* detection of a transient corresponds
	matches something known by the MinorPlanetCenter.
	"""

	# run only on last datapoint by default
	# NB: this assumes that docs are created by DualPointT2Ingester
	eligible: ClassVar[Optional[T2IngestOptions]] = {"filter": "PPSFilter", "sort": "jd", "select": "first"}

	# Search radius passed to MPC (arcminutes!!!)
	searchradius : float = 1
	# V-band magnitude limit passed to MPC
	maglim : float = 22
	# Potential filter for photopoint selection
	filters: Optional[Union[Dict, List[Dict]]] = None

	def process(self, pp: DataPoint) -> Union[UBson, UnitResult]:
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
		NEO_URL = "https://cgi.minorplanetcenter.net/cgi-bin/mpcheck.cgi"

		try:
			ra = pp['body']['ra']
			dec = pp['body']['dec']
			jd = pp['body']['jd']
		except KeyError:
			return UnitResult(code=DocumentCode.T2_MISSING_INFO)

		self.logger.debug('Checking MPC',extra={'ra':ra, 'dec':dec,'jd':jd } )

		# Convert date for HTTP request
		t = Time(jd, format="jd", scale="utc")
		year = t.strftime("%Y")
		month = t.strftime("%m")
		day = t.strftime("%d")
		daydecimal = t.mjd - np.fix(t.mjd)
		daydecimal = str(daydecimal).split(".")[1]
		day = day + "." + daydecimal
		day = np.around(np.float64(day), decimals=2)

		# Convert coordinates for HTTP request
		radec_skycoord = SkyCoord(ra, dec, unit=(u.deg, u.deg))
		ra = radec_skycoord.ra.to_string(u.hour, sep=" ",pad=True)
		dec = radec_skycoord.dec.to_string(u.deg, sep=" ",pad=True)


		request_data = {"year": f"{year}", "month": f"{month}", "day": f"{day}", "which": "pos",
			"ra": f"{ra}", "decl": f"{dec}", "TextArea": "", "radius": f"{self.searchradius}", 
			"limit": f"{self.maglim}", "oc": "500", "sort": "d", "mot": "h", 
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

		return result
