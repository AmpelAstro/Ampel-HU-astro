#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/t3/HealpixTokenGenerator.py
# License:             BSD-3-Clause
# Author:              jnordin
# Date:                24.03.2023
# Last Modified Date:  24.03.2023
# Last Modified By:    jnordin <jnordin@physik.hu-berlin.de>

import random
import time
from typing import Any
from astropy.time import Time  # type: ignore
from datetime import datetime
from requests_toolbelt.sessions import BaseUrlSession

from ampel.types import UBson
from ampel.struct.T3Store import T3Store
from ampel.struct.Resource import Resource
from ampel.struct.UnitResult import UnitResult
from ampel.secret.NamedSecret import NamedSecret
from ampel.abstract.AbsT3PlainUnit import AbsT3PlainUnit

from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix, deres



class HealpixTokenGenerator(AbsT3PlainUnit):
	'''
		Based on a URL to a Healpix map:
		- find pixels given requested prob contour.
		- request archive token for this stream.
        '''

	# Process pixels with p-values lower than this limit
	pvalue_limit: float = 0.9

	# Name (signifier)
	map_name: str

	# URL for healpix retrieval
	map_url: str
	map_dir: str    # Local dir where map is saved. File with this name del


	archive_token: NamedSecret[str] = NamedSecret(label="ztf/archive/token")
	#: Base URL of archive service
	archive: str = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/"


	date_str: None | str = None
	date_format: str = "%Y-%m-%d"
	delta_time: None | float = None  


	#: seconds to wait for query to complete
	timeout: float = 60

	debug: bool = False



	def process(self, t3s: T3Store) -> UBson | UnitResult:
	
		# Retrieve and process map
		ah = AmpelHealpix(map_name=self.map_name, map_url=self.map_url, save_dir = self.map_dir)
		map_hash = ah.process_map()

		# Get list of pixels within requested significance contour
		pixels = ah.get_pixelmask(self.pvalue_limit)
		self.logger.info('', extra={'map':self.map_name, 'hash':map_hash, 'size':len(pixels), 'nside':ah.nside})
				
		# JD time range		
		if self.delta_time:
			if self.date_str:
				end_jd = Time(
					str(datetime.strptime(self.date_str, self.date_format)),
					format="iso", scale="utc"
				).jd
			else:
				end_jd = Time.now().jd
			start_jd = end_jd - self.delta_time
		else:
			start_jd = 2459898.9
			end_jd = 2459899.



		session = BaseUrlSession(self.archive if self.archive.endswith("/") else self.archive + "/")
		session.headers["authorization"] = f"bearer {self.archive_token.get()}"

		# Combine pixels when possible
		deresdict = deres(ah.nside, pixels)
		healpix_regions = [ {"nside": nside, "pixels": members} for nside, members in deresdict.items() ]
		count = sum( [len(region['pixels']) for region in healpix_regions] )

		response = session.post(
			"streams/from_query",
			json = {
				"jd": {"$gt": start_jd, "$lt": end_jd},
				"regions": healpix_regions,
			}
		)
		response.raise_for_status()

		rd = response.json()
		try:
			token = rd.pop("resume_token")
		except KeyError as exc:
			raise ValueError(f"Unexpected response: {rd}") from exc

		# wait for query to finish - is this needed, or handled by alert consumer?
		t0 = time.time()
		delay = 1
		while time.time() - t0 < self.timeout:
			response = session.get(f"stream/{token}")
			if response.status_code != 423:
				break
			time.sleep(random.uniform(0, delay))
			delay *= 2
		else:
			raise RuntimeError(f"{session.base_url}stream/{token} still locked after {time.time() - t0:.0f} s")
		self.logger.info("Stream created", extra=response.json())

	
		# Package resource needed
		resource = {
		    'map_name': self.map_name,
		    'hash': map_hash,
		    'token': token	
		    }
	
		r = Resource(name=self.map_name, value=resource)
		t3s.add_resource(r)
		# Also add just the token for direct use by ZTFArchiveAlertLoader
		r = Resource(name=self.map_name+'_token', value=token)
		t3s.add_resource(r)
		

		if self.debug:
			return r.dict()

		return None


		


