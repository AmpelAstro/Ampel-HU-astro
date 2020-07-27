#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2Observability.py
# License           : BSD-3-Clause
# Author            : matteo.giomi@desy.de
# Date              : 19.09.2018
# Last Modified Date: 05.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

import astropy.units as u
from astropy.time import Time
from typing import Dict, Any
from ampel.util.Observatory import Observatory
from ampel.view.LightCurve import LightCurve
from ampel.type import T2UnitResult
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit


class T2Observability(AbsLightCurveT2Unit):
	"""
	cross match the position of a transient to those of sources in a set
	of catalogs and attach the required information to the transient.
	"""

	# empty dict of named EarthLocation objetcs.
	observatories: Dict = {}

	# default parameters for LightCurve.get_pos method
	lc_get_pos_defaults: Dict[str, Any] = {'ret': "brightest", 'filters': None}

	# default parameters for visibility
	visibility_defaults: Dict[str, Any] = {'airmass_th': 2, 'sun_alt_th': -12, 'min_moon_dist': 30}


	def init_observatory(self, name, position) -> Observatory:
		"""
		Return the earth location object corresponding to the desired observatory.
		Repeated requests to the same instance will not cause new duplicates.
		"""

		# check if the observatory has already been init
		obs = self.observatories.get(name)
		if obs is None:
			self.logger.debug(
				f"Observatory {name} not previously instantiated. Doing it now."
			)

			# init your obs depending on the position
			obs = Observatory(name, *position, logger=self.logger)
			self.observatories[name] = obs
			return obs

		self.logger.debug(f"Observatory {name} already exists.")
		return obs


	def run(self, light_curve: LightCurve) -> T2UnitResult:
		"""
		:param run_config: configuration parameter for this job.
		Eg:
		run_config = {
			'get_lc_pos_kwargs': None, # optional see ampel.view.LightCurve doc
			'observatories': {
				'SEDm': {
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
					{.. },
				...
				}
			}

		:returns: dict with the keys to append to each transient.
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

		self.logger.debug(
			f"getting transient position from lightcurve using args: {lc_get_pos_kwargs}"
		)

		transient_ra, transient_dec = light_curve.get_pos(**lc_get_pos_kwargs)

		self.logger.debug(
			f"Transient position (ra, dec): {transient_ra:.4f}, {transient_dec:.4f} deg"
		)

		# initialize the catalog quer(ies). Use instance variable to aviod duplicates
		out_dict: Dict[str, Any] = {}
		observatories = run_config.get('observatories')

		if observatories is None:
			raise KeyError("run_config missing list of observatories.")

		for name, observatory in observatories.items():

			my_obs = self.init_observatory(name, observatory['pos'])

			for k in range(3):
				trange = [
					(Time.now() + k * u.day).iso[:10],
					(Time.now() + (k + 1) * u.day).iso[:10]
				]

			ret = my_obs.compute_observability(
				transient_ra, transient_dec,
				trange, **observatory['constraints']
			)

			out_dict[name][f'night{k+1}'] = {
				'start': ret[0].iso,
				'end': ret[-1].iso
			}

		# return the info as dictionary
		return out_dict
