#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2MinorPlanetCenter.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 10.01.2019
# Last Modified Date: 05.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import Optional, Dict, List, Union, Any
from pydantic import BaseModel, BaseConfig
from astropy.time import Time
from ampel.view.LightCurve import LightCurve
from ampel.abstract.AbsT2Unit import AbsT2Unit


class T2MinorPlanetCenter(AbsT2Unit):
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

		# Will only match the latest photopoint
		only_latest: bool = True

		# Potential filter for photopoint selection
		filters: Optional[Union[Dict, List[Dict]]] = None


	def post_init(self):
		"""
		"""
		self.logger.debug('Initiated T2MinorPlanetCenter ')


	def run(self, light_curve: LightCurve, run_config: Dict):
		"""
		:returns: dict with entries as in class doc string.
		{'ndet' : 3, ...}
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
		mpc_checks: Dict[str, Any] = {}
		for pp in pps:
			print(
				'%s %s %s' % (
					pp.body.get('ra'),
					pp.body.get('dec'),
					pp.body.get('obs_date')
				)
			)
			# Convert date to UT
			t = Time(pp.body.get('obs_date'), format='jd')
			print(t)
			t.format = 'iso'
			print(t)
#			print(t.to_value('ymdhms', scale='ut1') )

		return mpc_checks
