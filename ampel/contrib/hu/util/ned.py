#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/util/ned.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                14.09.2021
# Last Modified Date:  14.09.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from typing import Tuple
from ampel.protocol.LoggerProtocol import LoggerProtocol

def check_ned_res(
	cat_res: dict,
	logger: LoggerProtocol,
	spectroscopic: bool = False,
	z_range: None | tuple[float, float] = None
) -> bool:

	if not cat_res.get('z'):
		logger.info("No redshift found in NED result")
		return True

	if spectroscopic and cat_res.get('n_spectra', 0) == 0 and cat_res["zflag"] != "SPEC":
		logger.info("Not a spectroscopic redshift")
		return True

	if z_range and (cat_res['z'] < z_range[0] or cat_res['z'] > z_range[1]):
		logger.info("Redshift exceeds allowed values")
		return True

	return False
