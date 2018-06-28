#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t0/DecentFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 06.06.2018
# Last Modified Date: 27.06.2018
# Last Modified By  : m. giomi <matteo.giomi@desy.de>


import numpy as np
import logging
from urllib.parse import urlparse
from astropy.coordinates import SkyCoord
from astropy.table import Table
from catsHTM import cone_search
from ampel.abstract.AbsAlertFilter import AbsAlertFilter


class DecentFilter(AbsAlertFilter):
	"""
		General-purpose filter with ~ 0.6% acceptance. It selects alerts based on:
			* numper of previous detections
			* positive subtraction flag
			* loose cuts on image quality (fwhm, elongation, number of bad pixels, and the
			difference between PSF and aperture magnitude)
			* distance to known SS objects
			* real-bogus
			* detection of proper-motion and paralax for coincidence sources in GAIA DR2
		
		The filter has a very weak dependence on the real-bogus score and it is independent
		on the provided PS1 star-galaxy classification.
	"""

	# Static version info
	version = 1.0
	resources = ('catsHTM',)

	def __init__(self, on_match_t2_units, base_config=None, run_config=None, logger=None):
		"""
		"""
		if run_config is None or len(run_config) == 0:
			raise ValueError("Please check you run configuration")

		self.on_match_t2_units = on_match_t2_units
		self.logger = logger if logger is not None else logging.getLogger()
		
		config_params = (
			'MIN_NDET',					# number of previous detections
			'MIN_RB',					# real bogus score
			'MAX_FWHM',					# sexctrator FWHM (assume Gaussian) [pix]
			'MAX_ELONG',				# Axis ratio of image: aimage / bimage 
			'MAX_MAGDIFF',				# Difference: magap - magpsf [mag]
			'MAX_NBAD',					# number of bad pixels in a 5 x 5 pixel stamp
			'MIN_DIST_TO_SSO',			# distance to nearest solar system object [arcsec]
			'GAIA_RS',					# search radius for GAIA DR2 matching [arcsec]
			'GAIA_PM_SIGNIF',			# significance of proper motion detection of GAIA counterpart [sigma]
			'GAIA_PLX_SIGNIF'			# significance of parallax detection of GAIA counterpart [sigma]
			)
		for el in config_params:
			if el not in run_config:
				raise ValueError("Parameter %s missing, please check your channel config" % el)
			if run_config[el] is None:
				raise ValueError("Parameter %s is None, please check your channel config" % el)
			self.logger.info("Using %s=%s" % (el, run_config[el]))
		
		
		# ----- set filter proerties ----- #
		
		# history
		self.min_ndet 					= run_config['MIN_NDET'] 
		
		# Image quality
		self.max_fwhm					= run_config['MAX_FWHM']
		self.max_elong					= run_config['MAX_ELONG']
		self.max_magdiff				= run_config['MAX_MAGDIFF']
		self.max_nbad					= run_config['MAX_NBAD']
		self.min_rb						= run_config['MIN_RB']
		
		# astro
		self.min_ssdistnr	 			= run_config['MIN_DIST_TO_SSO']
		self.gaia_rs					= run_config['GAIA_RS']
		self.gaia_pm_signif				= run_config['GAIA_PM_SIGNIF']
		self.gaia_plx_signif			= run_config['GAIA_PLX_SIGNIF']

		# technical
		self.catshtm_path 			= urlparse(base_config['catsHTM']).path
		self.logger.info("using catsHTM files in %s"%self.catshtm_path)
		self.keys_to_check = (
			'fwhm', 'elong', 'magdiff', 'nbad', 
			'isdiffpos', 'ra', 'dec', 'rb', 'ssdistnr')


	def _alert_has_keys(self, photop):
		"""
			check that given photopoint contains all the keys needed to filter
		"""
		
		for el in self.keys_to_check:
			if el not in photop:
				self.logger.debug("rejected: '%s' missing" % el)
				return False
			if photop[el] is None:
				self.logger.debug("rejected: '%s' is None" % el)
				return False
		return True


	def is_star_in_gaia(self, transient):
		"""
			match tranient position with GAIA DR2 and uses parallax 
			and proper motion to evaluate star-likeliness
			
			returns: True (is a star) or False otehrwise.
		"""
		transient_coords = SkyCoord(transient['ra'], transient['dec'], unit='deg')
		srcs, colnames, colunits = cone_search(
											'GAIADR2',
											transient_coords.ra.rad, transient_coords.dec.rad,
											self.gaia_rs,
											catalogs_dir=self.catshtm_path)
		my_keys = ['RA', 'Dec', 'Mag_G', 'PMRA', 'ErrPMRA', 'PMDec', 'ErrPMDec', 'Plx', 'ErrPlx']
		if len(srcs) > 0:
			gaia_tab					= Table(srcs, names=colnames)
			gaia_tab					= gaia_tab[my_keys]
			gaia_coords					= SkyCoord(gaia_tab['RA'], gaia_tab['Dec'], unit='rad')
			
			# compute distance
			gaia_tab['DISTANCE']		= transient_coords.separation(gaia_coords).arcsec
			gaia_tab['DISTANCE_NORM']	= (
				1.8 + 0.6 * np.exp( (20 - gaia_tab['Mag_G']) / 2.05) > gaia_tab['DISTANCE'])
			gaia_tab['FLAG_PROX']		= [True if x['DISTANCE_NORM'] == True and 11 <= x['Mag_G'] <= 19 else False for x in gaia_tab]

			# check for proper motion and parallax conditioned to distance
			gaia_tab['FLAG_PMRA']		= abs(gaia_tab['PMRA']  / gaia_tab['ErrPMRA']) > self.gaia_pm_signif
			gaia_tab['FLAG_PMDec']		= abs(gaia_tab['PMDec'] / gaia_tab['ErrPMDec']) > self.gaia_pm_signif
			gaia_tab['FLAG_Plx']		= abs(gaia_tab['Plx']   / gaia_tab['ErrPlx']) > self.gaia_plx_signif
			if (any(gaia_tab['FLAG_PMRA'] == True) or 
				any(gaia_tab['FLAG_PMDec'] == True) or
				any(gaia_tab['FLAG_Plx'] == True)) and any(gaia_tab['FLAG_PROX'] == True):
				return True
		return False


	def apply(self, alert):
		"""
		Mandatory implementation.
		To exclude the alert, return *None*
		To accept it, either return
			* self.on_match_t2_units
			* or a custom combination of T2 unit names
		"""
		
		# --------------------------------------------------------------------- #
		#					CUT ON THE HISTORY OF THE ALERT						#
		# --------------------------------------------------------------------- #
		
		npp = len(alert.pps)
		if npp < self.min_ndet:
			self.logger.debug("rejected: %d photopoints in alert (minimum required %d)"% 
				(npp, self.min_ndet))
			return None
		
		# --------------------------------------------------------------------- #
		#							IMAGE QUALITY CUTS							#
		# --------------------------------------------------------------------- #
		
		latest = alert.pps[0]
		if not self._alert_has_keys(latest):
			return None
		
		# ---- image quality cuts ----- #
		if (latest['isdiffpos'] == 'f' or latest['isdiffpos'] == '0'):
			self.logger.debug("rejected: 'isdiffpos' is %s", latest['isdiffpos'])
			return None
		
		if latest['rb'] < self.min_rb:
			self.logger.debug("rejected: RB score %.2f below threshod (%.2f)"%
				(latest['rb'], self.min_rb))
			return None
		
		if latest['fwhm'] > self.max_fwhm:
			self.logger.debug("rejected: fwhm %.2f above threshod (%.2f)"%
				(latest['fwhm'], self.max_fwhm))
			return None
		
		if latest['elong'] > self.max_elong:
			self.logger.debug("rejected: elongation %.2f above threshod (%.2f)"%
				(latest['elong'], self.max_elong))
			return None
		
		if abs(latest['magdiff']) > self.max_magdiff:
			self.logger.debug("rejected: magdiff (AP-PSF) %.2f above threshod (%.2f)"% 
				(latest['magdiff'], self.max_magdiff))
			return None
		
		# --------------------------------------------------------------------- #
		#								ASTRONOMY								#
		# --------------------------------------------------------------------- #
		
		# check for closeby ss objects
		if (0 < latest['ssdistnr'] < self.min_ssdistnr):
			self.logger.debug("rejected: solar-system object close to transient (max allowed: %d)."%
				(self.min_ssdistnr))
			return None
		
		# check with gaia
		if self.is_star_in_gaia(latest):
			self.logger.debug("rejected: within %.2f arcsec from a GAIA start (PM of PLX)" % 
				(self.gaia_rs))
			return None
		
		# congratulation alert! you made it!
		return self.on_match_t2_units
	
