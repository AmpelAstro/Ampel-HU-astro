#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2CatalogMatch.py
# License           : BSD-3-Clause
# Author            : matteo.giomi@desy.de
# Date              : 11.09.2018
# Last Modified Date: 11.09.2018
# Last Modified By  : matteo.giomi@desy.de


import logging
logging.basicConfig()

from astropy.table import Table
from scipy.interpolate import interp1d


from ampel.pipeline.common.ZTFUtils import ZTFUtils
from ampel.base.abstract.AbsT2Unit import AbsT2Unit
from ampel.core.flags.T2RunStates import T2RunStates


class T2LightcurveQuality(AbsT2Unit):
	"""
		determine the 'quality' of the light curve by computing ratios between
		the number of detection and that of upper limits.
		
		The LC 'quality' is measured by two numbers:
		
			* 'detection strenght' = n_det / n_obs
			
			* 'detection purity' = n_det / n_det + n_strong_ulims
		
		where:
			n_det: 
				total number of detections
			
			n_obs:
				number of observations (detections + upper lims) computed
				from the time of the first detection.
			
			n_strong_ulims: 
				number of upper limits which are below (higher magnitude) than
				what expected from a simple interpolation between the detections.
				
				That is, be interp_lc the function that iterpolates the detections
				(returning magnitude and accepting a time), an upper limit at 
				time jd_ul of magnitude mag_ul is considered 'strong' if:
					
					interp_lc(jd_ul) < mag_ul
	"""
	
	version = 0.1

	def __init__(self, logger, base_config):
		"""
		"""
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
		# set some default values
		self.default_filter_ids = [1, 2, 3]
		self.filter_names = {1: 'g', 2: 'r', 3: 'i'}
		self.default_lc_filter = [
									{
									'attribute': 'isdiffpos', 
									'operator': '!=', 
									'value': 'f'
									},
									{
									'attribute': 'isdiffpos', 
									'operator': '!=', 
									'value': '0'
									}
								]

	def _count_strong_ulims(self, det_tab, ulim_tab):
		"""
			compute the number of strong upper limts in the light curve. This is
			defined as the number of upper limits which are below (higher magnitude) than
			what expected from a simple interpolation between the detections.
		"""
		
		# interpolate detections
		interp_lc = interp1d(det_tab['jd'], det_tab['mag'], kind='zero', fill_value="extrapolate")
		
		# loop on uls and count the strong ones
		n_strong = 0
		for ul in ulim_tab:
			expected_mag = interp_lc(ul['jd'])
#			self.logger.debug("upper limit at jd %f is at %f, source should be at %f"%
#				(ul['jd'], ul['mag'], expected_mag))
			if ul['mag'] > expected_mag:
				n_strong+=1
		return n_strong
	
	def test_plot(dets, ulims)
		"""
			never used, but useful for debugging
		"""
		
		import matplotlib.pyplot as plt
		import numpy as np
		
		interp_lc = interp1d(dets['jd'], dets['mag'], kind='zero', fill_value="extrapolate")
		
		min_t = min([dets['jd'].min(), ulims['jd'].min()])
		max_t = max([dets['jd'].max(), ulims['jd'].max()])
		jd_int = np.arange(min_t, max_t, 0.1)
		plt.plot(dets['jd'], dets['mag'], 'o', label="data")
		plt.plot(ulims['jd'], ulims['mag'], 'o', label="ulims")
		plt.plot(jd_int, interp_lc(jd_int), label="interp")
		plt.gca().invert_yaxis()
		plt.legend()
		plt.xlabel("JD")
		plt.ylabel("Mag")
		
		# add text
		textstr = '$n_{det}$: %d, $n_{lim}$: %d, $n_{lim}^{strong}$: %d, purity: %.3f, strength: %.3f'%(
								len(dets),len(ulims),n_strong_ulims, purity, strength)
		plt.title("filter: %d "%fid+textstr)
		plt.show()

	def run(self, light_curve, run_config):
		""" 
			Parameters
			-----------
				light_curve: "ampel.base.LightCurve" instance. 
				 See the LightCurve docstring for more info.
			
				run_parameters: `dict`
						configuration parameter for this job. 
						
			
			Returns
			-------
				dict with the keys to append to each transient.
		"""
		
		# parse some run config params and use default if not given
		lc_filter = run_config.get('lc_filter')
		if lc_filter is None:
			lc_filter = self.default_lc_filter
		filter_ids = run_config.get('filter_ids')
		if filter_ids is None:
			filter_ids = self.default_filter_ids

		# run on the single bands individually
		out_dict = {}
		for fid in filter_ids:
			
			self.logger.debug("computing light curve quality for filter id %d (%s-band)"%
				(fid, self.filter_names[fid]))
			
			# concatenate the filters
			filters = [{'attribute': 'fid', 'operator': '==', 'value': fid}]
			if type(lc_filter) in [list, tuple]:
				filters += list(lc_filter)
			elif type(lc_filter) is dict:
				filters += [lc_filter]
			else:
				raise ValueError("parameter 'lc_filter' must be either list or tuple. got %s instead"%(type(lc_filter)))
			self.logger.debug("applying filter: %s"%repr(filters))
			
			# get upper limits and detections time series
			pps = light_curve.get_tuples('obs_date', 'magpsf', filters=filters)
			uls = light_curve.get_tuples('obs_date', 'diffmaglim', filters=filters, upper_limits=True)
			
			# if you have no detections, you're done
			if len(pps) == 0:
				self.logger.debug("No detections in light curve for this band.")
				out_dict[self.filter_names[fid]] = {'strength':0, 'purity':0, 'ndet': 0}
				continue
			
			# also easy
			if len(uls) == 0:
				self.logger.debug("No upper limits in light curve for this band.")
				out_dict[self.filter_names[fid]] = {'strength':1, 'purity':1, 'ndet': len(pps)}
				continue
			
			# cast to tables for convenience
			dets = Table(rows = pps, names=('jd', 'mag'))
			ulims = Table(rows = uls, names=('jd', 'mag'))
			self.logger.debug("got %d detections and %d ulims for lightcurve."%
				(len(dets), len(ulims)))
			
			# compute time of first and last detection and mask out upper limts happening before
			det_start, det_end = dets['jd'].min(), dets['jd'].max()
			ulims = ulims[ulims['jd']>det_start]
			self.logger.debug("retained %d upper limit from start of detection (at %f jd)"%
				(len(ulims), det_start))
			
			# if you don't have any upper limit after first detection, easy
			if len(ulims) == 0:
				out_dict[self.filter_names[fid]] = {'strength':1, 'purity':1, 'ndet': len(dets)}
				continue
			
			# compute number of detections, total observations, and upper limts.
			# for the strong upper limits, exclude those which happends after the last detection
			n_strong_ulims = self._count_strong_ulims(dets, ulims[ulims['jd']<det_end])
			strength = len(dets) / (len(dets) + len(ulims))
			purity = len(dets) / (len(dets) + n_strong_ulims)
			out_dict[self.filter_names[fid]] = {'strength':strength, 'purity':purity, 'ndet': len(dets)}
			
		# return the info as dictionary
		return out_dict
