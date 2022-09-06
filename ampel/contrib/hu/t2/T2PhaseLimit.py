#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2PhaseLimit.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 28.04.2021
# Last Modified Date: 21.03.2022
# Last Modified By  : mf@physik.hu-berlin.de

from typing import Dict, List, Optional, Sequence, Any, Union, Iterable
from astropy.coordinates import SkyCoord
import numpy as np
from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.view.T2DocView import T2DocView
import matplotlib.pyplot as plt


from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit

class T2PhaseLimit(AbsStateT2Unit, AbsTabulatedT2Unit):
	"""
	Lightcurve analysis tools assume a certain transient life-time (order weeks for SNe).
	These can be confused by spurious early/late detections and/or background
	variability at the position of the SN.
	This unit tries to identify the most likely time-range for a SN explosion and
	evaluate whether sufficient data for continued exploration exist.
	"""

	# *Conservative* estimates of how fast the transient rises and declines (in days)
	# These times are derived relative to the *median* detection date,
	# and are thus taken as symmetric around this
	half_time : float

	# Min number of detections remaining in the target range for subsequent analysis
	min_det : int = 3

	# Rejection sigma.
	rej_sigma : float = 5.  # A lower number will more aggressively reject data

	# Spurious light in references typically manifests as negative flux.
	# Will require the fraction of obs with negative to be lower than this fraction, in each filter.
	neg_frac_lim : float = 0.05

	# Max magnitude to consider (constant low flux can correspond to flux in reference)
	max_flux : float

	# Plot
	plot : bool = False

	def _magtoflux(self, mag: float) -> float:
		return 10 ** (-((mag) - 25) / 2.5)
	
	def __init__(self, *args, **kwargs):
		if not "max_flux" in kwargs:
			if not "max_mag" in kwargs:
				max_mag = 22.5 #previous default value
			else:
				max_mag = kwargs.pop("max_mag")
			kwargs["max_flux"] = self._magtoflux(max_mag)
		super().__init__(*args, **kwargs)

	def process(self,
		compound: T1Document,
		datapoints: Iterable[DataPoint],
	) -> Union[UBson, UnitResult]:
		"""
		Iteratively reject datapoints with outlying phase estimate.
		Once completed, check whether remaining data is compatible with expectations.
		Note: This unit does not differentiate between filters.
		Note: Upper limits are currently not taken into consideration.
		:returns: dict with estimate of likely peak and median time,
			tentative start and end-date and whether the phase looks like a sn.
		E.g.:
		{
			't_start' : 243101,
			't_end' : 243199,
			't_masked_duration' : 40,
			't_peak' : 243120,
			't_median' : 243124,
			'mag_peak' : 19.,
			'max_neg_frac' : 0.02,
			'n_remaining' : 5,
			't2eval' : 'OK' | 'fail:duration' | 'fail:detections' | 'fail:neg_flux' | 'warning:neg_iband'
		}
		"""

		# Starting set of dates (with positive flux )
		# jd_ispos = light_curve.get_tuples( 'jd', 'isdiffpos',
		# 	filters = [{'attribute' : 'magpsf', 'operator' : '<', 'value' : self.max_mag},
		# 		{'attribute' : 'magpsf', 'operator' : '>', 'value' : 0} ]
		# 	 )
		flux_table = self.get_flux_table(datapoints)
		fflux_table = flux_table[(flux_table["flux"]>self.max_flux) & (flux_table["flux"]<1e10)]
		
		# Enough data to even start evaluation?
		if len(fflux_table) < self.min_det :
			return {'t_start' : None, 't_end' : None,
				't_masked_duration' : None, 't_peak':None,
				't_median' : None, 'flux_peak' : None,
				'n_remaining' : None, 't2eval' : 'fail:detections' }

		jd = fflux_table["time"]
		t_median = np.median( jd )
		mask = (jd>0)

		# Iteratively reject datapoints until we reach a minimum (or find no more)
		# Using the median absolute deviation
		while( sum(mask)>= self.min_det):
			sig_est = 1.48 * np.median( np.abs(jd[mask]-t_median) )
			new_mask = ( np.abs(jd-t_median) < self.rej_sigma * sig_est )
			if sum(new_mask)<sum(mask):
				mask = new_mask
				t_median = np.median( jd[mask] )
			else:
				break

		if sum(mask)>0:
			t_masked_duration = np.max( jd[mask] ) - np.min( jd[mask] )
		else:
			t_masked_duration = 0

		# Based on the median, define course phase range
		t_start = t_median - 2 * self.half_time
		t_end = t_median + 2 * self.half_time

		# (Re) retrieve data and magnitudes in this range
		# dps = light_curve.get_tuples('jd','magpsf',filters=
		# 	[{'attribute': 'jd', 'operator': '>', 'value': t_start},
		# 	{'attribute': 'jd', 'operator': '<', 'value': t_end},
		# 	{'attribute': 'magpsf', 'operator': '>', 'value': 0},
		# 	] )
		dps = flux_table[(flux_table["time"]>t_start) & (flux_table["time"]<t_end) & ( flux_table["flux"] < 1e10 )] # do we need to check for magpsf>0?
		if len(dps) > 0:
			dps.sort("flux")
			dps.reverse()
			t_peak = dps["time"][0]
			flux_peak = dps["flux"][0]
			n_remaining = len(dps)
		else:
			n_remaining = 0
			t_peak = None
			flux_peak = None

		# Tighter phase range based on peak time
		if t_peak is not None:
			t_start = min( [t_median, t_peak] ) - self.half_time
			t_end = max( [t_median, t_peak] ) + self.half_time
		else:
			t_start = t_median - self.half_time
			t_end = t_median + self.half_time

		# Investigate negative flux
		neg_frac_bands = []
		for filtid in np.unique(flux_table["band"]):
			in_band = flux_table[flux_table["band"]==filtid]
			len_ispos = len(in_band[in_band["flux"]>0])
			filter_frac = 1 - float( len_ispos ) / len(in_band)
			if filter_frac > self.neg_frac_lim:
				neg_frac_bands.append( filtid )

		# Evaluate
		if any(band.endswith('g') for band in neg_frac_bands):
			t2eval = 'fail:neg_flux'
		elif any(band.endswith('r') for band in neg_frac_bands):
			t2eval = 'fail:neg_flux'
		elif t_masked_duration > 2*self.half_time:
			t2eval = 'fail:duration'
		elif n_remaining < self.min_det:
			t2eval = 'fail:detections'
		elif any(band.endswith('i') for band in neg_frac_bands):
			t2eval = 'warning:neg_iband'
		else:
			t2eval = 'OK'

		# Do plot?
		if self.plot:

			fig = plt.figure(figsize=(6,5) )

			all_jd = flux_table["time"]
			_1, bins, _2 = plt.hist( all_jd, bins=100, label = 'All alerts' )

			plt.hist( jd, bins=bins, label = 'Clipped det.' )
			plt.hist( jd[mask], bins=bins, label = 'Retained det.' )

			plt.axvline( x=t_start, color='k', linestyle='dashed',
				label='Start/end' )
			plt.axvline( x=t_end, color='k', linestyle='dashed')

			if t_peak is not None:
				plt.axvline( x=t_peak, label = 't(peak)' )
			plt.axvline( x=t_median, label= 't(median)' )
			stock_id = '-'.join(str(self.get_stock_id(datapoints)))
			plt.title( f'{stock_id} {t2eval}')

			plt.legend(loc='best')
			plt.savefig( f'./t2phase/{stock_id}.svg')
			plt.close()


		return {'t_start' : t_start, 't_end' : t_end,
			't_masked_duration' : t_masked_duration, 't_peak':t_peak,
			't_median' : t_median, 'flux_peak' : flux_peak,
			'n_remaining' : n_remaining, 't2eval' : t2eval }
