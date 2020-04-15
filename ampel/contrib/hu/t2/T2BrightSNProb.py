#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2GuessSN.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 06.04.2020
# Last Modified Date: 14.04.2020
# Last Modified By  : jnordin@physik.hu-berlin.de


import sys, os
import logging
logging.basicConfig()

from astropy.table import Table, unique
import numpy as np
#import xgboost as xgb
import pkgutil

from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils
from ampel.base.abstract.AbsT2Unit import AbsT2Unit
from ampel.core.flags.T2RunStates import T2RunStates
import ampel.contrib.hu.t2.xgb_trees as xgb_trees


class T2BrightSNProb(AbsT2Unit):
	"""
		Derive a number of simple metrics describing the rise, peak and decline of a lc.
		Run a XGB tree trained to check whether this transient are likely to be an RCF SN.
		
		Derived values:
		* t_predetect : time between final good upper limit and first detection
		* t_lc : duration (time between first and most recent detection)
		* jd_max : jd of peak light. None unless bool_peaked
		* jd_det : jd of first detection
		* jd_last : jd of last detection
		* ndet : number of detections
                * bool_peaked : is the lc estimated to be declining?
		* bool_pure : has there been no significant non-detections after first detection?
		* bool_rising : was the peak light within X days of the most recent detection?
		* bool_norise : was the first detection NOT significantly fainter than
					mag_peak IF bool_peaked, ELSE mag_lst
		* bool_hasgaps : The lc has a gap between detections of at least 30 days, 
				indicating either a recurrent event or a chance coincidental detection. 
		* mag_peak : magnitude at peak light (any band). Only calculated if bool_peaked
		* mag_det : detection magnitude (any band)
		* mag_last : magnitude of last detection (any band)
		* col_peak : color at peak. None unless bool_peaked AND g+r obs made within X days of jd_max
		* col_det : color at detection. None unless g&r obs made within X days of jd_det
		* col_last : color at last detection. None unless g&r obs made within X days of jd_last
		* slope_rise_{g,r} : magnitude slope between jd_det and jd_max. None if bool_norise
		* slope_decline_{g,r} : magnitude slope between jd_max and jd_lst. None unless bool_peaked
		* rb_med : median Real Bogus
		* drb_med : median dReal Bogus (if available)
		Additionally, will try to save the following host related properties:
		'distnr','magnr','classtar','sgscore1','distpsnr1','sgscore2','distpsnr2','neargaia','maggaia'
		

		Note 1:
		What is considered "simultaneus" for calculating e.g. colors is determined by the X parameter.
		The default value of X is 2 days.

		Note 2:
		Upper limits are only considered for the calculation of t_predect and bool_pure. 
		To guarantee that upper limits would not have affected entries only use evaluations 
		with bool_pure=True.
		Upper limits below minUL (default 20) are considered weak and are not considered.  

		Note 3:
		Only positive detections (isdiffpos) are included. 
		Further control of what observations are included can be obtained through the lc filter
		run parameter.

		Selection code taken from T2RiseDeclineStat

		A subset of these are then used as input to an XGB classifier
		trained to single out objects which are likely to evolve
		into a nearby SN.

		
	"""
	
	version = 0.1

	def __init__(self, logger, base_config):
		"""
		"""
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
		# Set some default values (Note that these are ZTF specific)
		self.default_filter_names = {1: 'g', 2: 'r', 3: 'i'}
		self.default_filter_ids = [1, 2]
		self.default_max_tsep = 2
		self.default_min_UL = 20
		self.default_max_tgap = 30
		self.default_lc_filter 	= [
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

		# For some reason we have duplicate photopoints. Why!!!
		# Through setting this we manually just keep the first of each occurance 
		# Check ZTF18aacbccj,ZTF17aaaekyn for example. Has something to do whether a detection
		# is from a prv history or is new. Will for now use the one with highest rb
		# ZTF18aaaorhy is another funny case with a repeated mag.
		self.del_duplicate_rows = True
	
		# Files and selection parameters for tree. Dict keys correspond to nbr of detections.
		# The other parameters were requirements for inclusion in training (and thus usage)
		self.model_params = ['mag_det','mag_last','t_lc','rb_med','col_det','t_predetect','distnr_med', 'magnr_med', 'classtar_med','sgscore1_med', 'distpsnr1_med', 'neargaia_med', 'maggaia_med']
		self.xgb_tree_param = {
					2 : {
						'max_duration' : 3.5,
						'max_predetect' : 3.5,
						'min_detmag' : 16
					},
					3 : {
						'max_duration' : 6.5,
						'max_predetect' : 3.5,
						'min_detmag' : 16
					},
					4 : {
						'max_duration' : 6.5,
						'max_predetect' : 3.5,
						'min_detmag' : 16
					},
					5 : {
						'max_duration' : 10,
						'max_predetect' : 3.5,
						'min_detmag' : 16
					},
					6 : {
						'max_duration' : 10,
						'max_predetect' : 3.5,
						'min_detmag' : 16
					},
					100 : {
						'max_duration' : 90,
						'max_predetect' : 10,
						'min_detmag' : 16
					}
				}
		# Load the (large) set of trees
		self.xgb_tree = xgb_trees.xgboost_tree()


	def run(self, light_curve, run_config):
		""" 
			Parameters
			-----------
				light_curve: `ampel.base.LightCurve` instance. 
				 	See the LightCurve docstring for more info.
			
				run_config: `dict` or None
					configuration parameter for this job. If none is given, the
					default behaviour would be to compute the metrics for the light
					curve based on g and r with filter ID as for ZTF, to exclude points 						with negative detectiions (having isdiffpos in ['f', 0]), 
					to consider upper limits larger than 20 and using detections within 2
					days when determining roughly same time claculations.
					
					These defaults can be changed by the following keys of the 
					run_config dictionary:
					
						lc_filter: `dict` or `list`
							to be passed to ampel.base.LightCurve.get_tuples.
							if list, the items must be dicts and they'll be combined 
							with a logical and. Pass an empy list to disable the filter
							completely (filtering on the ztf bands will still be applied).
							A sample usage would be RB / dRB cuts.
							This cut is not applied to upper limits.
						
						filter_ids: `list` or `tuple`
						filter_names: `dict`
							Pair of IDs of filter for which calculations are done
							with their names provided by filter_names
						
						min_UL : `float`
							Min upper limit mag for cosndieration.

						max_tsep : `float`
							Max obs time difference for "simulataneus" calc.
 						
						max_tgap : `float`
							Max time between detections for no gaps
			Returns
			-------
				dict with entries as in class doc string.
					
					{
						'ndet' : 3,
						...
					}
		"""
		
		# parse some run config params and use default if not given
		lc_filter = run_config.get('lc_filter')
		if lc_filter is None:
			lc_filter = self.default_lc_filter
		filter_ids = run_config.get('filter_ids')
		if filter_ids is None:
			filter_ids = self.default_filter_ids
		filter_names = run_config.get('filter_names')
		if filter_names is None:
			filter_names = self.default_filter_names
		min_UL =  run_config.get('min_UL')
		if min_UL is None:
			min_UL = self.default_min_UL
		max_tsep =  run_config.get('max_tsep')
		if max_tsep is None:
			max_tsep = self.default_max_tsep
		max_tgap =  run_config.get('max_tgap')
		if max_tgap is None:
			max_tgap = self.default_max_tgap

		# Output dict that we will start to populate
		o = {}

		### Step 1. Base determinations based on combined detections
		self.logger.debug('Starting joint band RiseDeclineStat estimations')

		# Detection photo-points
		pps = light_curve.get_ntuples(['obs_date', 'fid','mag','magerr','rb'], filters=lc_filter)
		# Check if no datapoints fulfill initial lc criteria
		if len(pps)==0:
			return {'success':False, 'cause':'No data survive selection criteria'}


		# Get observation times for non-detection limits
		ulfilter = [ {'attribute':'diffmaglim', 'operator':'>=', 'value':min_UL} ]
		uls = light_curve.get_tuples('obs_date','diffmaglim', filters=ulfilter, upper_limits=True)


		# cast to tables for convenience
		try:
			dets = Table(rows = pps, names=('jd', 'filter', 'mag', 'magerr','rb'))
			o['cut_pp']=0
			if self.del_duplicate_rows:
				unique_jd, counts = np.unique(dets['jd'],return_counts=True)
				double_jd = list(unique_jd[ (counts>1)])
				if not len(double_jd)==0:
					self.logger.info('Cuting duplicate jd photopoints at %s'%(double_jd) )
					for jd in double_jd:
						bDoubleJD = ( dets['jd']==jd )
						MaxRb = dets['rb'][bDoubleJD].max()
						bCut = ( dets['rb'][bDoubleJD] != MaxRb )
						iCut = np.arange(dets['jd'].size)[ bDoubleJD ][bCut]
						o['cut_pp'] += len(iCut)
						dets.remove_rows( iCut )
		
		except ValueError:
			dets = Table(rows = pps, names=('jd', 'filter', 'mag', 'magerr'))


		
		# First set of properties to derive
		o['jd_det'] = dets['jd'].min()
		o['jd_last'] = dets['jd'].max()
		o['ndet'] = dets['jd'].size
		try:
			o['mag_det'] = float( dets['mag'][ dets['jd']==o['jd_det'] ] )
		except TypeError:
			o['mag_det'] = float( dets['mag'][ dets['jd']==o['jd_det'] ] )

		o['mag_last'] = float( dets['mag'][ dets['jd']==o['jd_last'] ] )
		o['t_lc'] = o['jd_last'] - o['jd_det']

		# Check if (d)real bogus present for any of these + host values
		for proptype in ['rb','drb','distnr','magnr','classtar','sgscore1','distpsnr1','sgscore2','distpsnr2','neargaia','maggaia']:
			propvalues = list( filter(None, light_curve.get_values(proptype, filters=lc_filter) ) )
			if len(propvalues)>0:
				o['%s_med'%(proptype)] = np.nanmedian( propvalues )
			else:
				o['%s_med'%(proptype)] = None



		# Look at upper limits
		if len(uls)>0:
			ulims = Table(rows = uls, names=('jd','diffmaglim'))
			# Check for presence of upper limits after first detection 
			if ( ulims['jd']>o['jd_det'] ).sum() > 0:
				o['bool_pure'] = False
			else:
				o['bool_pure'] = True
			# Latest upper limit prior to detection
			if np.any((ulims['jd']<o['jd_det'])):
				o['t_predetect'] = o['jd_det'] - ulims['jd'][ (ulims['jd']<o['jd_det']) ].max()
			else:
				o['t_predetect'] = None
		else:
			ulims = Table(( [None], [None] ),names=('jd','diffmaglim'))
			o['bool_pure'] = True
			o['t_predetect'] = None
			# Will ignore upper limits from now

		# Has the lightcurve peaked?
		# Requires the most recent detection to be significantly fainter than the brightest one
		# and more than max_tsep to have gone since peak light
		min_mag = dets['mag'].min()
		min_mag_jd = float( dets['jd'][ dets['mag']==min_mag ] )
		min_mag_err = float( dets['magerr'][ dets['mag']==min_mag ] )
		try:
			last_mag_err = float( dets['magerr'][ dets['jd']==o['jd_last'] ] )
		except TypeError:
			print('last mag err')
			print(dets)
			print(o)
			print(dets['mag']==o['mag_last'])
			last_mag_err = float( dets['magerr'][ dets['mag']==o['mag_last'] ] )

		det_mag_err = float( dets['magerr'][ dets['jd']==o['jd_det'] ] )
		if (o['jd_last']-min_mag_jd)<max_tsep:
			self.logger.info('Latest detection too close to peak light to calculate peak stats')
			o['bool_peaked'] = False
		else:
			if (o['mag_last']-min_mag)>np.sqrt( min_mag_err**2 + last_mag_err**2 ):
				o['bool_peaked'] = True
			else:
				o['bool_peaked'] = False

	
		# If we concluded a peak was there, collect info
		if o['bool_peaked']:
			self.logger.debug('Calculating peak based statistics')
			o['jd_max'] = min_mag_jd
			o['mag_peak'] = min_mag
			o['bool_rising'] = False # If it has peaked it is def not rising
			
			# Other chracteristics
			# Did it not rise at all?
			if (o['mag_det']-o['mag_peak'])>np.sqrt( det_mag_err**2 + min_mag_err**2 ):
				o['bool_norise'] = False
			else:
				o['bool_norise'] = True

		else:
			self.logger.debug('Calculating statistics assuming no peak reached')
			o['jd_max'] = None
			o['mag_peak'] = None

			# Did it not rise at all from first detection detection?
			if (o['mag_det']-o['mag_last'])>np.sqrt( det_mag_err**2 + last_mag_err**2 ):
				o['bool_norise'] = False
			else:
				o['bool_norise'] = True

			# Here it makese sense to check whether it is still rising
			if (o['jd_last']-min_mag_jd)<max_tsep:
				o['bool_rising'] = True
			else:
				o['bool_rising'] = False

		# Are there long gaps among the detections?
		jdsorted = np.unique( dets['jd'] )
		if len(jdsorted)>1:
			if (jdsorted[1:]-jdsorted[0:-1]).max() > max_tgap:
				o['bool_hasgaps'] = True
			else:
				o['bool_hasgaps'] = False
		else:
				o['bool_hasgaps'] = None




		# Look at filter dependent qualities

		# Slopes
		for filtid in filter_ids:
			self.logger.debug('Starting slope fit filt %s'%(filtid))
			filter_det = dets[ dets['filter'] == filtid ]		

			if o['bool_peaked']:
				filter_det_rise = filter_det[ filter_det['jd']<= o['jd_max'] ]
				filter_det_fall = filter_det[ filter_det['jd']>= o['jd_max'] ]

				# Rise
				# Check that lc had rise with sufficient detections
				if o['bool_norise'] or filter_det_rise['jd'].size<2:
					o['slope_rise_%s'%(filter_names[filtid])] = None
				else:
					p = np.polyfit( filter_det_rise['jd'], filter_det_rise['mag'], 1, 
						w = 1./ filter_det_rise['magerr'] )
					o['slope_rise_%s'%(filter_names[filtid])] = p[0]
	
				# Decline
				# Only makes sense to check this if lc peaked and declined for significant time
				if filter_det_fall['jd'].size>1 and (filter_det['jd'].max()-o['jd_max'])>max_tsep:
					p = np.polyfit( filter_det_fall['jd'], filter_det_fall['mag'], 1, 
							w = 1./ filter_det_fall['magerr'] )
					o['slope_fall_%s'%(filter_names[filtid])] = p[0]
				else:	
					o['slope_fall_%s'%(filter_names[filtid])] = None
			else:
				# Will use all the data to fit rise parameter, set others to none
				if o['bool_norise'] or filter_det['jd'].size<2:
					o['slope_rise_%s'%(filter_names[filtid])] = None
				else:
					p = np.polyfit( filter_det['jd'], filter_det['mag'], 1, 
						w = 1./ filter_det['magerr'] )
					o['slope_rise_%s'%(filter_names[filtid])] = p[0]

		# Colors at specific phases
		for coljd, colname in zip( [o['jd_det'], o['jd_last'], o['jd_max']], ['col_det','col_last','col_peak'] ):
			self.logger.debug('Checking col %s at jd %s'%(colname, coljd))
			# Check if time defined (e.g. if peak not known)
			if coljd is None:
				o[colname] = None
				continue  
			dets_attime = dets[ np.abs(dets['jd']-coljd)<=max_tsep ]
			if  np.any( dets_attime['filter']==filter_ids[0] ) and  np.any( dets_attime['filter']==filter_ids[1] ):
				col = np.mean( dets_attime['mag'][dets_attime['filter']==filter_ids[0]] )
				col -= np.mean( dets_attime['mag'][dets_attime['filter']==filter_ids[1]] )
				if np.isnan(col):
					col = None
				o[colname] = col
			else:
				o[colname] = None
			
		## This is where the data collection stops and evaluation starts. 
		## Even though not all of the properties above are used we keep them for future compatibility.

		# Did not train for gap objects or objects with intervening upper limits
		if o['bool_hasgaps'] or not o['bool_pure']:
			self.logger.info("Abort due to LC quality",
				extra={'bool_hasgaps':o['bool_hasgaps'],'bool_pure':o['bool_pure']})
			o['SNGuess'] = 0.5   # Model could not run
			o['SNGuessBool'] = None
			o['success'] = False
			return o


		# Determine which lightcurve requirements to use
		if o['ndet']<2 or o['ndet']>100:
			self.logger.info("LC size not covered", extra={'nbr_det':o['ndet']})
			o['SNGuess'] = 0.5   # Model could not run
			o['SNGuessBool'] = None
			o['success'] = False
			return o
	
		# Which tree param key to use?
		if o['ndet']<7:
			treekey = o['ndet']
		elif o['ndet']<=100:
			treekey = 100

		# One large bin for all longer lightcurves at this time
		max_lc_time = self.xgb_tree_param[ treekey ]['max_duration']
		max_predetect_time = self.xgb_tree_param[ treekey ]['max_predetect']
		min_det_mag = self.xgb_tree_param[ treekey ]['min_detmag']

		# Verify lightcurve properties
		if o['t_lc'] >max_lc_time or o['t_predetect'] is None or o['t_predetect']>max_predetect_time or o['mag_det']<min_det_mag:
			self.logger.info("LC prop outside model training", extra={'t_lc':o['t_lc'],'t_predetect':o['t_predetect'],'mag_det':o['mag_det']} )
			o['SNGuess'] = 0.5   # Model could not run
			o['SNGuessBool'] = None
			o['success'] = False
			return o


		# The xgb tree cannot handle None, but np.nan works
		nandict = {k:(np.nan if v is None else v) for (k,v) in o.items()}
		fitresult = self.xgb_tree.predict_risedecline( nandict, treekey)
		o['SNGuess'] = fitresult

		if fitresult>0.5:
			o['SNGuessBool'] = 1
		else:
			o['SNGuessBool'] = 0
		
		# In either case the fit succeeded.
		o['success'] = True
		return o
