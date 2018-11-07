#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2CatalogMatch.py
# License           : BSD-3-Clause
# Author            : matteo.giomi@desy.de
# Date              : 24.08.2018
# Last Modified Date: 24.08.2018
# Last Modified By  : matteo.giomi@desy.de

import logging
from pymongo import MongoClient
from astropy.coordinates import SkyCoord
from astropy.table import Table
from urllib.parse import urlparse

from ampel.base.abstract.AbsT2Unit import AbsT2Unit
from ampel.core.flags.T2RunStates import T2RunStates

from extcats import CatalogQuery
from extcats.catquery_utils import get_closest
from numpy import asarray
import zerorpc

class T2CatalogMatch(AbsT2Unit):
	"""
		cross match the position of a transient to those of sources in a set
		of catalogs and attach the required information to the transient.
	"""
	
	version = 0.1
	resources = ('extcats.reader', 'catsHTM.default')

	def __init__(self, logger, base_config):
		"""
		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		
		# empty dict of suppoerted (AS WELL AS REQUESTED) extcats catalog query objects
		self.catq_objects = {}
		
		# initialize the catsHTM paths and the extcats query client.
		self.catshtm_client 			= zerorpc.Client(base_config['catsHTM.default'])
		self.catq_client 			= MongoClient(base_config['extcats.reader'])
		self.catq_kwargs_global 		= {
										'logger': self.logger,
										'dbclient': self.catq_client,
										'ra_key': 'ra',
										'dec_key': 'dec'
										}
		
		# default parameters for LightCurve.get_pos method
		self.lc_get_pos_defaults = {'ret': "brightest", 'filters': None}
		
		# mandatory keys
		self.mandatory_keys = ['use', 'rs_arcsec']

	def init_extcats_query(self, catalog, catq_kwargs=None):
		"""
			Return the extcats.CatalogQuery object corresponding to the desired
			catalog. Repeated requests to the same catalog will not cause new duplicated
			CatalogQuery instances to be created.
			
			Returns:
			--------
				
				extcats.CatalogQuery instance.
			
		"""
		
		# check if the catalog exist as an extcats database
		if not catalog in self.catq_client.database_names():
			raise ValueError("cannot find %s among installed extcats catalogs"%(catalog))
		
		# check if you have already init this peculiar query
		catq = self.catq_objects.get(catalog)
		if catq is None:
			self.logger.debug("CatalogQuery object not previously instantiated. Doing it now.")
			
			# add catalog specific arguments to the general ones
			if catq_kwargs is None:
				catq_kwargs = self.catq_kwargs_global
			else:
				merged_kwargs = self.catq_kwargs_global.copy()
				merged_kwargs.update(catq_kwargs)
			self.logger.debug("Using arguments: %s", merged_kwargs)
			
			# init the catalog query and remember it
			catq = CatalogQuery.CatalogQuery(catalog, **merged_kwargs)
			self.catq_objects[catalog] = catq
			return catq
		else:
			self.logger.debug("CatalogQuery object for catalog %s already exists."%catalog)
			return catq

	def run(self, light_curve, run_config):
		""" 
			Parameters
			-----------
				light_curve: "ampel.base.LightCurve" instance. 
				 See the LightCurve docstring for more info.
			
				run_parameters: `dict`
						configuration parameter for this job. There is provision to
						pass arguments to the LightCurve.get_pos method used to derive
						the position of the transient from the lightcure. 
						
						Most importantly, the catalogs key correspond to a nested dictionary
						in which each entry specify a catalog in extcats or catsHTM format 
						and the parameters used for the queries. 
						
						for each entry in the 'catalogs' configuration dictionary, the 
						following keys are MANDATORY:
							
							'use': `str`
								either extcats or catsHTM, depending on how the catalog is set up.
							
							'rs_arcsec': `float`
								search radius for the cone search, in arcseconds
							
						Optional keys includes:
							
						'catq_kwargs': `dict`
							parameters to pass to the catalog query routine. Two cases arise:
								if 'use' == 'extcats':
									the 'catq_kwargs' can (or MUST I don't ) contain the names of the ra and dec
									keys in the catalog (se example below), all valid arguments to
									extcats.CatalogQuert.findclosest can be given, such as pre- and post 
									cone-search query filters can be passed.
								if 'use' == 'catsHTM':
									the 'catq_kwargs' SHOULD contain the the names of the ra and dec
									keys in the catalog if those are different from 'ra' and 'dec'
						
						the 'keys_to_append' parameters is OPTIONAL and specifies wich fileds from 
							the catalog should be returned in case of positional match:
								
								if not present or if equal to 'all':
									all the fields in the given catalog will be returned.
								if `list`
									just take this subset of fields.
						
						Eg:
				
						run_config = 
							{
							'get_lc_pos_kwargs': None, # optional see ampel.base.LightCurve doc
							'catalogs':
								{
								'sdss_spec': 
									{
										'use': 'extcats',
										'catq_kwargs': {
														'ra_key': 'ra',
														'dec_key': 'dec'
													},
										'rs_arcsec': 3,
										'keys_to_append': ['z', 'bptclass', 'subclass', 'dist2alert']
									},
								'NED': 
									{
										'use': 'catshtm',
										'rs_arcsec': 20,
										'keys_to_append': ['fuffa1', 'fuffa2', ..],
										'catq_kwargs': 
									},
								...
								}
							}
			
			Returns
			-------
				dict with the keys to append to each transient. For example, with the
				above run-config (and in case of a match in SDSS but not in NED), the 
				returned dict would be:
				
				{
					'SDSS_spec': 
						{
						'z': 0.08820018172264099, 
						'bptclass': 2.0, 
						'subclass': '', 
						'dist2transient': 1.841666956181802e-09}
						},
					'NED': False
				}
				
				Note that, when a match is found, the distance of the lightcurve object
				to the catalog counterpart is also returned as the 'dist2transient' key.
		"""
		
		# get ra and dec from lightcurve object
		lc_get_pos_kwargs = run_config.get('lc_get_pos_kwargs')
		if lc_get_pos_kwargs is None:
			lc_get_pos_kwargs = self.lc_get_pos_defaults
		self.logger.debug("getting transient position from lightcurve using args: %s", lc_get_pos_kwargs)
		transient_ra, transient_dec = light_curve.get_pos(**lc_get_pos_kwargs)
		self.logger.debug("Transient position (ra, dec): %.4f, %.4f deg"%(transient_ra, transient_dec))
		
		# initialize the catalog quer(ies). Use instance variable to aviod duplicates
		out_dict = {}
		catalogs = run_config.get('catalogs')
		for catalog, cat_opts in catalogs.items():
			src, dist = None, None
			self.logger.debug("Loading catalog %s using options: %s"%(catalog, str(cat_opts)))
			
			# check options:
			for opt_key in self.mandatory_keys:
				if not opt_key in cat_opts.keys():
					message = ("options for catalog %s are missing mandatory %s argument. Check your run config."%
						(catalog, opt_key))
					raise KeyError(message)
			
			# how do you want to support the catalog?
			use = cat_opts.get('use')
			if use == 'extcats':
				
				# get the catalog query object and do the query
				catq = self.init_extcats_query(catalog, catq_kwargs=cat_opts.get('catq_kwargs'))
				src, dist = catq.findclosest(transient_ra, transient_dec, cat_opts['rs_arcsec'])
				
			elif use == 'catsHTM':
				
				# catshtm needs coordinates in radians
				transient_coords = SkyCoord(transient_ra, transient_dec, unit='deg')
				srcs, colnames, colunits = self.catshtm_client.cone_search(
													catalog,
													transient_coords.ra.rad, transient_coords.dec.rad,
													cat_opts['rs_arcsec'])
				if len(srcs) > 0:
					
					# format to astropy Table
					srcs_tab = Table(asarray(srcs), names=colnames)
					
					# find out how ra/dec are called in the catalog
					catq_kwargs = cat_opts.get('catq_kwargs')
					if catq_kwargs is None:
						ra_key, dec_key = 'ra', 'dec'
					else:
						ra_key, dec_key = catq_kwargs.get('ra_key', 'ra'), catq_kwargs.get('dec_key', 'dec')
					
					# get the closest source and its distance (catsHTM stuff is in radians)
					src, dist = get_closest(transient_coords.ra.rad, transient_coords.dec.rad, srcs_tab, ra_key, dec_key)
			else:
				message = "use option can not be %s for catalog %s. valid are 'extcats' or 'catsHTM'"%(use, catalog)
				raise ValueError(message)
			
			# now add the results to the output dictionary
			out_dict_catalog = {}
			if not src is None:
				self.logger.debug("found counterpart %.2f arcsec away from transient."%dist)
				# if you found a cp add the required field from the catalog:
				# if keys_to_append argument is given or if it is equal to 'all'
				# then take all the columns in the catalog. Otherwise only add the 
				# requested ones.
				keys_to_append = cat_opts.get('keys_to_append', 'all')
				if keys_to_append == 'all':
					keys_to_append = src.colnames

				out_dict[catalog] = {field: src[field] for field in keys_to_append}
				out_dict[catalog]['dist2transient'] = dist
			else:
				self.logger.debug("no match found in catalog %s within %.2f arcsec from transient"%
					(catalog, cat_opts['rs_arcsec']))
				out_dict[catalog] = False
			
		# return the info as dictionary
		return out_dict
