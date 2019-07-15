#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File				: ampel/contrib/hu/t3/T3MarshalMonitor
# License			: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date				: 17.11.2018
# Last Modified Date: 06.02.2019
# Last Modified By	: matteo.giomi@desy.de

import re
import logging
import numpy as np
from pydantic import BaseModel, BaseConfig
from typing import Dict, List

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.base.dataclass.JournalUpdate import JournalUpdate
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils

from ampel.contrib.hu.t3.ampel_tns import sendTNSreports, get_tnsname, TNSFILTERID, tnsInternal

# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
	# For item i in a range that is a length of l,
	for i in range(0, len(l), n):
		# Create an index range for l of n items:
		yield l[i:i+n]

# get the science records for the catalog match
def get_catalogmatch_srecs(tran_view, logger):
	cat_res = tran_view.get_science_records(t2_unit_id="CATALOGMATCH")
	if len(cat_res) == 0 or cat_res is None or cat_res[-1].get_results() is None:
		logger.info("NO CATALOG MATCH FOR THIS TRANSIENT")
		return {}
	return cat_res[-1].get_results()[-1]['output']


class TNSTalker(AbsT3Unit):
	"""
		Get TNS name if existing, and submit selected candidates
	"""
	
	version = 0.1
	
	class RunConfig(BaseModel):
		"""
 		Necessary class to validate configuration.
		"""
		class Config(BaseConfig):
			"""
			Raise validation errors if extra fields are present
			"""
			allow_extra = False
			ignore_extra = False
		
		# TNS config
		tns_api_key			: str	= None		# Bot api key frm TNS
		get_tns_force		: bool	= False		# Check for TNS for names even if internal name is known
		submit_tns	 		: bool	= True		# Submit candidates passing criteria (False gives you a 'dry run')
		resubmit_tns_nonztf	: bool	= True		# Resubmit candidate submitted w/o the same ZTF internal ID 
		resubmit_tns_ztf	: bool	= False		# Resubmit candidates even if they have been added with this name before
		sandbox				: bool	= True		# Submit to TNS sandbox only
		ext_journal			: bool	= True		# weather journal will go to separate collection.
		
		# AT report config
		base_at_dict		: Dict = {
			"groupid":"48", 
			"reporter": "J. Nordin, V. Brinnel, M. Giomi, J. van Santen (HU Berlin), A. Gal-Yam, O. Yaron, S. Schulze (Weizmann) on behalf of ZTF",
			"at_type":"1"
		}
		ztf_tns_at			: Dict = {	# Default values to tag ZTF detections / ulims
										"flux_units":"1",
										"instrument_value":"196",
										"exptime":"30",
										"Observer":"Robot"
								}
		max_maglim			: float = 20	# Limiting magnitude to consider upper limits as 'significant'
		nphot_submit		: int	= 2		# Number of photometric detection we include in the TNS AT report
		
		# cuts on T2 catalogs
		needed_catalogs	: List[str]	= []# reject candidates if they don't have matching in this list of T2CATALOGMATCH catalogs
		require_catalogmatch : bool = True
		max_redshift	: float	= 1.15	# maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_redshift	: float	= 0		# minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		start_dist		: float = 1.5 	# arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
		max_gaia_neighbour_gmag : float = 11 # reject transient if the GAIA source brighter than this is nearby.
		
		# cut on alert properties
		min_ndet		: int	= 2		# A candidate need to have at least this many detections
		min_ndet_postul	: int	= 2		# and if it has this minimum nr of detection after the last significant (max_maglim) UL.
		max_age			: float = 5		# days, If a detection has an age older than this, skip (stars,age).
		min_age			: float = 0		# Min age of detection history
		min_peak_mag	: float	= 19.5	# range of peak magnitudes for submission
		max_peak_mag	: float = 13	#
		min_n_filters	: int	= 1		# Reported detections in at least this many filters
		min_gal_lat		: float = 14	# Minimal galactic latitide
		ssdistnr_max	: float = 1		# reject alert if ssdistnr smaller than this value for any pp
		ps1_sgveto_rad	: float = 1		# reject alert if PS1 star for any pp
		ps1_sgveto_sgth	: float = 0.8	# 
		rb_minmed		: float = 0.3	# Minimal median RB.
		
		
		# Cut to apply to all the photopoints in the light curve.
		# This will affect all operations, i.e. evaluating the position, 
		# computing number of detections ecc.
		lc_filters		: List[Dict]= [
										{
										'attribute': 'sharpnr',
										'operator': '>=', 
										'value': -10.15
										}, 
										{
										'attribute': 'magfromlim',
										'operator': '>',
										'value': 0
										}
									]
		
		# parameters for adding remarks to AT reports
		nuclear_dist	: float = 1.	# Tag objects this close to SDSS galaxies as nuclear. Use negative to disable
		aav_dist		: float = 1.	# Required distance to match with aav catalog. TODO: move?
		max_gaia_noise	: float = 2.	# (sigma!) if GAIA match is noisier than this, add a remark

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		self.run_config = self.RunConfig() if run_config is None else run_config
		self.name = "TNSTalker"
		self.logger.info("Initialized T3 TNSTalker instance %s"%self.name)
		self.logger.info("base_config: %s"%self.base_config)
		self.logger.info("run_config: %s"%self.run_config)
		
#		# TODO: do we want to leave this
#		if self.api_key is None:
#			#raise KeyError("No TNS api_key, cannot run.")
#			self.logger.info("No TNS api_key, using default + sandbox")	
#			self.api_key = "a3f9bcbbe6a26a1ae97a0eeefe465932e68cba83"
#			self.sandbox = True

	def search_journal_tns(self, tran_view):
		"""
			Look through the journal for a TNS name.
                        Assumes journal entries came from this unit, that the TNS name is saved as "tnsName"
                        and internal names as "tnsInternal"
		"""
		# Find the latest tns name (skipping previous)
		jentry = tran_view.get_journal_entries(
			filterFunc=lambda x: x.get('t3unit')==self.name and 'tnsInternal' in x.keys(),
			latest=True)
		tns_name = None if jentry is None else jentry.get("tnsName", None)

		# Find internal names
		jentries = tran_view.get_journal_entries(
			filterFunc=lambda x: x.get('t3unit')==self.name and 'tnsInternal' in x.keys())
		tns_internals = [] if jentries is None else [j.get('tnsInternal',None) for j in jentries]
		self.logger.info('Journal search',extra={'tranId':tran_view.tran_id,'tnsName':tns_name,'tnsInternals':tns_internals})

		return tns_name, tns_internals

	def _query_tns_names(self, tran_view):
		"""
			query the TNS for names and internals at the position
			of the transient.
		"""
		# query the TNS for transient at this position. Note that we check the real TNS for names for compatibility...
		ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean", filters=self.run_config.lc_filters)
		tns_name, tns_internals = get_tnsname(
							ra=ra, dec=dec, 
							api_key=self.run_config.tns_api_key, 
							logger=self.logger, 
							sandbox=False
						)
		
		# Skip the AT SN prefix if present
		if tns_name is not None:
			tns_name = re.sub('^AT','',tns_name)
			tns_name = re.sub('^SN','',tns_name)
		
		# be nice and then go
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		self.logger.info(
			"looking for TNS name in the TNS.",
			extra={
				'ZTFname': ztf_name, 
				'ra': ra,
				'dec': dec,
				'tnsName': tns_name,
				'tnsInternals': tns_internals,
			})
		return tns_name, tns_internals


	def _find_tns_tran_names(self, tran_view):
		"""
			search for TNS name in tran_view.tran_names. If found, 
			look in the TNS for internal names and return them
		"""
		
		tns_name, tns_internals = None, []
		for tname in tran_view.tran_names:

			if 'TNS' in tname and (not self.run_config.get_tns_force):
				self.logger.info(
					"found TNS name in tran_names.", extra={'TNSname': tname, 'TransNames': tran_view.tran_names})
				# as TNS to give you the internal names.
				# we remove the 'TNS' part of the name, since this has been 
				# added by the TNSMatcher T3, plus we skip the prefix
				tns_name = tname.replace("TNS", "")	# We here assume that the AT/SN suffix is cut				
				# Not using sandbox (only checking wrt to full system)
				tns_internals, runstatus = tnsInternal(
								tns_name,   
								api_key=self.run_config.tns_api_key,
								sandbox=False
							)

		# be nice with the logging
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		self.logger.info(
			"looked for TNS name in self.tran_names",
			extra={
				'ZTFname': ztf_name,
				'tnsName': tns_name,
				'tnsInternals': tns_internals,
				'TransNames': tran_view.tran_names
			})
		
		#if you make it till here, no match was found
		return tns_name, tns_internals

	def find_tns_name(self, tran_view):
		"""
			extensive search for TNS names in:
				- tran_view.tran_names (if added by TNSMatcher)
				- the journal of tran_view (if added by this T3)
				- the TNS itself (if no name can be found with the above)
			
			Returns:
			--------
				tns_name, tns_internals, jup: tns_name, tns_internal, and journal update
		"""
		
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		self.logger.info("looking for TNS name for transient.", extra={'ZTFname': ztf_name})
		
		# first we look in the journal, this is the cheapest option. If we have 
		# a valid name from the journal and if you do not want to look again in 
		# the TNS, we are fine. NOTE: in this case you don't return a journal update.
		tns_name, tns_internals = self.search_journal_tns(tran_view)
		if (not tns_name is None) and (not self.run_config.get_tns_force):
			return tns_name, tns_internals, None
		
		# second option in case there is no TNS name in the journal: go and look in tran_names
		# and if you don't find any, go and ask TNS again.
		tns_name_new, tns_internals_new = self._find_tns_tran_names(tran_view)
		if tns_name_new is None:
			tns_name_new, tns_internals_new = self._query_tns_names(tran_view)
		
		# now, it is possible (if you set self.run_config.get_tns_force) that the
		# new TNS name is different from the one we had in the journal. We always
		# use the most recent one. In this case we also create a JournalUpdate 
		jup = None
		if not tns_name_new is None:
			
			# what happen if you have a new name that is different from the old one?
			if tns_name is not None and not tns_name==tns_name_new:
				self.logger.info("Adding new TNS name",extra={"tnsOld":tns_name,"tnsNew":new_tns_name})
			
				# create content of journal entry. Eventually 
				# update the list with the new internal names if any are found
				jcontent = {'t3unit': self.name, 'tnsName': new_tns_name}
				if tns_internals_new is not None:
					tns_internals.append(tns_internals_new)
					for tns_int in tns_internals_new:
						jcontent.update({'tnsInternal':tns_int})
				
				# create a journalUpdate and update the tns_name as well. TODO: check with JNo
				jup = JournalUpdate(tran_id=tran_view.tran_id, ext=self.run_config.ext_journal, content=jcontent)
				tns_name = tns_name_new
				tns_internals = tns_internals_new
		
		# bye!
		return tns_name, tns_internals, jup


	def accept_tview(self, tran_view):
		"""
			decide weather or not this transient is worth submitting to TNS 
			or not. 
			
			NOTE that even if many of these cuts could defined passed directly to
			the task/job config, some of them still require relatively non trivial
			computation (e.g. 'age' of the transient). This makes this selection method 
			necessary.
			
			FOR SAKE OF SIMPLICITY during the first period of TNS submission, 
			we here declare that all the transient selection logic should be 
			implemented in this method, even rejection based on catalog matching 
			that could be done by the select config param.
		"""
		
		# get the latest light curve
		lc = tran_view.get_latest_lightcurve()
		
		# apply cut on history: consider photophoints which are sharp enough
		pps = lc.get_photopoints(filters=self.run_config.lc_filters)
		self.logger.info("%d photop. passed filter %s"%(len(pps), self.run_config.lc_filters))
		
		# cut on number of detection
		if len(pps) < self.run_config.min_ndet:
			self.logger.info("not enough detections: got %d, required %d"%
				(len(pps), self.run_config.min_ndet))
			return False
		
		# cut on number of detection after last SIGNIFICANT UL
		ulims = lc.get_upperlimits(filters={'attribute':'diffmaglim', 'operator':'>=', 'value':self.run_config.max_maglim})
		if len(ulims) > 0:
			last_ulim_jd = sorted(ulims, key=lambda x: x.get_value('jd'))[-1].get_value('jd')
			pps_after_ndet = lc.get_photopoints( 
				filters = self.run_config.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_ulim_jd}])
			if len(pps_after_ndet) < self.run_config.min_ndet_postul:
				self.logger.info("not enough consecutive detections after last significant UL.",
					extra={'NDet': len(pps), 'lastUlimJD': last_ulim_jd})
				return False
		
		# cut on number of filters
		used_filters = set([pp.get_value('fid') for pp in pps])
		if len(used_filters) < self.run_config.min_n_filters:
			self.logger.info("requested detections in more than %d bands, got: %d"%
				(self.run_config.min_n_filters, len(used_filters)))
			return False
		
		# cut on range of peak magnitude
		mags = [pp.get_value('magpsf') for pp in pps]
		peak_mag = min(mags)
		if peak_mag > self.run_config.min_peak_mag or peak_mag < self.run_config.max_peak_mag:
			self.logger.info("peak magnitude of %.2f outside of range [%.2f, %.2f]"%
				(peak_mag, self.run_config.min_peak_mag, self.run_config.max_peak_mag))
			return False
		
		# cut on age
		jds = [pp.get_value('jd') for pp in pps]
		most_recent_detection, first_detection = max(jds), min(jds)
		#age = Time.now().jd - min(jds)
		age = most_recent_detection - first_detection
		if age > self.run_config.max_age or age < self.run_config.min_age:
			self.logger.info("age of %.2f days outside of range [%.2f, %.2f]"%
				(age, self.run_config.min_age, self.run_config.max_age))
			return False
		
		# cut on galactic coordinates
		ra, dec = lc.get_pos(ret="mean", filters=self.run_config.lc_filters)
		coordinates = SkyCoord(ra, dec, unit='deg')
		b = coordinates.galactic.b.deg
		if abs(b) < self.run_config.min_gal_lat:
			self.logger.info("transient at b=%.2f too close to galactic plane (cut at %.2f)"%
				(b, self.run_config.min_gal_lat))
			return False
		
		# cut on distance to closest solar system object
		# TODO: how to make this check: ('0.0' in list(phot["ssdistnr"])
		ssdist = np.array([pp.get_value('ssdistnr') for pp in pps])
		ssdist[ssdist==None] = -999
		#print (ssdist)

		close_to_sso = np.logical_and(ssdist < self.run_config.ssdistnr_max, ssdist > 0)
		if np.any(close_to_sso):
			self.logger.info("transient too close to solar system object", extra={'ssdistnr': ssdist.tolist()})
			return False
		
		# check PS1 sg for the full alert history
		distpsnr1, sgscore1 = zip(*lc.get_tuples('distpsnr1', 'sgscore1', filters=self.run_config.lc_filters))
		is_ps1_star = np.logical_and(
									np.array(distpsnr1) < self.run_config.ps1_sgveto_rad,
									np.array(sgscore1) > self.run_config.ps1_sgveto_sgth
								)
		if np.any(is_ps1_star):
			self.logger.info(
				"transient below PS1 SG cut for at least one pp.",
				extra={'distpsnr1': distpsnr1, 'sgscore1': sgscore1}
				)
			return False
		
		# cut on median RB score
		rbs = [pp.get_value('rb') for pp in pps]
		if np.median(rbs) < self.run_config.rb_minmed:
			self.logger.info(
				"Median RB %below limit.",
				extra={'median_rd': np.median(rbs), 'rb_minmed': self.run_config.rb_minmed}
				)
			return False
		
		# ----------------------------------------------------------------------#
		# 							CUTS ON T2 RECORDS							#
		# ----------------------------------------------------------------------#
		cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)
		
		# check that we got any catalogmatching results (that it was run)
		if self.run_config.require_catalogmatch and len(cat_res)==0:
			self.logger.info("no T2CATALOGMATCH results")
			return False

		# check that you have positive match in all of the necessary cataslogs:
		for needed_cat in self.run_config.needed_catalogs:
			if not cat_res.get(needed_cat, False):
				self.logger.info("no T2CATALOGMATCH results for %s"%needed_cat, extra={'catalog_matches': cat_res})
				return False
		
		nedz		= cat_res.get('NEDz', False)
		sdss_spec 	= cat_res.get("SDSS_spec", False)
		if ((nedz and not (self.run_config.min_redshift <  nedz['z'] < self.run_config.max_redshift)) or 
			(sdss_spec and not (self.run_config.min_redshift < sdss_spec['z'] > self.run_config.max_redshift))):
			self.logger.info("transient z above limit.", extra={'max_z': self.run_config.max_redshift, 'SDSSspec': sdss_spec, 'NEDz': nedz})
			return False
		
		# another battle in the endless war against stars.
		# here we define a dict to treat each catalog in the same way
		star_filters = {
			'SDSSDR10':		{'class_col': 'type', 'star_val': 6},
			'LAMOSTDr4':	{'class_col': 'class', 'star_val': 'STAR'},
			}
		for cat_name, sfilter in star_filters.items():
			cat = cat_res.get(cat_name, False)
			cname, sval = sfilter['class_col'], sfilter['star_val']
			if cat and cat[cname] == sval and cat['dist2transient'] < self.run_config.start_dist:
				self.logger.info("transient matched with star in catalog.", extra={'cat_name': cat_name, 'cat_res': cat})
				return False
		
		# cut matches with variable star catalog
		aavsovsx = cat_res.get('AAVSOVSX', False)
		if aavsovsx and aavsovsx['dist2transient'] < self.run_config.start_dist:
			self.logger.info("transient too close to AAVSOVSX sorce", extra=aavsovsx)
			return False
		
		# cut away bright stars. TODO: this considers just the closest matches...
		gaia_dr2 = cat_res.get('GAIADR2', None)
		if gaia_dr2 and gaia_dr2['Mag_G'] > 0 and gaia_dr2['Mag_G'] < self.run_config.max_gaia_neighbour_gmag:
			self.logger.info("transient close to bright GAIA source", extra=gaia_dr2)
			return False
		
		# congratulation TransientView, you made it!
		return True


	def add_atreport_remarks(self, tran_view):
		"""
			create additional remarks based, i.e. on catalog matching data
		"""
		
		# TODO: check values for the cuts and and put them in the config
		# TODO: can remarks be combined? e.g. nucler + noisy? 
		
		# get the science records for the catalog match
		cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)
		
		# tag AGNs
		milliquas = cat_res.get('milliquas', False)
		sdss_spec = cat_res.get('SDSS_spec', False)
		if (milliquas and milliquas['redshift'] > 0) or (sdss_spec and sdss_spec['bptclass'] in [4, 5]): #TODO: add distance cut?
			self.logger.info("Transient is SDSS BPT or Milliquas AGN.",
				extra={"tranId":tran_view.tran_id, 'milliquas': milliquas, 'SDSS_spec': sdss_spec})
			return {
					"remarks": "Known SDSS and/or MILLIQUAS QSO/AGN. ",
					"at_type": 3
				}
		
		# tag nuclear
		sdss_dr10 = cat_res.get('SDSSDR10', False)
		if sdss_dr10 and sdss_dr10['type'] == 3 and sdss_dr10['dist2transient'] < self.run_config.nuclear_dist:
			self.logger.info("Transient close to SDSS photometric galaxy - possibly nuclear",
				extra={"tranId":tran_view.tran_id, 'SDSSDR10': sdss_dr10})
			return {
					"remarks": "Close to core of SDSS DR10 galaxy",
					"at_type": 4
				}
		
		# tag noisy gaia
		lc = tran_view.get_latest_lightcurve()
		distpsnr1, sgscore1 = zip(*lc.get_tuples('distpsnr1', 'sgscore1', filters=self.run_config.lc_filters))
		galaxylike_ps1 = np.logical_and(np.array(distpsnr1)<1.5, np.array(sgscore1)<0.5)
		gaia_dr2 = cat_res.get('GAIADR2', False)
		nedz	 = cat_res.get('NEDz', False)
		if ( (gaia_dr2 and gaia_dr2['ExcessNoise'] > self.run_config.max_gaia_noise and gaia_dr2['dist2transient'] < 1) and 
			 (nedz and not (nedz['z']>0.01 and nedz['dist2transient'] < 1)) and 					#if it's extragalactic
			 (sdss_dr10 and not (sdss_dr10['type'] == 3 and sdss_dr10['dist2transient'] <3)) and	# and if it's not a galaxy
			 (not np.any(galaxylike_ps1))															# TODO: check the logic
			):
			self.logger.info("Significant noise in Gaia DR2 - variable star cannot be excluded.", 
				extra={"tranId":tran_view.tran_id, 'GAIADR2': gaia_dr2, 'NEDz': nedz, 'SDSSDR10': sdss_dr10})
			return {"remarks": "Significant noise in Gaia DR2 - variable star cannot be excluded." }


	def create_atreport(self,tran_view):
		"""
			Collect the data needed for the atreport. Return None in case 
			you have to skip this transient for some reason.
		"""
		
		self.logger.info("creating AT report for transient.")
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		lc = tran_view.get_latest_lightcurve()
		ra, dec = lc.get_pos(ret="mean", filters=self.run_config.lc_filters)
		
		# Start defining AT dict: name and position
		atdict = {}
		atdict.update(self.run_config.base_at_dict)
		atdict["internal_name"] = ztf_name
		atdict["ra"] = {"value": ra, "error" : 1., "units":"arcsec"}
		atdict["dec"] = {"value": dec, "error" : 1., "units":"arcsec"}
		
		# Add information on the latest SIGNIFICANT non detection. TODO: check!
		ulims = lc.get_upperlimits(filters={'attribute':'diffmaglim', 'operator':'>=', 'value':self.run_config.max_maglim})
		last_non_obs = 0
		if len(ulims) == 0:
			atdict["non_detection"] = {
								"archiveid"			:"0",
								"archival_remarks"	:"ZTF non-detection limits not available"}
		else:
			last_ulim = sorted(ulims, key=lambda x: x.get_value('jd'))[-1]
			last_non_obs = last_ulim.get_value('jd')
			filter_name = TNSFILTERID.get(last_ulim.get_value('fid'))
			atdict["non_detection"] = {
									"obsdate"		: last_ulim.get_value('jd'), 
									"limiting_flux"	: last_ulim.get_value('diffmaglim'),
									"filter_value"	: filter_name
								}
		atdict["non_detection"].update(self.run_config.ztf_tns_at)# Add the default ZTF values
		
		# now add info on photometric detections: consider only candidates which
		# have some consecutive detection after the last ulim
		pps = lc.get_photopoints(
			filters = self.run_config.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_non_obs}])
		#TODO: pps are those sorted in time?
		
		# Lets create a few photometry points: TODO: should they be the latest or the first?
		atdict["photometry"] = {"photometry_group":{}}
		atdict["discovery_datetime"] = 10**30
		for ipp, pp in enumerate(pps[:self.run_config.nphot_submit]):
			photdict = {	#TODO: do we need to round the numerical values?
				"obsdate"		: pp.get_value('jd'),
				"flux"			: pp.get_value('magpsf'),
				"flux_error"	: pp.get_value('sigmapsf'),
				"limiting_flux"	: pp.get_value('diffmaglim'),
				"filter_value"	: TNSFILTERID.get(pp.get_value('fid'))
				}
			if pp.get_value('jd')<atdict["discovery_datetime"]:
				atdict["discovery_datetime"] = pp.get_value('jd')
			photdict.update(self.run_config.ztf_tns_at)
			atdict["photometry"]["photometry_group"][ipp] = photdict
		
		# finally, add remarks based on catalogs adn return
		remarks = self.add_atreport_remarks(tran_view)
		if not remarks is None:
			atdict.update(remarks)
		return atdict

	def add(self, transients):
		"""
			Loop through transients and check for TNS names and/or candidates to submit
		"""
		
		if transients is None:
			self.logger.info("no transients for this task execution")
			return []
		
		# select the transients
		transients_to_submit = [tv for tv in transients if self.accept_tview(tv)]
		self.logger.info("of the %d transients presented to this task, %d passed selection criteria"%
			(len(transients), len(transients_to_submit)))
		
		journal_updates = []	# Will be saved to future journals
		atreports = {}			# Reports to be sent, indexed by the transient view IDs (so that we can check in the replies)
		for tran_view in transients_to_submit:
			
			ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
			self.logger.info("TNS start", extra={"tranId":tran_view.tran_id, 'ztfName': ztf_name})
			
			# find the TNS name, either from the journal, from tran_names, or 
			# from TNS itself. If new names are found, create a new JournalUpdate
			tns_name, tns_internals, jup = self.find_tns_name(tran_view)
			if not jup is None:
				journal_updates.append(jup)
			
			# Chech whether this ID has been submitted (note that we do not check 
			# whether the same candidate was submitted as different ZTF name) and
			# depending on what's already on the TNS we can chose to submit or not
			is_ztfsubmitted = ztf_name in tns_internals
			if not ( (is_ztfsubmitted and self.run_config.resubmit_tns_ztf) or 
					 (not is_ztfsubmitted and self.run_config.resubmit_tns_nonztf) ):
				self.logger.info(
					"we won't submit candidate.",
					extra = {
						'is_ztfsub': is_ztfsubmitted, 
						'tnsInternals': tns_internals
						})
				continue
			
			# create AT report
			atreport = self.create_atreport(tran_view)
			self.logger.info("Added to report list")
			atreports[tran_view.tran_id] = atreport
			
			# TODO MEGADEBUG REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			#break
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
			# REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT REMOVE IT
		
		# TODO: we save the atreports to send them to the TNS. 
		# This is just part of the tesing and will have to go away
#		atreports = {k: atreports[k] for k in list(atreports.keys())[:2]}
		self.atreports = atreports
		self.logger.info("collected %d AT reports to post"%len(atreports))
		
		# If we do not want to submit anything, or if there's nothing to submit
		if len(atreports) == 0 or (not self.run_config.submit_tns):
			self.logger.info("submit_tns config parameter is False or there's nothing to submit", 
				extra={'n_reports': len(atreports), 'submit_tns': self.run_config.submit_tns})
			return journal_updates
		
		# Send reports in chunks of size 90 (99 should work)
		atchunks = list(chunks([atr for atr in atreports.values()], 90))
		tnsreplies = sendTNSreports(atchunks, self.run_config.tns_api_key, self.logger, sandbox=self.run_config.sandbox)
		
		# Now go and check and create journal updates for the cases where SN was added
		for tran_id in atreports.keys():
			ztf_name = ZTFUtils.to_ztf_id(tran_id)
			if not ztf_name in tnsreplies.keys():
				self.logger.info("No TNS add reply",extra={"tranId":tran_id})
				continue
			
			# Create new journal entry 			#TODO: do we want to add to the journal a failed TNS submit?
			jup = JournalUpdate(
					tranId=tran_id,
					ext=self.run_config.ext_journal,
					content={
						't3unit': self.name,
						'tnsName': tnsreplies[ztf_name][1]["TNSName"],
						'tnsInternal': ztf_name,
						'tnsSubmitresult': tnsreplies[ztf_name][0]
					})
			journal_updates.append(jup)
		return journal_updates


	def done(self):
		"""
		"""
		self.logger.info("done running T3")
		
		if not hasattr(self, 'atreports'):
			self.logger.info("No atreports collected.")
			return
		
		#TODO: to help debugging and verification, we post the collected atreports
		# to the slack, so that we can compare them with what JNo script is doing
		# ALL THE CONTENT OF THIS METHOD SHOULD GO AWAY AS SOON AS WE TRUST THIS T3
		self.logger.warning("Posting collected ATreports to Slack. I'm still running as a test!")
		
		import datetime, io, json
		from slackclient import SlackClient
		from slackclient.exceptions import SlackClientError
		
		slack_token		= "xoxb-297846339667-549790069252-FLwKXkra0NL3FNnrvu9XYm4a"
		slack_channel 		= "#ampel_tns_test"
		slack_username		= "Ampel_TNS_test"
		max_slackmsg_size	= 200	# if you have more than this # of reports, send different files
		
		sc = SlackClient(slack_token)
		
		tstamp = datetime.datetime.today().strftime("%Y-%m-%d-%X")
		atlist = list(self.atreports.values())
		last = 0
		for ic, atrep in enumerate(chunks(atlist, max_slackmsg_size)):
			
			# add the atreport to a file
			self.logger.info("Posting chunk #%d"%ic)
			filename = "TNSTalker_DEBUG_%s_chunk%d.json"%(tstamp, ic)
			fbuffer = io.StringIO(filename)
			json.dump(atrep, fbuffer, indent=2)
		
			# upload the file with the at reports
			first = last
			last += len(atrep)
			msg = ("A total of %d atreports found by TNSTalker T3. Here's chunk #%d (reports from %d to %d)"%
				(len(self.atreports), ic, first, last))
			api = sc.api_call(
					'files.upload',
					channels = [slack_channel],
					title = "TNSTalker_DEBUG_%s_chunk%d"%(tstamp, ic),
					initial_comment = msg,
					username = slack_username,
					as_user = False,
					filename =  filename,
					filetype = 'javascript',
					file = fbuffer.getvalue()
				)
			if not api['ok']:
				raise SlackClientError(api['error'])
		
		self.logger.warning("DONE DEBUG Slack posting. Look at %s for the results"%slack_channel)
