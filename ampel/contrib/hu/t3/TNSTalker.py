#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File				: ampel/contrib/hu/t3/T3MarshalMonitor
# License			: BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date				: 17.11.2018
# Last Modified Date: 06.02.2019
# Last Modified By	: matteo.giomi@desy.de

import re
from pydantic import BaseModel, BaseConfig
from typing import Dict, List

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.base.dataclass.JournalUpdate import JournalUpdate
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils

from ampel.contrib.hu.t3.ampel_tns import sendTNSreports, get_tnsname, TNSFILTERID

# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
	# For item i in a range that is a length of l,
	for i in range(0, len(l), n):
		# Create an index range for l of n items:
		yield l[i:i+n]

def pedantic_get(rc, key):
	"""
		get a key from a RunConfig object, raising exception if the key is not found
	"""
	val = getattr(rc, key, 'fuffaccia')
	if val == 'fuffaccia':
		raise KeyError("parameter %s not found in the run_config. Avaliable are: %s"%
			(key, ", ".join(rc.dict().keys())))
	return val


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
		get_tns				: bool	= False		# Check for TNS names and add to journal
		get_tns_force		: bool	= False		# Check for TNS info even if internal name is known
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
		
		
		# cut on transients:
		max_redshift	: float	= 1.15	# maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_redshift	: float	= 0		# minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
		min_ndet		: int	= 2		# A candidate need to have at least this many detections
		min_ndet_postul	: int	= 2		# and if it has this minimum nr of detection after the last significant (max_maglim) UL.
		max_age			: float = 5		# days, If a detection has an age older than this, skip (stars,age).
		min_age			: float = 0		# Min age of detection history
		min_peak_mag	: float	= 19.5	# range of peak magnitudes for submission
		max_peak_mag	: float = 13	#
		min_n_filters	: int	= 1		# Reported detections in at least this many filters
		min_gal_lat		: float = 14	# Minimal galactic latitide
		needed_catalogs	: List[str]	= []# reject candidates if they don't have matching in this list of T2CATALOGMATCH catalogs
		
		# Cut to apply to all the photopoints in the light curve.
		# This will affect all operations, i.e. evaluating the position, computing number of detections ecc.
		lc_filters		: List[Dict]= [{
									'attribute': 'sharpnr',
									'operator': '>=', 
									'value': -10.15
								}]
		
		# parameters for adding remarks to AT reports
		nuclear_dist	: float = 1.	# Tag objects this close to SDSS galaxies as nuclear. Use negative to disable
		aav_dist		: float = 1.	# Required distance to match with aav catalog. TODO: move?
		max_gaia_noise	: float = 2.	# (sigma!) if GAIA match is noisier than this, add a remark
		
		#maxGaiaNoise = 2.0 # Significance!
		#minGaiaNoise = 0.0 # Significance!
		#if sandbox:
		#	maxGaiaNoise = 99.0 # Significance!
		#	minGaiaNoise = 2.0 # Significance!

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		self.run_config = {} if run_config is None else run_config
		self.name = "TNSTalker"
		self.logger.info("Initialized T3 TNSTalker instance %s"%self.name)
		self.logger.info("base_config: %s"%self.base_config)
		self.logger.info("run_config: %s"%self.run_config)
		
		# parse run_config - tns commands
		self.get_tns				= pedantic_get(run_config, 'get_tns')
		self.get_tns_force			= pedantic_get(run_config, 'get_tns_force')
		self.submit_tns				= pedantic_get(run_config, 'submit_tns')
		self.resubmit_tns_nonztf	= pedantic_get(run_config, 'resubmit_tns_nonztf')
		self.resubmit_tns_ztf 		= pedantic_get(run_config, 'resubmit_tns_ztf')
		self.sandbox				= pedantic_get(run_config, 'sandbox')
		self.tns_api_key				= pedantic_get(run_config, 'tns_api_key')
		self.ext_journal			= pedantic_get(run_config, 'ext_journal')
		
		# AT request parameters
		self.base_at_dict			= pedantic_get(run_config, 'base_at_dict')
		self.ztf_tns_at				= pedantic_get(run_config, 'ztf_tns_at')
		self.max_maglim				= pedantic_get(run_config, 'max_maglim')
		self.nphot_submit			= pedantic_get(run_config, 'nphot_submit')
		self.max_gaia_noise			= pedantic_get(run_config, 'max_gaia_noise')
		
#		# TODO: do we want to leave this
#		if self.api_key is None:
#			#raise KeyError("No TNS api_key, cannot run.")
#			self.logger.info("No TNS api_key, using default + sandbox")	
#			self.api_key = "a3f9bcbbe6a26a1ae97a0eeefe465932e68cba83"
#			self.sandbox = True
		
		# parse run_config - selection parameters
		self.max_redshift	= pedantic_get(run_config, 'max_redshift')
		self.min_redshift	= pedantic_get(run_config, 'min_redshift')
		self.min_ndet		= pedantic_get(run_config, 'min_ndet')
		self.min_ndet_postul= pedantic_get(run_config, 'min_ndet_postul')
		self.max_age		= pedantic_get(run_config, 'max_age')
		self.min_age		= pedantic_get(run_config, 'min_age')
		self.min_peak_mag	= pedantic_get(run_config, 'min_peak_mag')
		self.max_peak_mag	= pedantic_get(run_config, 'max_peak_mag')
		self.min_n_filters	= pedantic_get(run_config, 'min_n_filters')
		self.min_gal_lat	= pedantic_get(run_config, 'min_gal_lat')
		self.lc_filters		= pedantic_get(run_config, 'lc_filters')
		self.needed_catalogs= pedantic_get(run_config, 'needed_catalogs')


	def search_journal_tns(self, tran_view):
		"""
			Look through the journal for a TNS name.
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


	def get_tnsname_journal(self, tran_view):
		"""
			look for the TNS name of this transient, either in the 
			journal, or in the TNS directly.
			
			This function will create and return a JournalUpdate instance, together
			with the tns_name and tns_internals lists:
			
			tns_name, tns_internals, jup = self.get_tnsname_journal(tran_view)
		"""
		
		# first try to look for the TNS name in the journal
		tns_name, tns_internals = self.search_journal_tns(tran_view)
		
		# if you can't find it, or if you want to re-scan, go and look in the TNS
		jup = None
		if (tns_name is None and self.get_tns) or self.get_tns_force:
			
			# query the TNS for transient at this position
			ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean", filters=self.lc_filters)
			new_tns_name, new_internal = get_tnsname(
				ra=ra, dec=dec, api_key=self.tns_api_key, logger=self.logger, sandbox=self.sandbox)#self.get_tnsname(tran_view)
			
			# Create new journal entry if we have a tns_name
			if not new_tns_name is None:
				
				# what happen if you have a new name that is different from the old one?
				if tns_name is not None and not tns_name==new_tns_name:
					self.logger.info("Adding new TNS name",extra={"tnsOld":tns_name,"tnsNew":new_tns_name})
				
				# create content of journal entry. Eventually 
				# update the list with the new internal names if any are found
				jcontent = {'t3unit': self.name, 'tnsName': new_tns_name}
				if new_internal is not None:
					tns_internals.append(new_internal)
					jcontent.update({'tnsInternal':new_internal})
				
				# create a journalUpdate and update the tns_name as well. TODO: check with JNo
				jup = JournalUpdate(tran_id=tran_view.tran_id, ext=self.ext_journal, content=jcontent)
				tns_name = new_tns_name
		
		# if you have already your name, then don't add to journal
		return tns_name, tns_internals, jup

	def acctept_tview(self, tran_view):
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
		pps = lc.get_photopoints(filters=self.lc_filters)
		self.logger.debug("%d photop. passed filter %s"%(len(pps), self.lc_filters))
		
		# cut on number of detection
		if len(pps) < self.min_ndet:
			self.logger.debug("not enough detections: got %d, required %d"%
				(len(pps), self.min_ndet))
			return False
		
		# cut on number of detection after last SIGNIFICANT UL
		ulims = lc.get_upperlimits(filters={'attribute':'diffmaglim', 'operator':'>=', 'value':self.max_maglim})
		if len(ulims) > 0:
			last_ulim_jd = sorted(ulims, key=lambda x: x.get_value('jd'))[-1].get_value('jd')
			pps_after_ndet = lc.get_photopoints( 
				filters = self.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_ulim_jd}])
			if len(pps_after_ndet) < self.min_ndet_postul:
				self.logger.debug("not enough consecutive detections after last significant UL.",
					extra={'NDet': len(pps), 'lastUlimJD': last_ulim_jd})
				return False
		
		# cut on number of filters
		used_filters = set([pp.get_value('fid') for pp in pps])
		if len(used_filters) < self.min_n_filters:
			self.logger.debug("requested detections in more than %d bands, got: %d"%
				(self.min_n_filters, len(used_filters)))
			return False
		
		# cut on range of peak magnitude
		mags = [pp.get_value('magpsf') for pp in pps]
		peak_mag = min(mags)
		if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
			self.logger.debug("peak magnitude of %.2f outside of range [%.2f, %.2f]"%
				(peak_mag, self.min_peak_mag, self.max_peak_mag))
			return False
		
		# cut on age
		jds = [pp.get_value('jd') for pp in pps]
		age = Time.now().jd - min(jds)
		if age > self.max_age or age < self.min_age:
			self.logger.debug("age of %.2f days outside of range [%.2f, %.2f]"%
				(age, self.min_age, self.max_age))
			return False
		
		# cut on galactic coordinates
		ra, dec = lc.get_pos(ret="mean", filters=self.lc_filters)
		coordinates = SkyCoord(ra, dec, unit='deg')
		b = coordinates.galactic.b.deg
		if abs(b) < self.min_gal_lat:
			self.logger.debug("transient at b=%.2f too close to galactic plane (cut at %.2f)"%
				(b, self.min_gal_lat))
			return False
		
		# ----------------------------------------------------------------------#
		# 							CUTS ON T2 RECORDS							#
		# ----------------------------------------------------------------------#
		cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)
		
		# check that you have positive match in all of the necessary cataslogs:
		for needed_cat in self.needed_catalogs:
			if not cat_res.get(needed_cat, False):
				self.logger.debug("no T2CATALOGMATCH results", extra={'catalog_matches': cat_res})
				return False
		
		nedz		= cat_res.get('NEDz', False)
		sdss_spec 	= cat_res.get("SDSS_spec", False)
		if ((nedz and not (self.min_redshift <  nedz['z'] < self.max_redshift)) or 
			(sdss_spec and not (self.min_redshift < sdss_spec['z'] > self.max_redshift))):
			self.logger.debug("transient z above limit.", extra={'max_z': self.max_redshift, 'SDSSspec': sdss_spec, 'NEDz': nedz})
			return False
		
		# cut stars in SDSS	#TODO: veryfy!
		sdss_dr10 = cat_res.get('SDSSDR10', False)
		if sdss_dr10 and sdss_dr10['type'] == 6:
			self.logger.debug("transient matched with star in SDSS_DR10.", extra=sdss_dr10)
			return False
		
		# cut matches with variable star catalog
		aavsovsx = cat_res.get('AAVSOVSX', False)
		if aavsovsx and aavsovsx['dist2transient'] < 1:
			self.logger.debug("transient too close to AAVSOVSX sorce", extra=aavsovsx)
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
		if ( milliquas and milliquas['redshift'] > 0) or (sdss_spec and sdss_spec['bptclass'] in [4, 5]):
			self.logger.info("Transient is SDSS BPT or Milliquas AGN.",
				extra={"tranId":tran_view.tran_id, 'milliquas': milliquas, 'SDSS_spec': sdss_spec})
			return {
					"remarks": "Known SDSS and/or MILLIQUAS QSO/AGN. ",
					"at_type": 3
				}
		
		# tag nuclear
		sdss_dr10 = cat_res.get('SDSSDR10', False)
		if sdss_dr10 and sdss_dr10['type'] == 3 and sdss_dr10['dist2transient'] < self.nuclear_dist:
			self.logger.info("Transient close to SDSS photometric galaxy - possibly nuclear",
				extra={"tranId":tran_view.tran_id, 'SDSSDR10': sdss_dr10})
			return {
					"remarks": "Close to core of SDSS DR10 galaxy",
					"at_type": 4
				}
		
		# tag noisy gaia
		gaia_dr2 = cat_res.get('GAIADR2', False)
		nedz	 = cat_res.get('NEDz', False)
		if ( (gaia_dr2 and gaia_dr2['ExcessNoise'] > self.max_gaia_noise and gaia_dr2['dist2transient'] < 1) and 
			 (nedz and not (nedz['z']>0.01 and nedz['dist2transient'] < 1)) and 				#if it's not extragalactic
			 (sdss_dr10 and not (sdss_dr10['type'] == 3 and sdss_dr10['dist2transient'] <3))	# and if it's not a galaxy
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
		ra, dec = lc.get_pos(ret="mean", filters=self.lc_filters)
		
		# Start defining AT dict: name and position
		atdict = {}
		atdict.update(self.base_at_dict)
		atdict["internal_name"] = ztf_name
		atdict["ra"] = {"value": ra, "error" : 1., "units":"arcsec"}
		atdict["dec"] = {"value": dec, "error" : 1., "units":"arcsec"}
		
		# Add information on the latest SIGNIFICANT non detection. TODO: check!
		ulims = lc.get_upperlimits(filters={'attribute':'diffmaglim', 'operator':'>=', 'value':self.max_maglim})
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
		atdict["non_detection"].update(self.ztf_tns_at)# Add the default ZTF values
		
		# now add info on photometric detections: consider only candidates which
		# have some consecutive detection after the last ulim
		pps = lc.get_photopoints(
			filters = self.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_non_obs}])
		
		# Lets create a few photometry points: TODO: should they be the latest or the first?
		atdict["photometry"] = {"photometry_group":{}}
		atdict["discovery_datetime"] = 10**30
		for ipp, pp in enumerate(pps[:self.nphot_submit]):
			photdict = {	#TODO: do we need to round the numerical values?
				"obsdate"		: pp.get_value('jd'),
				"flux"			: pp.get_value('magpsf'),
				"flux_error"	: pp.get_value('sigmapsf'),
				"limiting_flux"	: pp.get_value('diffmaglim'),
				"filter_value"	: TNSFILTERID.get(pp.get_value('fid'))
				}
			if pp.get_value('jd')<atdict["discovery_datetime"]:
				atdict["discovery_datetime"] = pp.get_value('jd')
			photdict.update(self.ztf_tns_at)
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
		transients_to_submit = [tv for tv in transients if self.acctept_tview(tv)]
		self.logger.info("of the %d transients presented to this task, %d passed selection criteria"%
			(len(transients), len(transients_to_submit)))
		
		journal_updates = []	# Will be saved to future journals
		atreports = {}			# Reports to be sent, indexed by the transient view IDs (so that we can check in the replies)
		for tran_view in transients_to_submit:
			
			ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
			self.logger.info("TNS start", extra={"tranId":tran_view.tran_id, 'ztfName': ztf_name})
			
			# find the TNS name, either from the journal of from TNS itself,
			# if necessary, create a new JournalUpdate
			tns_name, tns_internals, jup = self.get_tnsname_journal(tran_view)
			if not jup is None:
				journal_updates.append(jup)
			
			# Chech whether this ID has been submitted (note that we do not check 
			# whether the same candidate was ubmitte as different ZTF name) and
			# depending on what's already on the TNS we can chose to submit or not
			is_ztfsubmitted = ztf_name in tns_internals
			if not ( (is_ztfsubmitted and self.resubmit_tns_ztf) or 
					 (not is_ztfsubmitted and self.resubmit_tns_nonztf) ):
				self.logger.debug("we won't submit candidate.", extra={'is_ztfsub': is_ztfsubmitted})
				continue
			
			# create AT report
			atreport = self.create_atreport(tran_view)
			self.logger.info("Added to report list")
			atreports[tran_view.tran_id] = atreport
		
		# TODO: we save the atreports to send them to the TNS. 
		# This is just part of the tesing and will have to go away
#		atreports = {k: atreports[k] for k in list(atreports.keys())[:2]}
		self.atreports = atreports
		self.logger.info("collected %d AT reports to post"%len(atreports))
		
		# If we do not want to submit anything, or if there's nothing to submit
		if len(atreports) == 0 or (not self.submit_tns):
			self.logger.info("submit_tns config parameter is False or there's nothing to submit", 
				extra={'n_reports': len(atreports), 'submit_tns': self.submit_tns})
			return journal_updates
		
		# Send reports in chunks of size 90 (99 should work)
		atchunks = list(chunks([atr for atr in atreports.values()], 90))
		tnsreplies = sendTNSreports(atchunks, self.tns_api_key, self.logger, sandbox=self.sandbox)
		
		# Now go and check and create journal updates for the cases where SN was added
		for tran_id in atreports.keys():
			ztf_name = ZTFUtils.to_ztf_id(tran_id)
			if not ztf_name in tnsreplies.keys():
				self.logger.info("No TNS add reply",extra={"tranId":tran_id})
				continue
			
			# Create new journal entry 			#TODO: do we want to add to the journal a failed TNS submit?
			jup = JournalUpdate(
					tranId=tran_id,
					ext=self.ext_journal,
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
		
		slack_token			= "xoxb-297846339667-549790069252-FLwKXkra0NL3FNnrvu9XYm4a"
		slack_channel 		= "#ampel_tns_test"
		slack_username		= "Ampel_TNS_test"
		max_slackmsg_size	= 200	# if you have more than this # of reports, send different files
		
		sc = SlackClient(slack_token)
		
		tstamp = datetime.datetime.today().strftime("%Y-%m-%d-%X")
		atlist = list(self.atreports.values())
		last = 0
		for ic, atrep in enumerate(chunks(atlist, max_slackmsg_size)):
			
			# add the atreport to a file
			self.logger.debug("Posting chunk #%d"%ic)
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
					token = slack_token,
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
