#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File			  : ampel/contrib/hu/t3/T3MarshalMonitor
# License		   : BSD-3-Clause
# Author			: jnordin@physik.hu-berlin.de
# Date			  : 17.11.2018
# Last Modified Date: 06.02.2019
# Last Modified By  : matteo.giomi@desy.de

import re
from pydantic import BaseModel, BaseConfig
from typing import Dict


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

def is_sdss_star(tran_view):
	"""
	Check whether the there is CATALOGMATCH (T2) info 
	suggesting this an SDSS photometric star
	"""
	# Look for catalog matching
	if tran_view.t2records is not None:
		for j, t2record in enumerate(tran_view.t2records):
			# Look for catalog matching
			if not t2record.get_t2_unit_id() == "CATALOGMATCH":
				continue

			foo = t2record.get_results()
			if foo is None:
				continue
			res = foo[-1]
			
			print (t2record)
			
			# Try to check for spectroscopic class in SDSS
			if 'SDSSDR10' in res.keys():
				sdssphot = res['SDSSDR10'].get("type",None)
				if sdssphot=="STAR":
					return True

	# If we got here we did not get a match
	return False


def is_agn(tran_view):
	"""
	Check whether the there is CATALOGMATCH (T2) info suggesting this is a QSO/AGN
	"""

	# Look for catalog matching
	if tran_view.t2records is not None:
		for j, t2record in enumerate(tran_view.t2records):
			# Look for catalog matching
			if not t2record.get_t2_unit_id() == "CATALOGMATCH":
				continue

			foo = t2record.get_results()
			if foo is None:
				continue
			res = foo[-1]

			# Try to check for spectroscopic class in SDSS
			if 'SDSS_spec' in res.keys():
				sdssclass = res['SDSS_spec'].get("bptclass",-1)
				if sdssclass==4 or sdssclass==5:
					return True

			# Check if there is a million quasar catalog z
			if 'milliquas' in res.keys():
				mz = res['milliquas'].get("redshift",-1)
				if mz>0:
					return True

				
	# If we got here we did not get a match
	return False





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
		api_key				: str	= None		# Bot api key frm TNS
		get_tns				: bool	= False		# Check for TNS names and add to journal
		get_tns_force		: bool	= False		# Check for TNS info even if internal name is known
		submit_tns	 		: bool	= True		# Submit candidates passing criteria 
		resubmit_tns_nonztf	: bool	= True		# Resubmit candidate submitted w/o the same ZTF internal ID 
		resubmit_tns_ztf	: bool	= False		# Resubmit candidates even if they have been added with this name before
		n_photo_submit		: int	= 2			# how many photometric points we add to the TNS report
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
		
		# TODO: dry run method / default API key?
		
		# cut on transients:
		max_redshift	: float	= 1.15	# maximum redshift from T2 CATALOGMATCH catalogs
		min_ndet		: int	= 2		# A candidate need to have at least this many detections
		max_age			: float = 5		# days, If a detection has an age older than this, skip (stars,age).
		min_age			: float = 0		# Min age of detection history
		min_peak_mag	: float	= 19.5	# range of peak magnitudes for submission
		max_peak_mag	: float = 13	#
		min_n_filters	: int	= 1		# Reported detections in at least this many filters
		min_gal_lat		: float = 14	# Minimal galactic latitide
		min_sharpness	: float = -10.15# Remove datapoints with sharpness below this value? TODO!
		cut_sdss_stars	: bool 	= True  # Cut objects with SDSS photometric definition as star
		
		
		# tag
		nuclear_dist: float = 1.	# Tag objects this close to SDSS galaxies as nuclear. Use negative to disable
		aav_dist	: float = 1.	# Required distance to match with aav catalog. TODO: move?
		
		
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
		
		self.name = "TNSTalker"		#used to tag entries in the journal
		
		# parse run_config - tns commands
		self.get_tns				= pedantic_get(run_config, 'get_tns')
		self.get_tns_force			= pedantic_get(run_config, 'get_tns_force')
		self.submit_tns				= pedantic_get(run_config, 'submit_tns')
		self.resubmit_tns_nonztf	= pedantic_get(run_config, 'resubmit_tns_nonztf')
		self.resubmit_tns_ztf 		= pedantic_get(run_config, 'resubmit_tns_ztf')
		self.sandbox				= pedantic_get(run_config, 'sandbox')
		self.api_key				= pedantic_get(run_config, 'api_key')
		self.ext_journal			= pedantic_get(run_config, 'ext_journal')
		
		# AT request parameters
		self.base_at_dict			= pedantic_get(run_config, 'base_at_dict')
		self.ztf_tns_at				= pedantic_get(run_config, 'ztf_tns_at')
		self.max_maglim				= pedantic_get(run_config, 'max_maglim')
		self.nphot_submit			= pedantic_get(run_config, 'nphot_submit')
		
#		# TODO: do we want to leave this
#		if self.api_key is None:
#			#raise KeyError("No TNS api_key, cannot run.")
#			self.logger.info("No TNS api_key, using default + sandbox")	
#			self.api_key = "a3f9bcbbe6a26a1ae97a0eeefe465932e68cba83"
#			self.sandbox = True
		
		# parse run_config - selection parameters
		self.max_redshift	= pedantic_get(run_config, 'max_redshift')
		self.min_ndet		= pedantic_get(run_config, 'min_ndet')
		self.max_age		= pedantic_get(run_config, 'max_age')
		self.min_age		= pedantic_get(run_config, 'min_age')
		self.min_peak_mag	= pedantic_get(run_config, 'min_peak_mag')
		self.max_peak_mag	= pedantic_get(run_config, 'max_peak_mag')
		self.min_n_filters	= pedantic_get(run_config, 'min_n_filters')
		self.min_gal_lat	= pedantic_get(run_config, 'min_gal_lat')
		self.min_sharpness	= pedantic_get(run_config, 'min_sharpness')
		
		# TODO: can we include in the select- run config part of the parameters?
		self.cut_sdss_stars	= pedantic_get(run_config, 'cut_sdss_stars')

		# Using external journal
		self.ext_journal = True

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
			ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean")
			new_tns_name, new_internal = get_tnsname(
				ra=ra, dec=dec, api_key=self.api_key, logger=self.logger, sandbox=self.sandbox)#self.get_tnsname(tran_view)
			
			# Create new journal entry if we have a tns_name
			if not new_tns_name is None:
				
				# what happen if you have a new name that is different from the old one?
				if tns_name is not None and not tns_name==new_tns_name:
					self.logger.info("Adding new TNS name",extra={"tnsOld":tns_name,"tnsNew":new_tns_name})
				
				# create content of journal entry
				jcontent = {'t3unit': self.name, 'tnsName': new_tns_name}
				
				# update the list with the new internal names if any are found
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
		"""
		
		# get the latest light curve
		lc = tran_view.get_latest_lightcurve()
		
		# apply cut on history: consider photophoints which are sharp enough
		pp_filter = {'attribute': 'sharpnr', 'operator': '>=', 'value': self.min_sharpness}
		pps = lc.get_photopoints(filters=pp_filter)
		self.logger.debug("%d photop. passed filter %s"%(len(pps), pp_filter))
		
		# cut on number of detection
		if len(pps) < self.min_ndet:
			self.logger.debug("not enough detections: got %d, required %d"%
				(len(pps), self.min_ndet))
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
		ra, dec = lc.get_pos(ret="mean", filters=pp_filter)
		coordinates = SkyCoord(ra, dec, unit='deg')
		b = coordinates.galactic.b.deg
		if abs(b) < self.min_gal_lat:
			self.logger.debug("transient at b=%.2f too close to galactic plane (cut at %.2f)"%
				(b, self.min_gal_lat))
			return False
		
		# congratulation TransientView, you made it!
		return True
#		self.cut_sdss_stars	= pedantic_get(run_config, 'cut_sdss_stars')

	def get_atreport(self,tran_view):
		"""
			Collect the data needed for the atreport. Return None in case 
			you have to skip this transient for some reason.
		"""
		
		self.logger.info("creating AT report for transient.")
		ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
		lc = tran_view.get_latest_lightcurve()
		ra, dec = lc.get_pos(ret="mean")
		
		# Start defining AT dict: name and position
		atdict = {}
		atdict.update(self.base_at_dict)
		atdict["internal_name"] = ztf_name
		atdict["ra"] = {"value": ra, "error" : 1., "units":"arcsec"}
		atdict["dec"] = {"value": dec, "error" : 1., "units":"arcsec"}
		
#		# now add remarks based on catalog matching
#		atdict["remarks"] = ""
#		if sne[1]["ztf_name"] in list(agn_sne):
#			atdict["remarks"] = atdict["remarks"] + "Known SDSS and/or MILLIQUAS QSO/AGN. "
#			atdict["at_type"] = 3
#			#sys.exit('have an agn') 
#		elif sne[1]["ztf_name"] in list(nuclear_sne):
#			print("close to core of sdss gal, but do mark as nuclear for now")
#	#		atdict["remarks"] = "Close to core of SDSS DR10 galaxy"
#	#		atdict["at_type"] = 4
#			#sys.exit('have a nuclar') 
#		elif sne[1]["ztf_name"] in AAVSOVSX.keys():
#			atdict["remarks"] = atdict["remarks"] + "AAVSOVSX match: %s "%(AAVSOVSX[sne[1]["ztf_name"]])
#		if sne[1]["ztf_name"] in list(gaianoisy_sne):
#			atdict["remarks"] = atdict["remarks"] + "Significant noise in Gaia DR2 - variable star cannot be excluded. " 
#			#sys.exit("Time to check whether the gaia noise remark works.")

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
		pps = lc.get_photopoints(filters = {'attribute': 'jd', 'operator': '>=', 'value': last_non_obs})
		if len(pps) < self.min_ndet:	#TODO: all the cut should be together! 
			self.logger.info("too few detections after last significant upper limit.", 
				extra = {'NDet': len(pps), 'lastUlimJD': last_non_obs})
			return None
		
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
			
			# If we do not want to submit anything we can continue
			if not self.submit_tns:
				self.logger.debug("submit_tns config parameter is False, candidate won't be submitted")
				continue
			
			# Chech whether this ID has been submitted (note that we do not check 
			# whether the same candidate was ubmitte as different ZTF name) and
			# depending on what's already on the TNS we can chose to submit or not
			is_ztfsubmitted = ztf_name in tns_internals
			if not ( (is_ztfsubmitted and self.resubmit_tns_ztf) or 
					 (not is_ztfsubmitted and self.resubmit_tns_nonztf) ):
				 self.logger.debug("we won't submit candidate.", extra={'is_ztfsub': is_ztfsubmitted})

			# Check if transient is a star
#			if True and self.is_sdss_star(tran_view):
#				self.logger.info("SDSS photstar")

			# Check if transient is an agn
#			tran_is_agn = is_agn(tran_view)
#			self.logger.info("",extra={"isAgn":is_agn})

			atreport = self.get_atreport(tran_view)
			if not atreport is None:
				self.logger.info("Added to report list")
				atreports[tran_view.tran_id] = atreport
		
		self.logger.info("collected %d AT reports to post"%len(atreports))
		if len(atreports) == 0:		#TODO: do we want to return them even though we haven't submitted anything?
			return journal_updates
		
		# Send reports in chunks of size 90 (99 should work)
		atchunks = list(chunks([atr for atr in atreports.values()], 90))
#		tnsreplies = sendTNSreports(atchunks, self.api_key, self.logger, sandbox=self.sandbox)
		
		# Now go and check and create journal updates for the cases where SN was added
		for tran_id in atreports.keys():
			ztf_name = ZTFUtils.to_ztf_id(tran_id)
			
			#TODO: do we want to add to the journal a failed TNS submit?
			if not ztf_name in tnsreplies.keys():
				self.logger.info("No TNS add reply",extra={"tranId":tran_id})
				continue
			
			# Create new journal entry
			jup = JournalUpdate(
					tranId=tran_id,
					ext=self.ext_journal,
					content={
						't3unit': self.name,
						'tnsName': tnsreplies[ztf_name][1]["TNSName"],
						'tnsInternal': ztf_name,
						'tnsSubmitresult':tnsreplies[ztf_name][0]
					})
			journal_updates.append(jup)

		return journal_updates


	def done(self):
		"""
		"""
		self.logger.info("done running T3")

