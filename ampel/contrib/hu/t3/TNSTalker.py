#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/T3MarshalMonitor
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 17.11.2018
# Last Modified Date: 06.02.2019
# Last Modified By  : matteo.giomi@desy.de

import re
from pydantic import BaseModel, BaseConfig

from astropy.time import Time
from astropy.coordinates import SkyCoord

from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.base.dataclass.JournalUpdate import JournalUpdate
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils

from ampel.contrib.hu.t3.ampel_tns import sendTNSreports, get_tnsname # tnsName, tnsInternal, 

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
		
		api_key				: str	= None		# Bot api key frm TNS
		get_tns				: bool	= False		# Check for TNS names and add to journal
		get_tns_force		: bool	= False		# Check for TNS info even if internal name is known
		submit_tns	 		: bool	= True		# Submit candidates passing criteria 
		resubmit_tns_nonztf	: bool	= True		# Resubmit candidate submitted w/o the same ZTF internal ID 
		resubmit_tns_ztf	: bool	= False		# Resubmit candidates even if they have been added with this name before
		n_photo_submit		: int	= 2			# how many photometric points we add to the TNS report
		sandbox				: bool	= True		# Submit to TNS sandbox only
		ext_journal			: bool	= True		# weather journal will go to separate collection.
		
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


	def _find_latest_journal_entry(self, tran_view, t3unit=None, required_keys=[]):
		"""
			return the latest journal entry, possibly belonging to a specific 
			t3 unit and containing all in list of keys.
			
			replace with tran_view method
		"""
		if t3unit is None:
			t3unit = self.name
		# TODO: find decent schema for naming of journal entry
		journal_entry, journal_dt = {}, -10**10
		for je in tran_view.journal:

			# Only look at later journals - if key dt is not guaranteed to exist we have to add timestamps for this
			if je['dt']<journal_dt :
				continue
			# Match saved t3unit name if this is set
			if t3unit is not None and not je.get('t3unit', 'fuffa') == t3unit:
				continue

			# Make sure all the required keys are there
			key_lacking = False
			for reqkey in required_keys:
				if not reqkey in je.keys():
					key_lacking = True
			if key_lacking : continue
			
			# Made it - this is a correct entry and the current last one
			journal_entry = je
			journal_dt = je['dt']

		return journal_entry

	def _find_all_journal_entries(self, tran_view, t3unit=None, required_keys=[]):
		"""
			return the all journal entry, possibly belonging to a specific t3 unit
			and containing all in list of keys.
			replace with tran_view method
		"""
		if t3unit is None:
			t3unit = self.name
		# TODO: find decent schema for naming of journal entry
		journal_entries = {}
		for je in tran_view.journal:

			# Match saved t3unit name if this is set
			if t3unit is not None and not je.get('t3unit', 'fuffa') == t3unit:
				continue

			# Make sure all the required keys are there
			key_lacking = False
			for reqkey in required_keys:
				if not reqkey in je.keys():
					key_lacking = True
			if key_lacking : continue
			
			# Made it - this is a correct entry and the current last one
			journal_entries[je['dt']] = je
		return journal_entries


	def search_journal_tns(self, tran_view):
		"""
			Look through the journal for a TNS name
		"""

		# Find the latest tns name (skipping previous)
#		jentry = self._find_latest_journal_entry(tran_view, t3unit=self.name, required_keys=['tnsName'])
#		tns_name = jentry.get("tnsName",None)
		jentry = tran_view.get_journal_entries(
			filterFunc=lambda x: x.get('t3unit')==self.name and 'tnsInternal' in x.keys(),
			latest=True)
		tns_name = None if jentry is None else jentry.get("tnsName", None)

		# Find internal names
#		jentries = self._find_all_journal_entries(tran_view, t3unit=self.name, required_keys=['tnsInternal'])
#		if jentries is not None:
#			tns_internals = [j.get('tnsInternal',None) for dt, j in jentries.items()]
#		else:
#			tns_internals = []
		jentries = tran_view.get_journal_entries(
			filterFunc=lambda x: x.get('t3unit')==self.name and 'tnsInternal' in x.keys())
		tns_internals = [] if jentries is None else [j.get('tnsInternal',None) for j in jentries]
		self.logger.info('Journal search',extra={'tranId':tran_view.tran_id,'tnsName':tns_name,'tnsInternals':tns_internals})
		return tns_name, tns_internals


	def get_tnsname_journal(self, tran_view):
		"""
			look for the TNS name of this transient, either in the 
			journal, or in the TNS directly.
		"""
		
		# first try to look for the TNS name in the journal
		tns_name, tns_internals = self.search_journal_tns(tran_view)
		
		# if you can't find it, or if you want to re-scan, go and look in the TNS
		if (tns_name is None and self.get_tns) or self.get_tns_force:
			ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean")
			new_tns_name, new_internal = get_tnsname(
				ra=ra, dec=dec, api_key=self.api_key, logger=self.logger, sandbox=self.sandbox)#self.get_tnsname(tran_view)
			
			# update / create the journal entries
			if new_internal is not None:
				tns_internals.append(new_internal)

			# Create new journal entry if we have a tns_name
			if new_tns_name is not None:
				if tns_name is not None and not tns_name==new_tns_name:
					self.logger.info("Adding new TNS name",extra={"tnsOld":tns_name,"tnsNew":new_tns_name})
					pass

				jcontent = {'t3unit': self.name, 'tnsName': new_tns_name}
#						jcontent = OrderedDict([ ('t3unit',self.name), ('tns_name',new_tns_name) ] )
				if new_internal is not None:
#							jcontent['tns_internal'] = new_internal
					jcontent.update( {'tnsInternal':new_internal} )
				journal_updates.append(
					# Constructing journal update. The module ID and timestamp are added automatically
					# ext saves it to the resiliant db, so persistent across DB resets
					JournalUpdate(
						tran_id= getattr(tran_view, "tran_id"),
						ext=self.ext_journal,
						content= jcontent
					)
#							{"tran_id":getattr(tran_view, "tran_id"),
#								"ext":True,
#								"content": jcontent }
				)
#						print(journal_updates[-1])




	def is_sdss_star(self,tran_view):
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


	def is_agn(self,tran_view):
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
			Collect the data needed for the atreport
		"""
		
		if self.tag_agn:
			pass
			
		
		
		return {}


	def add(self, transients):
		"""
			Loop through transients and check for TNS names and/or candidates to submit
		"""
		
		
		
		
#		# check the find jentries
#		for tv in transients:
#			here = self._find_all_journal_entries(tv, t3unit=self.name, required_keys=['tnsInternal'])
#			there = tv.get_journal_entries(
#				filterFunc=lambda x: x.get('t3unit') == self.name and 'tnsInternal' in x.keys())
#			
#			my_keys = ['marshalProgram','marshalSave']
#			my_keys = []
#			
#			here = self._find_all_journal_entries(tv, t3unit="MarshalPublisher")#, required_keys=my_keys)
#			there = tv.get_journal_entries(
#				filterFunc=lambda x: x.get('t3unit')=="MarshalPublisher", latest=True)# and set(my_keys).issubset(x.keys()))
#			
#			
#			if here != there:
#				print (here)
#				print (there)
#				input()
#		
#		input()
		
		if transients is None:
			self.logger.info("no transients for this task execution")
			return []
		
		# select the transients
		transients_to_submit = [tv for tv in transients if self.acctept_tview(tv)]
		self.logger.info("of the %d transients presented to this task, %d passed selection criteria"%
			(len(transients), len(transients_to_submit)))
		
#		# Still turned off
#		return []


		journal_updates = []  # Will be saved to future journals
		atreports = []        # Reports to be sent

		for tran_view in transients_to_submit:
			self.logger.info("TNS start", extra={"tranId":tran_view.tran_id})

			# ----- get TNS name, either from the journal or from TNS itself -----#
			tns_name, tns_internals = self.search_journal_tns(tran_view)
			if (tns_name is None and self.get_tns) or self.get_tns_force:
				ra, dec = tran_view.get_latest_lightcurve().get_pos(ret="mean")
				new_tns_name, new_internal = get_tnsname(
					ra=ra, dec=dec, api_key=self.api_key, logger=self.logger, sandbox=self.sandbox)#self.get_tnsname(tran_view)
				
				# update / create the journal entries
				if new_internal is not None:
					tns_internals.append(new_internal)

				# Create new journal entry if we have a tns_name
				if new_tns_name is not None:
					if tns_name is not None and not tns_name==new_tns_name:
						self.logger.info("Adding new TNS name",extra={"tnsOld":tns_name,"tnsNew":new_tns_name})
						pass

					jcontent = {'t3unit': self.name, 'tnsName': new_tns_name}
#						jcontent = OrderedDict([ ('t3unit',self.name), ('tns_name',new_tns_name) ] )
					if new_internal is not None:
#							jcontent['tns_internal'] = new_internal
						jcontent.update( {'tnsInternal':new_internal} )
					journal_updates.append(
						# Constructing journal update. The module ID and timestamp are added automatically
						# ext saves it to the resiliant db, so persistent across DB resets
						JournalUpdate(
							tran_id= getattr(tran_view, "tran_id"),
							ext=self.ext_journal,
							content= jcontent
						)
#							{"tran_id":getattr(tran_view, "tran_id"),
#								"ext":True,
#								"content": jcontent }
					)
#						print(journal_updates[-1])


			# If we do not want to submit anything we can continue
			if not self.submit_tns:
				continue
				
			# Chech whether this ID has been submitted (note that we do not check whether the same candidate was ubmitte as different ZTF name)
			is_ztfsubmitted = False
			ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
			for iname in tns_internals:
				if iname==ztf_name:
					is_ztfsubmitted = True

			# Should we submit
			if not is_ztfsubmitted and self.resubmit_tns_nonztf:
				# Ok!
				pass
			elif is_ztfsubmitted and self.resubmit_tns_ztf:
				# Ok!
				pass
			else:
				# Not ok
				continue

			# Check if transient is a star
			if self.is_sdss_star(tran_view):
				self.logger.info("SDSS photstar")
				if self.cut_sdsstar:
					continue

			# Check if transient is an agn
			is_agn = self.is_agn(tran_view)
			self.logger.info("",extra={"isAgn":is_agn})


#			# Does the transient fulfill submit requirements
#			if not self.check_submit_criteria(tran_view):
#				self.logger.info("Not fulfill submit criteria")
#				continue

			
			atreport = self.get_atreport(tran_view)
			self.logger.info("Added to report list")
			#atreports.append(atreport)
		
			if len(atreports)>0:
				# Send reports in chunks of size 90 (99 should work)
				atchunks = list(chunks(atreports,90))
#				tnsreplies = sendTNSreports(atchunks, self.api_key, self.logger, sandbox=self.sandbox)
				
				# Create journal updates for the cases where SN was added
				for tran_view in transients:
					ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
					if not ztf_name in tnsreplies.keys():
						self.logger.info("No TNS add reply",extra={"tranId":tran_view.tran_id})
						continue
					# Create new journal entry
					journal_updates.append(
						# Constructing journal update. The module ID and timestamp are added automatically
						# ext saves it to the resiliant db, so persistent across DB resets
						JournalUpdate(
							tranId=tran_view.tran_id,
							ext=self.ext_journal,
							content={
								't3unit': self.name,
								'tnsName': tnsreplies[ztf_name][1]["TNSName"],
								'tnsInternal': ztf_name,
								'tnsSubmitresult':tnsreplies[ztf_name][0]
							}
						)
					)


		return journal_updates


	def done(self):
		"""
		"""
		self.logger.info("done running T3")


