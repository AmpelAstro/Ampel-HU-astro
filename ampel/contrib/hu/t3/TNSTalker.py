#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/T3MarshalMonitor
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 17.11.2018
# Last Modified Date: 17.11.2018
# Last Modified By  : jnordin@physik.hu-berlin.de

import time


from ampel.base.abstract.AbsT3Unit import AbsT3Unit
from ampel.pipeline.logging.AmpelLogger import AmpelLogger
from typing import Union, List, Dict, Any
from ampel.base.dataclass.JournalUpdate import JournalUpdate

# Changed location
#from ampel.pipeline.common.ZTFUtils import ZTFUtils
from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils

def tstamp(): return time.time()

from ampel.contrib.hu.t3.ampel_tns import sendTNSreports, tnsName, tnsInternal

# Create a function called "chunks" with two arguments, l and n:
def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


class TNSTalker(AbsT3Unit):
	"""
		Get TNS name if existing, and submit selected candidates
	"""
	
	version = 0.1

	def __init__(self, logger, base_config=None, run_config=None, global_info=None):
		"""
		
		run_config: `dict`
					configuration parameter for this job. An example for this dictionary
					is provided below:
					
					
						run_config = {
								"api_key": "xxx",               # Bot api key frm TNS
								"get_tns": True,                # Check for TNS names and add to journal
								"get_tns_force": False,         # Check for TNS info even if internal name is known
								"submit_tns": True,             # Submit candidates passing criteria 
								"resubmit_tns_nonztf": True,    # Resubmit candidate submitted w/o the same ZTF internal ID 
								"resubmit_tns_ztf": False,      # Resubmit candidates even if they have been added with this name before
								"sandbox": True                 # Submit to TNS sandbox only
							}

		"""
		
		self.logger = logger if logger is not None else logging.getLogger()
		self.base_config = {} if base_config is None else base_config
		self.run_config = {} if run_config is None else run_config
		
		self.name = "TNSTalker"		#used to tag entries in the journal
		
		# parse run_config - tns commands
		self.get_tns = self.run_config.get('get_tns', True)
		self.get_tns_force = self.run_config.get('get_tns_force', False)
		self.submit_tns = self.run_config.get('submit_tns', False)
		self.resubmit_tns_nonztf = self.run_config.get('resubmit_tns_nonztf', False)
		self.resubmit_tns_ztf = self.run_config.get('resubmit_tns_ztf', False)
		self.sandbox = self.run_config.get('sandbox', True)

		self.api_key = self.run_config.get('api_key')
		if self.api_key is None:
			#raise KeyError("No TNS api_key, cannot run.")
			self.logger.info("No TNS api_key, using default + sandbox")	
			self.api_key = "a3f9bcbbe6a26a1ae97a0eeefe465932e68cba83"
			self.sandbox = True		

		# Selecting parameters
		self.minmagpdf = self.run_config.get('minmagpdf', 18.) 
		self.mindet = self.run_config.get('mindet', 4.)
		self.maxage = self.run_config.get('maxage', 30.)
		self.check_minref = self.run_config.get('check_minref', True) 
		self.cut_agn = self.run_config.get('cut_agn', True)
		self.cut_sdsstar = self.run_config.get('cut_sdsstar', True)

	



		
	def _find_latest_journal_entry(self, tran_view, t3unit=self.name, required_keys=[]):
		"""
			return the latest journal entry, possibly belonging to a specific t3 unit and containing all in list of keys
			replace with tran_view method
		"""
		
		# TODO: find decent schema for naming of journal entry
		journal_entry, journal_dt = None, -10**10
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

	def _find_all_journal_entries(self, tran_view, t3unit=self.name, required_keys=[]):
		"""
			return the all journal entry, possibly belonging to a specific t3 unit and containing all in list of keys
			replace with tran_view method
		"""
		
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



	def get_tnsname(self, tran_view):
		"""
			look for names registered at tns
		"""
		tns_name = None

		# Should use the direct coord method, but ok
		phot = tran_view.get_photopoints()
		ras = [pp.get_value("ra") for pp in phot]
		decs = [pp.get_value("dec") for pp in phot]

		# Look for TNS name
		tnsnames, runstatus = tnsName( np.mean(ras), np.mean(decs), self.api_key, sandbox=self.sandbox )
		if re.match('Error', runstatus):
			self.logger.info("TNS get error", extra={"tns_request":runstatus} )
			return None, []
		if len(tnsnames)>1:
			self.logger.info("Multipe TNS names, choosing first", extra={"tns_names":tnsnames} )
			tns_name = tnsnames[0]
		elif len(tnsnames)==1:
			tns_name = tnsnames[0]
		elif len(tnsnames)==0:
			# No TNS name, then no need to look for internals
			return tns_name, tns_internal
		self.logger.info("TNS search",extra={"tns_name":tns_name})


		# Look for internal name (note that have to skip the prefix)
		internal_name, runstatus = tnsInternal( tns_name[2:], self.api_key, sandbox=self.sandbox )
		if internal_name is not None:
			self.logger.info("",extra={"tns_internal_name":internal_name})
			if re.search('ZTF',internal_name):		
				if not internal_name == sne[1]["ztf_name"]:
					self.logger.info("TNS registered under other ZTF name %s", extra={"tns_other_internal":internal_name} )

		return tns_name, internal_name


	def search_journal_tns(self, tran_view):
		"""
			Look through the journal for a TNS name
		"""

		# Find the latest tns name (skipping previous)
		jentry = self._find_latest_journal_entry(tran_view, t3unit=self.name, required_keys=['tns_name'])
		tns_name = jentry['tns_name']

		# Find internal names
		jentries = self._find_all_journal_entries(tran_view, t3unit=self.name, required_keys=['tns_internal'])
		tns_internals = [j['tns_internal'] for j in jentries]

		self.logger.info('',extra={'tran_id':tran_view.tran_id,'tns_name':tns_name,'tns_internals':tns_internals})


		return tns_name, tns_internals
			





	def check_submit_criteria(self,tran_view):
		"""
		Check whether the additional criteria for submition are fulfilled.
		These include:
			minmagpdf [float] (we allow this to be different from channel acceptance)
			mindet [int]
			maxage [float]    (do not submit if detections span a range longer than this)
			check_minref [bool] 
			cut_agn [bool]   (if T2 info available)
			cut_sdsstar [bool]
		"""

		# Look for catalog matching
		if self.cut_agn or self.cut_sdsstar:
			if transient.t2records is not None:
				for j, t2record in enumerate(transient.t2records):
					if not t2record.results:
						continue
					res = (t2record.results[-1])
					if not "output" in res:
						continue
					# Check for AGN
					if self.cut_agn:
						if t2record.get_t2_unit_id()=='X':	
							pass
		

		return False
			



	def add(self, transients):
		"""
			Loop through transients and check for TNS names and/or candidates to submit
		"""

		self.logger.inf("Turned off")
		return []

		journal_updates = []  # Will be saved to future journals
		atreports = []        # Reports to be sent

		if transients is not None:
			for tran_view in transients:

				tns_name, tns_internals = self.search_journal_tns(tran_view.tran_id)

				# Search TNS for a name
				if (tns_name is not None and self.get_tns) or self.get_tns_force:
					new_tns_name, new_internal = self.get_tnsname(tran_view)
					if tns_name is not None and not tns_name==new_tns_name:
						self.logger.info("Adding new TNS name",extra={"tns_old":tns_name,"tns_new":new_tns_name})
					tns_internals.append(new_internal)
					# Create new journal entry
					journal_updates.append(
						# Constructing journal update. The module ID and timestamp are added automatically
						# ext saves it to the resiliant db, so persistent across DB resets
						JournalUpdate(
							tranId=tran_view.tran_id,
							ext=True,
							content={
								't3unit': self.name,
								'tns_name': tns_name,
								'tns_internal': new_internal
							}
						)
					)
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
	
				# Does the transient fulfill submit requirements
				if not self.check_submit_criteria(tran_view):
					continue
					atreport = self.get_atreport(tran_view)
				atreports.append(atreport)
		
			if len(atreports)>0:
				# Send reports in chunks of size 90 (99 should work)
				atchunks = list(chunks(atreports,90))
				tnsreplies = sendTNSreports(atchunks, self.api_key, self.logger, sandbox=self.sandbox)
				
				# Create journal updates for the cases where SN was added
				for tran_view in transients:
					ztf_name = ZTFUtils.to_ztf_id(tran_view.tran_id)
					if not ztf_name in tnsreplies.keys():
						self.logger.info("No TNS add reply",extra={"tran_id":tran_view.tran_id})
						continue
					# Create new journal entry
					journal_updates.append(
						# Constructing journal update. The module ID and timestamp are added automatically
						# ext saves it to the resiliant db, so persistent across DB resets
						JournalUpdate(
							tranId=tran_view.tran_id,
							ext=True,
							content={
								't3unit': self.name,
								'tns_name': tnsreplies[ztf_name][1]["TNSName"],
								'tns_internal': ztf_name,
								'tns_submitresult':tnsreplies[ztf_name][0]
							}
						)
					)

		return journal_updates


						
	def done(self):
		"""
		"""
		self.logger.info("done running T3")


