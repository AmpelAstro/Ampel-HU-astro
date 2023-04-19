#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/alert/ResourceDependentConsumer.py
# License           : BSD-3-Clause
# Author            : jno
# Date              : 28.03.2023
# Last Modified Date: 28.03.2023
# Last Modified By  : jno

import re
import json

from ampel.log import AmpelLogger
from ampel.alert.AlertConsumer import AlertConsumer
from ampel.model.ingest.IngestDirective import IngestDirective
from ampel.ingest.ChainedIngestionHandler import ChainedIngestionHandler
from ampel.mongo.update.DBUpdatesBuffer import DBUpdatesBuffer
from ampel.model.ingest.CompilerOptions import CompilerOptions

class ResourceDependentConsumer(AlertConsumer):
	"""
	
	Extension of standard AlertConsumer where the directives provided to the ingestion handler
	can be updated based on (dynamic) resources available to the EventHandler.
	
	Use case is when a unit configuration parameter can not be determined until the job is run. 
	Since updated values are present during configuration building, these will still be part of the
	provenance/hash sequence.
	
	:param directives_map:
	
	For each key of directives_map with a value which exists in the resources,  
	which corresponds to a resource entry, any occurance of the value
	in the directies will be replaced with the resource value.
	

	"""
	
	directives_map: dict[str,str]

	# Overload
	def get_ingestion_handler(self,
		run_id: int,
		updates_buffer: DBUpdatesBuffer,
		logger: AmpelLogger
	) -> ChainedIngestionHandler:


		# Replacement pattern for provided resources
		pattern = re.compile("|".join( 
				[
					k  
					for k, v in self.directives_map.items()
					if v in self.alert_supplier.resources
				]
			)) 
		print('pattern', pattern)
		print('directives pre', self.directives)
		directives = [
			IngestDirective(
				**json.loads( 
					pattern.sub(
						lambda m: self.alert_supplier.resources[self.directives_map[m.group(0)]], 
						json.dumps(el.dict()) 
					)
				)
			)
			for el in self.directives
		]
		print('directives post', directives)

		return ChainedIngestionHandler(
			self.context, self.shaper, directives, updates_buffer,
			run_id, tier = 0, logger = logger, database = self.database,
			trace_id = {'alertconsumer': self._trace_id},
			compiler_opts = self.compiler_opts or CompilerOptions()
		)

