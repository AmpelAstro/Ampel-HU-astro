#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/ingest/ZiDataPointShaper.py
# License:             BSD-3-Clause
# Author:              jno <firstname.lastname@gmail.com>
# Date:                26.04.2013
# Last Modified Date:  28.04.2023
# Last Modified By:    jno <firstname.lastname@gmail.com>

from typing import Any
from bson import encode
from ampel.types import StockId
from ampel.abstract.AbsT0Unit import AbsT0Unit
from ampel.content.DataPoint import DataPoint
from ampel.util.hash import hash_payload
from ampel.ztf.ingest.ZiDataPointShaper import ZiDataPointShaperBase


class ZiGWDataPointShaper(ZiDataPointShaperBase, AbsT0Unit):
	"""
	Shaper used when a (ZTF) alert is ingested based on a match with a GW contour.
	Will, besides the ZTF photopoints also add a datapoint representing the Healpix map it was generated from.
	
	
	Q: How to cretae a uniq id? Hash of healpix_info? Should be unique wrt other healpix maps, but will it be os
	with respect to other datapoints?
	Fix origin through first digit. also propagate e.g. to PPSFilter.
	For now we force id to be neg to make sure it is not selected by PPSFilter. 
	
	"""

	map_name: str         # Used to generate the ID, will also be added to the dp tags
	healpix_info: dict    # Directly propagated as the dp body

	#: Byte width of datapoint ids
	digest_size: int = 8
	
	def process(self, arg: Any, stock: None | StockId = None) -> list[DataPoint]:
		assert stock is not None
		
		dplist = super().process(arg, stock)
		
		# Create new datapoint 
		sorted_body = dict(sorted(self.healpix_info.items()))
		dplist.append(
			{    # type: ignore[typeddict-item]
				'id': -abs(hash_payload(encode(sorted_body), size=-self.digest_size*8)),
				'tag': [self.map_name],
				'stock': stock,
				'body': self.healpix_info
			}
		)		
		
		return dplist




