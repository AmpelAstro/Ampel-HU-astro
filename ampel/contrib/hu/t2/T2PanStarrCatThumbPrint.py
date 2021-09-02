#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2PanStarrCatThumbPrint.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 31.01.2021
# Last Modified Date: 24.02.2021
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import List, Sequence, Dict, Any
from ampel.contrib.hu.t2.T2PanStarrThumbPrint import T2PanStarrThumbPrint
from ampel.abstract.AbsTiedPointT2Unit import AbsTiedPointT2Unit
from ampel.util.collections import ampel_iter
from ampel.content.DataPoint import DataPoint
from ampel.plot.utils import mplfig_to_svg_dict1
from ampel.enum.T2RunState import T2RunState
from ampel.type import T2UnitResult
from ampel.model.PlotProperties import PlotProperties
from ampel.view.T2DocView import T2DocView


class T2PanStarrCatThumbPrint(AbsTiedPointT2Unit, T2PanStarrThumbPrint): # type: ignore[misc]
	"""
	Retrieve panstarrs images at datapoint location and save
	these as compressed svg into the returned dict.

	:param band: example: ["g", "r", "i", "z", "y"]
	"""

	ingest: Dict[str, Any] = {'eligible': {'pps': 'first'}} # type: ignore
	cmaps: Sequence[str] = ["gray", "viridis"]
	plot_props: PlotProperties = PlotProperties(
		tags = ["THUMBPRINT", "PANSTARRS"],
		file_name = {
			"format_str": "%s_%s_ps1_thumb.svg",
			"arg_keys": ["stock", "catalog"]
		},
		title = {
			"format_str": "%s %s (%s band) ",
			"arg_keys": ["stock", "catalog", "band"]
		},
		id_mapper = "ZTFIdMapper"
	)

	# mandatory
	@classmethod
	def get_tied_unit_names(cls) -> List[str]:
		return ["T2CatalogMatch"]


	def run(self, datapoint: DataPoint, t2_views: Sequence[T2DocView]) -> T2UnitResult: # type: ignore[override]
		"""
		:param light_curve: see "ampel.view.LightCurve" docstring for more info.
		"""
		# That would be a config error
		if not t2_views[0].is_point_type():
			return T2RunState.UNEXPECTED_DEPENDENCY

		payload = t2_views[0].get_payload()
		if payload is None:
			return T2RunState.UNEXPECTED_DEPENDENCY
		cat_results = payload[-1] if isinstance(payload, list) else payload

		pt = self.get_ps1_target(datapoint)

		return {
			'plots': [
				mplfig_to_svg_dict1(
					pt.show(
						ellipse=False, band=band, cmap=cmap, show=False,
						show_target = False, show_coord = (cat_res['ra'], cat_res['dec'])
					),
					self.plot_props,
					extra = {
						"band": band,
						"stock": datapoint["stock"][0], # type: ignore[index]
						"catalog": cat_name,
						"cmap": cmap
					},
					logger = self.logger
				)
				for cat_name, cat_res in cat_results.items() if cat_res
				for band in ampel_iter(self.band)
				for cmap in ampel_iter(self.cmaps)
			]
		}
