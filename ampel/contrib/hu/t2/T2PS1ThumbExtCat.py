#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/contrib/hu/t2/T2PS1ThumbExtCat.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 31.01.2021
# Last Modified Date: 14.09.2021
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import Sequence, Union, Literal
from ampel.types import UBson
from ampel.contrib.hu.t2.T2PanStarrThumbPrint import T2PanStarrThumbPrint
from ampel.abstract.AbsTiedPointT2Unit import AbsTiedPointT2Unit
from ampel.util.collections import ampel_iter
from ampel.content.DataPoint import DataPoint
from ampel.plot.utils import mplfig_to_svg_dict1
from ampel.struct.UnitResult import UnitResult
from ampel.model.PlotProperties import PlotProperties
from ampel.model.UnitModel import UnitModel
from ampel.view.T2DocView import T2DocView
from ampel.enum.DocumentCode import DocumentCode


class T2PS1ThumbExtCat(AbsTiedPointT2Unit):
	"""
	Retrieve panstarrs images at datapoint location and for each tied extcat catalog matching result:
	- create a new image
	- mark the datapoint location
	- mark the matched location from the catalog

	A dict structure containing each image as a compressed svg is returned.
	"""

	t2_dependency: Sequence[UnitModel[Literal['T2CatalogMatch']]]

	cmaps: Sequence[str] = ["cividis"]
	band: Union[str, Sequence[str]] = "g"
	plot_props: PlotProperties = PlotProperties(
		tags = ["THUMBPRINT", "PANSTARRS"],
		file_name = {
			"format_str": "%s_%s_%s_ps1_thumb.svg",
			"arg_keys": ["stock", "band", "catalog"]
		},
		title = {
			"format_str": "%s %s (%s band) ",
			"arg_keys": ["stock", "catalog", "band"]
		},
		id_mapper = "ZTFIdMapper"
	)


	def process(self, datapoint: DataPoint, t2_views: Sequence[T2DocView]) -> Union[UBson, UnitResult]:
		""" """

		# That would be a config error
		if not t2_views[0].is_point_type():
			return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

		cat_results = t2_views[0].get_payload()
		if not isinstance(cat_results, dict):
			return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

		if 'data' not in cat_results:
			return UnitResult(code=DocumentCode.T2_MISSING_INFO)

		pt = T2PanStarrThumbPrint.get_ps1_target(datapoint, self.band)
		plots = []

		for cat_name, cat_res in cat_results['data'].items():

			if not cat_res:
				continue

			for band in ampel_iter(self.band):
				for cmap in ampel_iter(self.cmaps):
					plots.append(
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
					)

		return {'data': {'plots': plots}}
