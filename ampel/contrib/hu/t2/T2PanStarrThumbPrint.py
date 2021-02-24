#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2ThumbPrinter.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 31.01.2021
# Last Modified Date: 08.02.2021
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import Dict, Union, Sequence
from pymage.panstarrs import PS1Target # type: ignore[import]

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.util.collections import ampel_iter
from ampel.content.DataPoint import DataPoint
from ampel.plot.SVGUtils import SVGUtils
from ampel.type import T2UnitResult
from ampel.model.PlotProperties import PlotProperties


class T2PanStarrThumbPrint(AbsPointT2Unit):
	"""
	Retrieve panstarrs images at datapoint location and save
	these as compressed svg into the returned dict.

	:param band: example: ["g", "r", "i", "z", "y"]
	"""

	ingest: Dict = {"eligible": {"pps": "first"}}
	cmaps: Sequence[str] = ["gray", "viridis"]
	band: Union[str, Sequence[str]] = "g"
	plot_props: PlotProperties = PlotProperties(
		tags = ["THUMBPRINT", "PANSTARRS"],
		file_name = {
			"format_str": "%s_%s_ps1_thumb.svg",
			"arg_keys": ["stock", "band"]
		},
		title = {
			"format_str": "%s (%s band) ",
			"arg_keys": ["stock", "band"]
		},
		id_mapper = "ZTFIdMapper"
	)

	# No need to download data for thumbs
	# def post_init(self):
	#	# Either this or an env var with the dangerous name DATAPATH need to be set
	#	pymage.io.DATAPATH = self.pymage_data_path


	def run(self, datapoint: DataPoint) -> T2UnitResult:
		"""
		:param light_curve: see "ampel.view.LightCurve" docstring for more info.
		"""

		pt = self.get_ps1_target(datapoint)

		return {
			'plots': [
				SVGUtils.mplfig_to_svg_dict1(
					pt.show(ellipse=False, band=band, show_target=False, cmap=cmap, show=False),
					self.plot_props,
					extra = {"band": band, "stock": datapoint["stock"][0], "cmap": cmap}, # type: ignore
					logger = self.logger
				)
				for band in ampel_iter(self.band)
				for cmap in ampel_iter(self.cmaps)
			]
		}


	def get_ps1_target(self, datapoint: DataPoint) -> PS1Target:

		pt = PS1Target(None)
		pt.set_coordinate(
			datapoint["body"]["ra"],
			datapoint["body"]["dec"]
		)
		pt.download_cutout(filters=self.band)
		return pt
