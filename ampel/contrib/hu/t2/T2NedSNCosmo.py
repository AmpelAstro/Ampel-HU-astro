#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-HU-astro/ampel/contrib/hu/t2/T2NedSNCosmo.py
# License           : BSD-3-Clause
# Author            : vb <vbrinnel@physik.hu-berlin.de>
# Date              : 10.03.2021
# Last Modified Date: 14.09.2021
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import List, Any, Optional, Tuple, Union, Sequence, Literal
from ampel.types import UBson, Tag
from ampel.contrib.hu.t2.T2SNCosmo import T2SNCosmo
from ampel.contrib.hu.util.ned import check_ned_res
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.model.operator.AnyOf import AnyOf
from ampel.model.operator.AllOf import AllOf
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.view.T2DocView import T2DocView
from ampel.view.LightCurve import LightCurve


class T2NedSNCosmo(AbsTiedLightCurveT2Unit, T2SNCosmo):
	"""
	Fits lightcurves using SNCOSMO (using SALT2 defaultwise)
	with redshift constrained by catalog matching results.

	config example:
	{
		"model": "salt2",
		"jd_reject_sigma": 3,
		"plot_props": {
			"tags": ["SALT", "SNCOSMO"],
			"file_name": {
				"format_str": "%s_%s_%s_fit.svg",
				"arg_keys": ["stock", "model", "catalog"]
			},
			"title": {
				"format_str": "%s %s lightcurve fit",
				"arg_keys": ["stock", "catalog"]
			},
			"fig_text": {
				"format_str": "%s %s %s\nchisq %.2f\nndof %s\nfit_ok %s",
				"arg_keys": ["stock", "model", "catalog", "chisq", "ndof", "fit_ok"]
			},
			"width": 10,
			"height": 6,
			"id_mapper": "ZTFIdMapper",
			"disk_save": "/path/to/plots/dir/"
		},
		"t2_dependency": [
			{
				"unit": "T2NEDTap",
				"link_override": {"pps": "first"}
			}
		]
	}
	"""

	t2_dependency: Sequence[StateT2Dependency[Literal["T2NedTap"]]]
	z_range: Optional[Tuple[float, float]]
	fit_all: bool = False # otherwise, fit only the first matching catalog result
	spectroscopic: bool = True
	merge_tags: bool = True
	require_tags: Union[None, AnyOf[Tag], AllOf[Tag]] = None


	# mandatory
	def process(self, light_curve: LightCurve, t2_views: List[T2DocView]) -> Union[UBson, UnitResult]: # type: ignore[override]
		"""
		:param light_curve: see "ampel.view.LightCurve" docstring for more info.
		"""

		if not light_curve.photopoints:
			self.logger.error("Lightcurve has no photopoint")
			return UnitResult(code=DocumentCode.ERROR)

		if not t2_views: # Should not happen actually, T2Worker should catch that case
			self.logger.error("Missing tied t2 views")
			return UnitResult(code=DocumentCode.T2_MISSING_INFO)

		t2_view = t2_views[0]

		# That would be a config error
		if not t2_view.is_point_type():
			return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

		# Unsure this can happen as T2Processor might not allow code execution to go that far
		if (ned_tap_res := t2_view.get_payload()) is None:
			self.logger.error("Tied catalog match t2 record has no content")
			return UnitResult(code=DocumentCode.T2_MISSING_INFO)

		if not isinstance(ned_tap_res, dict) or 'data' not in ned_tap_res:
			return UnitResult(
				code = DocumentCode.T2_MISSING_INFO,
				body = {"msg": "Invalid/incompatible tied T2NedTap result"}
			)

		if self.require_tags:
			if isinstance(self.require_tags, AllOf):
				for t in self.require_tags.all_of:
					if t not in t2_view.tag:
						return UnitResult(code=30)
			else:
				if not [t for t in self.require_tags.any_of if t not in t2_view.tag]:
					return UnitResult(code=30)

		ret: List[Any] = []
		# Assimilate T2NedTag tags into the T2SNCosmo doc tags for convenience
		tags = list(t2_view.tag) if self.merge_tags else []

		for i, cat_res in enumerate(ned_tap_res['data']):

			if check_ned_res(cat_res, self.logger, self.spectroscopic, self.z_range):
				self.logger.info(f"Skipping cat result with index {i}")
				continue

			self.z_discrete = [cat_res['z']]

			try:
				r = T2SNCosmo.process(self, light_curve)
			except Exception as e:
				from traceback import format_exc
				if "No bands in data overlap the model" in format_exc():
					ret.append(
						{
							'name': cat_res['prefname'],
							'z': cat_res['z'],
							'error': "No bands in data overlap the model"
						}
					)
					continue
				raise e # populate "trouble" collection

			if isinstance(r, int):
				ret.append(r)

			elif isinstance(r, (dict, UnitResult)):

				if isinstance(r, UnitResult):
					d: dict = r.body # type: ignore
					if isinstance(r.tag, (int, str)):
						tags.append(r.tag)
					elif isinstance(r.tag, list):
						tags += r.tag
					if r.code:
						d['code'] = r.code
				else:
					d = r

				if "flux_dict" in d:
					del d["flux_dict"] # remove superfluous info

				# Append ned cat result to the sncosmo result for convenience
				d['catalog'] = cat_res
				ret.append(d)

			else:
				self.logger.error("Unsupported T2SNCosmo result, please update T2NedSNCosmo")
				return UnitResult(code=DocumentCode.T2_OUTDATED_CODE)

			if not self.fit_all:
				break

		return UnitResult(tag=tags, body={'data': ret})
