#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2PS1ThumbNedSNCosmo.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                12.09.2021
# Last Modified Date:  14.10.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from collections.abc import Mapping, Sequence
from typing import Literal

from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.contrib.hu.t2.T2PanStarrThumbPrint import T2PanStarrThumbPrint
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.PlotProperties import FormatModel, PlotProperties
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.util.collections import ampel_iter
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView


class T2PS1ThumbNedSNCosmo(AbsTiedLightCurveT2Unit):
    """
    This state t2 unit is tied with the state T2 unit T2NedSNCosmo.
    It retrieves panstarrs images at a location encoded in the first datapoint and for each SNCosmo result:
    - creates a new image
    - marks the (first) datapoint location
    - marks the matched location from the NED (encoded in the NedSNComso result under key 'ned')
    """

    #: band: example: ["g", "r", "i", "z", "y"]
    band: str | Sequence[str] = "g"
    cmaps: Sequence[str] = ["cividis"]
    plot_all: bool = False
    z_range: None | tuple[float, float]
    spectroscopic: bool = True
    merge_tags: bool = True
    # Copy selected keys from T2NedCosmo results such as 'fit_results' or 'sncosmo_info'
    copy_keys: list[str] = []
    t2_dependency: Sequence[StateT2Dependency[Literal["T2NedSNCosmo"]]]

    plot_props: PlotProperties = PlotProperties(
        tags=["THUMBPRINT", "PANSTARRS"],
        file_name=FormatModel(
            format_str="%s_%s_%s_ps1_%s_%i_thumb.svg",
            arg_keys=["stock", "obj_name", "band", "cmap", "index_pos"],
        ),
        title=FormatModel(
            format_str="%s - %s\nz=%s (%s band, pos %i) ",
            arg_keys=["stock", "obj_name", "z", "band", "index_pos"],
        ),
        id_mapper="ZTFIdMapper",
    )

    def process(
        self, light_curve: LightCurve, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:
        """
        :param light_curve: see "ampel.view.LightCurve" docstring for more info.
        """

        if not light_curve.photopoints:
            self.logger.info("Lightcurve does not contain any photopoints")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        # That would be a config error
        if not t2_views[0].is_state_type():
            return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

        snc_result = t2_views[0].get_payload()
        if not isinstance(snc_result, Mapping):
            return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

        if "data" not in snc_result:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        if not isinstance(snc_result["data"], list):
            return 1

        ret: dict = {"data": []}
        datapoint = light_curve.photopoints[0]
        pt = T2PanStarrThumbPrint.get_ps1_target(datapoint, self.band)

        for i, d in enumerate(snc_result["data"]):
            for k in ("model_analysis", "fit_results"):
                if k not in d:
                    self.logger.info(
                        f"Skipping incompatible sncosmo with index {i} ('{k}' missing)"
                    )
                    continue

            for k in ("has_premax_data", "has_postmax_data"):
                if not d.get("model_analysis", {}).get(k, False):
                    self.logger.info(
                        f"Skipping sncosmo with index {i}: {k} is False or missing"
                    )
                    continue

            for band in ampel_iter(self.band):
                for cmap in ampel_iter(self.cmaps):
                    d2 = {
                        "plot": [
                            create_plot_record(
                                pt.show(
                                    ellipse=False,
                                    band=band,
                                    cmap=cmap,
                                    show=False,
                                    show_target=False,
                                    show_coord=(
                                        d["catalog"]["ra"],
                                        d["catalog"]["dec"],
                                    ),
                                ),
                                self.plot_props,
                                extra={
                                    "band": band,
                                    "stock": datapoint["stock"][0],  # type: ignore[index]
                                    "cmap": cmap,
                                    "obj_name": d["catalog"]["prefname"].replace(
                                        " ", "_"
                                    ),
                                    "z": d["catalog"]["z"],
                                    "index_pos": i,
                                },
                                logger=self.logger,
                            )
                        ]
                    }

                    for k in self.copy_keys:
                        if k in d:
                            d2[k] = d[k]

                    ret["data"].append(d2)

            if not self.plot_all:
                break

        # Assimilate T2NedSNCosmo tags into the T2PS1ThumbNedSNCosmo doc tags for convenience
        if self.merge_tags:
            return UnitResult(tag=t2_views[0].tag, body=ret if ret["data"] else None)

        return ret if ret["data"] else None
