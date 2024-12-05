#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2PS1ThumbExtCat.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                31.01.2021
# Last Modified Date:  14.09.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from collections.abc import Mapping, Sequence
from typing import Literal

from ampel.abstract.AbsTiedPointT2Unit import AbsTiedPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t2.T2PanStarrThumbPrint import T2PanStarrThumbPrint
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.PlotProperties import FormatModel, PlotProperties
from ampel.model.UnitModel import UnitModel
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.util.collections import ampel_iter
from ampel.view.T2DocView import T2DocView


class T2PS1ThumbExtCat(AbsTiedPointT2Unit):
    """
    Retrieve panstarrs images at datapoint location and for each tied extcat catalog matching result:
    - create a new image
    - mark the datapoint location
    - mark the matched location from the catalog

    A dict structure containing each image as a compressed svg is returned.
    """

    t2_dependency: Sequence[UnitModel[Literal["T2CatalogMatch"]]]

    cmaps: Sequence[str] = ["cividis"]
    band: str | Sequence[str] = "g"
    plot_props: PlotProperties = PlotProperties(
        tags=["THUMBPRINT", "PANSTARRS"],
        file_name=FormatModel(
            format_str="%s_%s_%s_ps1_thumb.svg", arg_keys=["stock", "band", "catalog"]
        ),
        title=FormatModel(
            format_str="%s %s (%s band) ", arg_keys=["stock", "catalog", "band"]
        ),
        id_mapper="ZTFIdMapper",
    )

    def process(
        self, datapoint: DataPoint, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:
        """ """

        # That would be a config error
        if not t2_views[0].is_point_type():
            return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

        cat_results = t2_views[0].get_payload()
        if not isinstance(cat_results, Mapping):
            return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

        if "data" not in cat_results:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        pt = T2PanStarrThumbPrint.get_ps1_target(datapoint, self.band)
        plots = [
            create_plot_record(
                pt.show(
                    ellipse=False,
                    band=band,
                    cmap=cmap,
                    show=False,
                    show_target=False,
                    show_coord=(cat_res["ra"], cat_res["dec"]),
                ),
                self.plot_props,
                extra={
                    "band": band,
                    "stock": datapoint["stock"][0],  # type: ignore[index]
                    "catalog": cat_name,
                    "cmap": cmap,
                },
                logger=self.logger,
            )
            for cat_name, cat_res in cat_results["data"].items()
            if cat_res
            for band in ampel_iter(self.band)
            for cmap in ampel_iter(self.cmaps)
        ]

        return {"data": {"plots": plots}}
