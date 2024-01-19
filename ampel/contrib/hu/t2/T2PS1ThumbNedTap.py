#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2PS1ThumbNedTap.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                10.03.2021
# Last Modified Date:  14.09.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from collections.abc import Sequence
from typing import Literal

from ampel.abstract.AbsTiedPointT2Unit import AbsTiedPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t2.T2PanStarrThumbPrint import T2PanStarrThumbPrint
from ampel.contrib.hu.util.ned import check_ned_res
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.PlotProperties import FormatModel, PlotProperties
from ampel.model.UnitModel import UnitModel
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.util.collections import ampel_iter
from ampel.view.T2DocView import T2DocView


class T2PS1ThumbNedTap(AbsTiedPointT2Unit):
    """
    This point t2 unit is tied with the point T2 unit T2NedTap.
    It retrieves panstarrs images at a datapoint location and for each NED catalog matching result:
    - creates a new image
    - marks the datapoint location
    - marks the matched location from the catalog

    If you use custom 'ingest options', please make sure that these are the same between T2NedTap and T2NedTapPS1ThumbPrint.
    Please note that the variant t2 class T2PS1ThumbNedSNCosmo exists and allows to restrict image retrieval
    based on SNCosmo convergence criteria.

    :param band: example: ["g", "r", "i", "z", "y"]
    """

    t2_dependency: Sequence[UnitModel[Literal["T2NedTap"]]]

    cmaps: Sequence[str] = ["cividis"]
    band: str | Sequence[str] = "g"
    plot_all: bool = False
    z_range: None | tuple[float, float]
    spectroscopic: bool = True

    plot_props: PlotProperties = PlotProperties(
        tags=["THUMBPRINT", "PANSTARRS"],
        file_name=FormatModel(
            format_str="%s_%s_%s_ps1_%s_thumb.svg",
            arg_keys=["stock", "obj_name", "band", "cmap"],
        ),
        title=FormatModel(
            format_str="%s - %s z=%s (%s band) ",
            arg_keys=["stock", "obj_name", "z", "band"],
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
        if not isinstance(cat_results, dict):
            return UnitResult(code=DocumentCode.T2_UNEXPECTED_DEPENDENCY)

        if "data" not in cat_results:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        if not isinstance(cat_results["data"], list):
            return 1

        plots = []

        for i, cat_res in enumerate(cat_results["data"]):
            if check_ned_res(cat_res, self.logger, self.spectroscopic, self.z_range):
                self.logger.info(f"Skipping cat result with index {i}")
                continue

            for band in ampel_iter(self.band):
                pt = T2PanStarrThumbPrint.get_ps1_target(datapoint, band)

                for cmap in ampel_iter(self.cmaps):
                    plots.append(
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
                                "cmap": cmap,
                                "obj_name": cat_res["prefname"].replace(" ", "_"),
                                "z": cat_res["z"],
                            },
                            logger=self.logger,
                        )
                    )

            if not self.plot_all:
                break

        return {"data": {"plots": plots}}
