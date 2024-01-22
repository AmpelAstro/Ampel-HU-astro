#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2PanStarrThumbPrint.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                31.01.2021
# Last Modified Date:  08.09.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from collections.abc import Sequence

from pymage.panstarrs import PS1Target  # type: ignore[import]

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.model.PlotProperties import FormatModel, PlotProperties
from ampel.plot.create import create_plot_record
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.util.collections import ampel_iter


class T2PanStarrThumbPrint(AbsPointT2Unit):
    """
    Retrieve panstarrs images at datapoint location and save
    these as compressed svg into the returned dict.

    :param band: example: ["g", "r", "i", "z", "y"]
    """

    cmaps: Sequence[str] = ["cividis"]
    band: str | Sequence[str] = "g"
    plot_props: PlotProperties = PlotProperties(
        tags=["THUMBPRINT", "PANSTARRS"],
        file_name=FormatModel(
            format_str="%s_%s_ps1_thumb.svg", arg_keys=["stock", "band"]
        ),
        title=FormatModel(format_str="%s (%s band) ", arg_keys=["stock", "band"]),
        id_mapper="ZTFIdMapper",
    )

    def process(self, datapoint: DataPoint) -> UBson | UnitResult:
        """ """

        pt = self.get_ps1_target(datapoint, self.band)
        stock_id = (
            datapoint["stock"]
            if isinstance(datapoint["stock"], int)
            else datapoint["stock"][0]
        )
        return {
            "data": {
                "plots": [
                    create_plot_record(
                        pt.show(
                            ellipse=False,
                            band=band,
                            show_target=False,
                            cmap=cmap,
                            show=False,
                        ),
                        self.plot_props,
                        extra={
                            "band": band,
                            "stock": stock_id,
                            "cmap": cmap,
                        },
                        logger=self.logger,
                    )
                    for cmap in ampel_iter(self.cmaps)
                    for band in ampel_iter(self.band)
                ]
            }
        }

    @staticmethod
    def get_ps1_target(datapoint: DataPoint, band: str | Sequence[str]) -> PS1Target:
        pt = PS1Target(None)
        pt.set_coordinate(datapoint["body"]["ra"], datapoint["body"]["dec"])
        pt.download_cutout(filters=band)
        return pt
