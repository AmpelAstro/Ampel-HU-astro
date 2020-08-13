#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/complement/AddTNSNames.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 13.12.2018
# Last Modified Date: 13.08.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>


from typing import Iterable

import numpy
from pydantic import Field

from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.core.AmpelBuffer import AmpelBuffer
from ampel.core.AmpelContext import AmpelContext
from ampel.t3.complement.AbsT3DataAppender import AbsT3DataAppender

from ..tns import TNSMirrorDB


class AddTNSNames(AbsT3DataAppender):
    """
    Add TNS names to transients
    """

    search_radius: float = Field(3, description="Matching radius in arcsec")

    def __init__(self, context: AmpelContext, **kwargs) -> None:

        AmpelBaseModel.__init__(self, **kwargs)

        self.tns = TNSMirrorDB(
            context.config.get("resource.extcats.writer"), logger=self.logger
        )

    def complement(self, records: Iterable[AmpelBuffer]) -> None:
        for record in filter(self.nameable, records):
            coords = [pp.get_tuple("ra", "dec") for pp in record.get("t0")]
            ra, dec = map(numpy.mean, zip(*coords))
            names = [
                "TNS" + str(n)
                for n in self.tns.get_names_for_location(ra, dec, self.search_radius)
            ]
            if not names:
                continue
            stock = record["stock"]
            if not stock.get("name"):
                stock["name"] = names
            else:
                stock["name"] = list(stock["names"]) + names

    @staticmethod
    def nameable(item: AmpelBuffer) -> bool:
        if not isinstance((stock := item.get("stock", None)), dict):
            return False
        if any(name.startswith("TNS") for name in stock.get("name", [])):
            return False
        if not item.get("datapoints"):
            return False
        return True
