#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/complement/AddTNSNames.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 13.12.2018
# Last Modified Date: 13.08.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>


from typing import Iterable, Tuple, Union

import numpy
from pydantic import Field
from pymongo import MongoClient

from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.contrib.hu.t3.tns.TNSMirrorDB import TNSMirrorDB
from ampel.core.AmpelBuffer import AmpelBuffer
from ampel.core.AmpelContext import AmpelContext
from ampel.model.Secret import Secret
from ampel.t3.complement.AbsT3DataAppender import AbsT3DataAppender
from ampel.view.ReadOnlyDict import ReadOnlyDict


class AddTNSNames(AbsT3DataAppender):
    """
    Add TNS names to transients
    """

    search_radius: float = Field(3, description="Matching radius in arcsec")
    auth: Secret[dict] = {"key": "extcats/reader"}  # type: ignore[assignment]

    def __init__(self, context: AmpelContext, **kwargs) -> None:

        AmpelBaseModel.__init__(self, **kwargs)

        self.tns = TNSMirrorDB(
            MongoClient(
                context.config.get(f"resource.extcats", str), **self.auth.get()
            ),
            logger=self.logger,
        )

    def complement(self, records: Iterable[AmpelBuffer]) -> None:
        for record in records:
            if not isinstance((stock := record["stock"]), dict):
                continue
            if (stock_name := stock["name"]) and any(
                (
                    isinstance(name, str) and name.startswith("TNS")
                    for name in stock_name
                )
            ):
                continue
            if not (photopoints := record["t0"]):
                continue
            coords = [
                (pp["body"]["ra"], pp["body"]["dec"])
                for pp in photopoints
                if pp["_id"] > 0
            ]
            ra, dec = map(numpy.mean, zip(*coords))
            names: Tuple[Union[int, str], ...] = tuple(
                "TNS" + str(n)
                for n in self.tns.get_names_for_location(ra, dec, self.search_radius)
            )
            if not names:
                continue
            elif stock_name:
                names = tuple(stock_name) + names
            record["stock"] = ReadOnlyDict({**stock, **{"name": names}})  # type: ignore[arg-type]
