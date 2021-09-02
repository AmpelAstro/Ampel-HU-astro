#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t0/LensedTransientFilter.py
# License           : BSD-3-Clause
# Author            : m. giomi <matteo.giomi@desy.de>
# Date              : 04.27.2018
# Last Modified Date: 19.03.2021
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

from ampel.ztf.base.CatalogMatchFilter import CatalogMatchFilter

class LensedTransientFilter(CatalogMatchFilter):

    ClusListSearchRadius: float
    MasterlensSearchRadius: float
    CastleQSOSearchRadius: float

    def __init__(self, **kwargs):
        kwargs["accept"] = {
            "any_of": [
                {
                    "name": catalog.lower(),
                    "use": "extcats",
                    "rs_arcsec": kwargs.get(f"{catalog}SearchRadius"),
                }
                for catalog in ("ClusList", "Masterlens", "CastleQSO")
            ]
        }
        super().__init__(**kwargs)
