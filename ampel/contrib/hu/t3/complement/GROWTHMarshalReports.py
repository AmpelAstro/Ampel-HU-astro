#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/complement/GROWTHMarshalReports.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 03.11.2020
# Date              : 03.11.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

from functools import cached_property
from typing import Any, Dict, Optional, Tuple

from ampel.contrib.hu.t3.complement.BaseCatalogRecordComplementer import (
    BaseCatalogRecordComplementer,
)
from ampel.contrib.hu.t3.tns.TNSName import TNSName


class GROWTHMarshalReports(BaseCatalogRecordComplementer):
    """
    Add GROWTH Marshal records from a local extcats mirror of the ProgramList.
    Though the GROWTH Marshal is no longer being updated, this is useful for
    looking up classifications of sources first discovered with ZTF I.
    """

    @cached_property
    def sources(self):
        return self.mongo_client.get_database("GROWTHMarshal").get_collection("srcs")

    def get_catalog_item(self, names: Tuple[str, ...]) -> Optional[Dict[str, Any]]:
        """Get catalog entry associated with the stock name"""
        for name in names:
            if name.startswith("ZTF") and (
                entry := self.sources.find_one(
                    {"name": name}, {"_id": False, "pos": False}
                )
            ):
                return entry
        return None
