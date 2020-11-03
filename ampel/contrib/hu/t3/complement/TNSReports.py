#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/complement/TNSReports.py
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


class TNSReports(BaseCatalogRecordComplementer):
    """
    Add TNS report from a local extcats mirror. This requires that TNSNames be
    run first.
    """

    @cached_property
    def sources(self):
        return self.mongo_client.get_database("TNS").get_collection("srcs")

    def get_catalog_item(self, names: Tuple[str, ...]) -> Optional[Dict[str, Any]]:
        """Get TNS entry associated with the first stock name that starts with 'TNS'"""
        for name in names:
            if name.startswith("TNS") and (
                entry := self.sources.find_one(
                    {"_id": int(TNSName.from_str(name[3:]))},
                    {"_id": False, "pos": False},
                )
            ):
                return entry
        return None
