#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/complement/BaseCatalogRecordComplementer.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 03.11.2020
# Date              : 03.11.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>


from functools import cached_property
from typing import Any, Dict, Iterable, Optional, Tuple

from pymongo import MongoClient

from ampel.base import abstractmethod
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.core.AmpelBuffer import AmpelBuffer
from ampel.model.Secret import Secret
from ampel.t3.complement.AbsT3DataAppender import AbsT3DataAppender


class BaseCatalogRecordComplementer(AbsT3DataAppender, abstract=True):
    """
    Add entries from an extcats catalog to transients, matched by canonical
    name. This is distinct from, and much cheaper than, general catalog
    matching, where entries are matched by cone search, and can be used for
    dynamic "catalogs" like the TNS or instrument-specific marshals.
    """

    auth: Secret[dict] = {"key": "extcats/reader"}  # type: ignore[assignment]

    @cached_property
    def mongo_client(self):
        return MongoClient(
            self.context.config.get(f"resource.extcats", str), **self.auth.get()
        )

    @abstractmethod
    def get_catalog_item(self, names: Tuple[str, ...]) -> Optional[Dict[str, Any]]:
        """Get catalog entry associated with the stock name"""
        ...

    def get_tag(self):
        """
        Key to use for extra items. Default assumes that class name ends in s.
        """
        return self.__class__.__name__[:-1]

    def complement(self, records: Iterable[AmpelBuffer]) -> None:
        for record in records:
            if (stock := record["stock"]) is None:
                raise ValueError(f"{type(self).__name__} requires stock records")
            item = self.get_catalog_item(
                tuple(name for name in (stock["name"] or []) if isinstance(name, str))
            )
            if record.get("extra") is None or record["extra"] is None:
                record["extra"] = {self.get_tag(): item}
            else:
                record["extra"][self.get_tag()] = item
