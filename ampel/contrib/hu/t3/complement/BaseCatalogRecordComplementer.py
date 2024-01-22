#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/complement/BaseCatalogRecordComplementer.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                03.11.2020
# Date:                03.11.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>


from collections.abc import Iterable
from functools import cached_property
from typing import Any

from pymongo import MongoClient

from ampel.abstract.AbsBufferComplement import AbsBufferComplement
from ampel.base.decorator import abstractmethod
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.AmpelBuffer import AmpelBuffer
from ampel.struct.T3Store import T3Store


class BaseCatalogRecordComplementer(AbsBufferComplement, abstract=True):
    """
    Add entries from an extcats catalog to transients, matched by canonical
    name. This is distinct from, and much cheaper than, general catalog
    matching, where entries are matched by cone search, and can be used for
    dynamic "catalogs" like the TNS or instrument-specific marshals.
    """

    auth: NamedSecret[dict] = NamedSecret(label="extcats/reader")

    @cached_property
    def mongo_client(self):
        return MongoClient(
            self.context.config.get("resource.extcats", str), **self.auth.get()
        )

    @abstractmethod
    def get_catalog_item(self, names: tuple[str, ...]) -> None | dict[str, Any]:
        """Get catalog entry associated with the stock name"""
        ...

    def get_tag(self):
        """
        Key to use for extra items.
        """
        return self.__class__.__name__

    def complement(self, records: Iterable[AmpelBuffer], t3s: T3Store) -> None:
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
