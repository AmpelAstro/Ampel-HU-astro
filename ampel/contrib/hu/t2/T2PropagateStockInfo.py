#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:               ampel/contrib/hu/t2/T2PropagateStockInfo.py
# License:            BSD-3-Clause
# Author:             jnordin@physik.hu-berlin.de
# Date:               04.01.2022
# Last Modified Date: 04.01.2022
# Last Modified By:   jnordin@physik.hu-berlin.de

from collections.abc import Generator
from typing import Any

from ampel.abstract.AbsStockT2Unit import AbsStockT2Unit


def get_recursively(search_dict: dict[str,Any], field: str) -> Generator[str, None, None]:
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """

    for key, value in search_dict.items():
        if key == field:
            yield f"{key}_found_{value}"

        elif isinstance(value, dict):
            for result in get_recursively(value, field):
                yield f"{key}.{result}"

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    for result in get_recursively(item, field):
                        yield f"{key}.{result}"

class T2PropagateStockInfo(AbsStockT2Unit):
    """
    Collects a set of info from the transient stock collection and propagates.
    Potentially useful as chained T2.

    All entries are converted to strings.

    """

    # Paths to properties to search for
    # dict keys corresponds to output dict labels
    # dict values are reduce path list to intended values
    #    prop_paths: dict[str,list[str]] = {'explosion_time':['journal','healpix','time']}
    prop_paths: dict[str, list[str]]

    def process(self, stock_doc):
        outd = {}
        for label, pathlist in self.prop_paths.items():
            for path in get_recursively(stock_doc, pathlist[-1]):
                (dictpath, dictval) = path.split("_found_")
                if ".".join(pathlist) == dictpath:
                    outd[label] = dictval

        return outd
