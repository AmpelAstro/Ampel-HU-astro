#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:               ampel/contrib/hu/t2/T2PropagateStockInfo.py
# License:            BSD-3-Clause
# Author:             jnordin@physik.hu-berlin.de
# Date:               04.01.2022
# Last Modified Date: 04.01.2022
# Last Modified By:   jnordin@physik.hu-berlin.de

from functools import reduce
from ampel.abstract.AbsStockT2Unit import AbsStockT2Unit

def get_recursively(search_dict, field):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    """
    fields_found = []

    for key, value in search_dict.items():

        if key == field:
            fields_found.append(key+'_found_'+value)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(key+'.'+result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(key+'.'+another_result)

    return fields_found


class T2PropagateStockInfo(AbsStockT2Unit):
    """
    Collects a set of info from the transient stock collection and propagates.
    Potentially useful as chained T2
    """

    # Paths to properties to search for
    # dict keys corresponds to output dict labels
    # dict values are reduce path list to intended values
#    prop_paths: dict[str,list[str]] = {'explosion_time':['journal','healpix','time']}
    prop_paths: dict[str,list[str]]



    def process(self, stock_doc):
        print(stock_doc)
        outd = {}
        for label, pathlist in self.prop_paths.items():
            paths = get_recursively(stock_doc, pathlist[-1])
            for path in paths:
                (dictpath, dictval) = path.split('_found_')
                if '.'.join(pathlist)==dictpath:
                    outd[label] = dictval

        return outd
