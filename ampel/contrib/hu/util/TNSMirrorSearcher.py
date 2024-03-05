#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/util/TNSMirrorSearcher.py
# License:             BSD-3-Clause
# Author:              jnordin
# Date:                29.02.2024
# Last Modified Date:  29.02.2024
# Last Modified By:    jnordin <jnordin@physik.hu-berlin.de>

from ampel.ztf.base.CatalogMatchUnit import CatalogMatchUnit


class TNSMirrorSearcher(CatalogMatchUnit):
    """
    Methods for checking whether TNS data, grabbed from the DESY mirror,
    exists for specific transients.
    """

    # DESY TNS matching keys
    search_radius: float = 3.0

    def tns_conesearch(self, ra: float, dec: float, keys_to_append=None) -> list:
        return self.cone_search_all(
            ra,
            dec,
            [
                {
                    "name": "TNS",
                    "use": "extcats",
                    "rs_arcsec": self.search_radius,
                    "keys_to_append": keys_to_append,
                }
            ],
        )

    def get_tns_name(self, ra: float, dec: float) -> None | str:
        if (
            matches := self.tns_conesearch(ra, dec, keys_to_append=["objname"])
        ) and matches[0] is not None:
            return matches[0][0]["body"]["objname"]
        return None

    def get_tns_reports(self, ra: float, dec: float) -> None | list:
        if matches := self.tns_conesearch(ra, dec):
            return matches
        return None

    def get_tns_name_internal(
        self, ra: float, dec: float
    ) -> tuple[str | None, list[str]]:
        # Return a TNS name if existing + any other internal names
        if matches := self.tns_conesearch(ra, dec):
            if matches[0] is None:
                return (None, [])
            tnsname = matches[0][0]["body"]["objname"]
            names = matches[0][0]["body"]["internal_names"].split(", ")
            return (tnsname, names)
        return (None, [])
