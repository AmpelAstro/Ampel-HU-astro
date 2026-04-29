#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2NuclearFilter.py
# License:             BSD-3-Clause
# Author:              jannis.necker@gmail.com
# Date:                29.04.2026
# Last Modified Date:  29.04.2026
# Last Modified By:    jannis.necker@gmail.com

from collections.abc import Sequence

import numpy as np
from astropy.coordinates.angles import angular_separation

from ampel.abstract.AbsTiedPointT2Unit import AbsTiedPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.util.catalog_match_position_units import (
    get_catalog_position_unit_map,
)
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2NuclearFilter(AbsTiedPointT2Unit):
    match_dist_arcsec: float
    group_matches_within_arcsec: float = 0.5

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # set up unit mapping for positions
        unit_mapping = get_catalog_position_unit_map()
        convert_to_rad = []
        for name, units in unit_mapping.items():
            if ((u := units["ra"]) is None) or (u.lower().startswith("deg")):
                convert_to_rad.append(name)
        self._known_mapping = list(unit_mapping.values())
        self._convert_to_rad = convert_to_rad

    def process(
        self, datapoint: DataPoint, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:
        md = self.match_dist_arcsec

        # Collect all matches from the match units that are within self.match_dist_arcsec
        # Attention! If there are multiple catalog match units looking at the same catalogs
        # with different configurations the matches will be overwritten!
        matches = {}
        for t2_view in t2_views:
            if (t2_view.unit in {"T2CatalogMatch", "T2LSPhotoZTap"}) and (
                body := t2_view.get_payload()
            ) is not None:
                matches.update(
                    {
                        k: v
                        for k, v in body.items()
                        if (v is not None) and (v["dist2transient"] <= md)
                    }
                )

        # set up result structure
        res = {"closest_matches": None, "pass": False, "dist": None}

        # If there are no matches within self.match_dist_arcsec, return False
        if not matches:
            return res

        # normalize keys
        matches = {
            k: {kk.lower(): vv for kk, vv in v.items()} for k, v in matches.items()
        }

        # PS1_photoz is special
        if "PS1_photoz" in matches:
            matches["PS1_photoz"]["ra"] = matches["PS1_photoz"]["ramean"]
            matches["PS1_photoz"]["dec"] = matches["PS1_photoz"]["decmean"]

        # remove matches with unknown units
        for m in matches:
            if m not in self._known_mapping:
                self.logger.info(f"Removing {m} from matches because unit is unknown")
                matches.pop(m)

        # convert to radians where necessary
        for name in matches:
            for k in ["ra", "dec"]:
                if name in self._convert_to_rad:
                    matches[name][k] = np.radians(float(matches[name][k]))
                else:
                    matches[name][k] = float(matches[name][k])

        # find the closest match
        try:
            match_map = np.array(
                [
                    (k, float(v["dist2transient"]), float(v["ra"]), float(v["dec"]))
                    for k, v in matches.items()
                ],
                dtype=[
                    ("name", "<U30"),
                    ("distance", "<f8"),
                    ("ra_rad", "<f8"),
                    ("dec_rad", "<f8"),
                ],
            )

        except KeyError as e:
            raise e

        best_match_id = np.argmin(match_map["distance"])
        dist = match_map["distance"][best_match_id]

        # find matches that are probably the same object
        separations = (
            np.degrees(
                angular_separation(
                    match_map["ra_rad"][best_match_id],
                    match_map["dec_rad"][best_match_id],
                    match_map["ra_rad"],
                    match_map["dec_rad"],
                )
            )
            * 3600
        )
        matched_catalogs = match_map["name"][
            separations < np.radians(self.group_matches_within_arcsec / 3600)
        ]
        res["closest_matches"] = matched_catalogs.tolist()
        res["pass"] = dist <= md
        res["dist"] = dist

        return res
