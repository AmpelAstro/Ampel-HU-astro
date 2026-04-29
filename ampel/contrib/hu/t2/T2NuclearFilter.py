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
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2NuclearFilter(AbsTiedPointT2Unit):
    match_dist_arcsec: float
    group_matches_within_arcsec: float = 0.5

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

        # find the closest match
        match_map = np.array(
            [(k, v["dist2transient"], v["ra"], v["dec"]) for k, v in matches.items()]
        )
        best_match_id = np.argmin(match_map[:, 1])
        dist = matches[best_match_id, 1]

        # find matches that are probably the same object
        separations = angular_separation(
            np.radians(match_map[best_match_id, 2]),
            np.radians(match_map[best_match_id, 3]),
            np.radians(match_map[:, 2]),
            np.radians(match_map[:, 3]),
        )
        matched_catalogs = match_map[separations < self.group_matches_within_arcsec, 0]
        res["closest_matches"] = matched_catalogs.tolist()
        res["pass"] = dist <= md
        res["dist"] = dist

        return res
