#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2NuclearFilter.py
# License:             BSD-3-Clause
# Author:              jannis.necker@gmail.com
# Date:                29.04.2026
# Last Modified Date:  29.04.2026
# Last Modified By:    jannis.necker@gmail.com

from collections.abc import Sequence
from typing import Literal, TypedDict

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import angular_separation
from scipy.stats import chi2

from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.util.catalog_column_info import (
    get_catalog_position_unit_map,
    get_type_and_redshift_columns,
)
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class NuclearFilterResult(TypedDict):
    passed: bool
    mean_ra: float
    mean_dec: float
    sep90: float
    sep1sig: float
    host_catalogs: list[str] | None
    host_ra: float | None
    host_dec: float | None
    host_dist_arcsec: float | None
    host_type: dict[str, dict[str, str]] | None
    template_flux: dict[str, tuple[float, float, float]] | None


class T2NuclearFilter(AbsTiedStateT2Unit):
    match_dist_arcsec: float
    group_matches_within_arcsec: float = 0.5

    t2_dependency: Sequence[
        StateT2Dependency[Literal["T2CatalogMatch", "T2LSPhotoZTap"]]
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # set up unit mapping for positions
        unit_mapping = get_catalog_position_unit_map()
        convert_to_rad = []
        for name, units in unit_mapping.items():
            if ((u := units["ra"]["unit"]) is None) or (u.lower().startswith("deg")):
                convert_to_rad.append(name)
        self._known_mapping = [*list(unit_mapping.keys()), "T2LSPhotoZTap"]
        self._convert_to_rad = [*convert_to_rad, "T2LSPhotoZTap"]
        self._redshift_columns, self._type_columns = get_type_and_redshift_columns()
        self._percentile_2dsig = chi2.cdf(1, 2)

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
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

        # calculate mean position and variance
        coords = SkyCoord(
            *np.array(
                [
                    [dp["body"][k] for k in ["ra", "dec"]]
                    for dp in datapoints
                    if ("diaSourceId" in dp["body"])
                ]
            ).T,
            unit="deg",
        )
        weights = np.array(
            [1] * len(coords)
        )  # TODO: figure weights out, (pos errors?, magnitiudes?)
        mean_pos = SkyCoord(
            (coords.represent_as("cartesian") * weights).sum() / sum(weights)
        )
        mean_pos.representation_type = "spherical"
        separations_to_mean = coords.separation(mean_pos)
        sep90 = np.quantile(separations_to_mean, 0.9).to_value("arcsec")
        sep1sig = np.quantile(separations_to_mean, self._percentile_2dsig).to_value(
            "arcsec"
        )

        mean_ra = mean_pos.ra.to_value("deg")
        mean_dec = mean_pos.dec.to_value("deg")

        # calculate template fluxes
        template_flux = {
            b: tuple(
                np.quantile(
                    [
                        dp["body"]["templateFlux"]
                        for dp in datapoints
                        if ("templateFlux" in dp) and (dp["body"]["band"] == b)
                    ],
                    [0.5, 0.05, 0.95],
                ).tolist()
            )
            for b in "ugrizy"
        }

        # If there are no matches within self.match_dist_arcsec, return False
        if not matches:
            return NuclearFilterResult(
                passed=False,
                mean_ra=mean_ra,
                mean_dec=mean_dec,
                sep90=sep90,
                sep1sig=sep1sig,
                host_ra=None,
                host_dec=None,
                host_dist_arcsec=None,
                host_catalogs=None,
                host_type=None,
                template_flux=template_flux,
            )

        # normalize keys
        matches = {
            k: {kk.lower(): vv for kk, vv in v.items()} for k, v in matches.items()
        }

        # PS1_photoz is special
        if "PS1_photoz" in matches:
            matches["PS1_photoz"]["ra"] = matches["PS1_photoz"]["ramean"]
            matches["PS1_photoz"]["dec"] = matches["PS1_photoz"]["decmean"]

        # remove matches with unknown units
        for m in [mm for mm in matches if mm not in self._known_mapping]:
            self.logger.info(f"Removing {m} from matches because unit is unknown")
            matches.pop(m)

        # convert to radians where necessary
        for name in matches:
            for k in ["ra", "dec"]:
                if k not in matches[name]:
                    raise KeyError(f"Key {k} not found in matches")
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
            separations <= self.group_matches_within_arcsec
        ]
        passed = bool(dist <= md)

        type_info = {
            k: {kk: v[kk] for kk in self._type_columns if kk in v}
            for k, v in matches.items()
        }

        if ("T2LSPhotoZTap" in type_info) and (
            type_info["T2LSPhotoZTap"]["type"] in {"DUP", "PSF"}
        ):
            passed = False

        return NuclearFilterResult(
            passed=passed,
            mean_ra=mean_ra,
            mean_dec=mean_dec,
            sep90=sep90,
            sep1sig=sep1sig,
            host_ra=np.degrees(match_map["ra_rad"][best_match_id]),
            host_dec=np.degrees(match_map["dec_rad"][best_match_id]),
            host_dist_arcsec=dist,
            host_catalogs=matched_catalogs.tolist(),
            host_type=type_info,
        )
