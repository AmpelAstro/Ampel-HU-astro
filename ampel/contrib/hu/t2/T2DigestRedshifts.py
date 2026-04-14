#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2DigestRedshifts.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                06.06.2021
# Last Modified Date:  19.10.2022
# Last Modified By:    atownsend@physik.hu-berlin.de

from collections.abc import Mapping, Sequence
from typing import Any, Literal

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.util import get_payload
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2DigestRedshifts(AbsTiedStateT2Unit):
    """

    Compare potential matches from different T2 units providing redshifts.

    Using table comparisons from (T3) CompareRedshifts to select best matches
    and provide a general uncertainty estimate.

    Available (studied) redshifts are assigned to one of seven redshift "groups",
    with decreasing average quality. The mean redshift of the lowest populated
    group is returned as 'ampel_z', together with info about this group
    ('group_z_nbr' and 'group_z_precision').

    (This is for now an AbsTiedLightCurveT2Unit, but lc info not used.)

    """

    # Max redshift uncertainty category: 1-7
    # (where 7 is any, and 1 only nearby spectroscopic matches)
    max_redshift_category: int = 3

    # Redshift estimates associated with each region ( only rough guideline!!! )
    category_precision: list[float] = [0.0003, 0.003, 0.01, 0.02, 0.04, 0.1, 0.3]

    # Maximum physical host offset in kpc
    max_host_offset_kpc: float = 25.0

    # Cosmology used to convert kpc <-> arcsec
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # CatalogMatch(Local) results might be overriden,
    # for example if specialized catalog is being used
    # Each override dict is assumed to be built asmed to be built according to
    # "catalog_name" : {
    #                    "z_keyword": "redshift field in catalog",
    #                    "max_distance": "max arcsec in which to allow match,
    #                    "max_redshift": "max redshift to use",
    #                    "z_group": "which redshift group to assign to" }
    catalogmatch_override: None | dict[str, Any]

    # Options for the get_redshift option
    # T2MatchBTS : Use the redshift published by BTS and  synced by that T2.
    # T2DigestRedshifts : Use the best redshift as parsed by DigestRedshift.
    # AmpelZ: equal to T2DigestRedshifts
    # T2ElasticcRedshiftSampler: Use a list of redshifts and weights from the sampler.
    # None : Use the fixed z value
    redshift_kind: (
        None
        | Literal[
            "T2MatchBTS", "T2DigestRedshifts", "T2ElasticcRedshiftSampler", "AmpelZ"
        ]
    ) = None

    # It is also possible to use fixed redshift whenever a dynamic redshift kind is not possible
    # This could be either a single value or a list
    fixed_z: None | float | Sequence[float] = None
    # Finally, the provided lens redshift might be multiplied with a scale
    # Useful for lensing studies, or when trying multiple values
    scale_z: None | float = None

    # These are the units through which we look for redshifts
    # Which units should this be changed to
    t2_dependency: Sequence[
        StateT2Dependency[
            Literal[
                "T2CatalogMatch", "T2LSPhotoZTap", "T2CatalogMatchLocal", "T2MatchBTS"
            ]
        ]
    ]

    def _is_close_match(self, dist: float, z: float, half: bool = False) -> bool:
        """
        Distance/redshift cut based on a maximum physical host offset.
        half=True corresponds to half the full cut.
        """

        if z <= 0 or dist < 0:
            return False

        max_offset_kpc = self.max_host_offset_kpc / 2 if half else self.max_host_offset_kpc

        kpc_per_arcsec = self.cosmo.kpc_proper_per_arcmin(z).value / 60.0
        max_dist_arcsec = max_offset_kpc / kpc_per_arcsec

        return dist < max_dist_arcsec

    def _append_group_with_unc(
        self,
        group_z: list[list[float]],
        group_dist: list[list[float]],
        group_unc: list[list[float]],
        z: float,
        dist: float,
        z_unc: float,
        close_group: int,
        medium_group: int,
        far_group: int = 6,
    ) -> None:
        """
        Append a redshift/distance/uncertainty triple to the appropriate quality group.
        Group indices are zero-based.
        """
        if self._is_close_match(dist, z, half=True):
            group_z[close_group].append(z)
            group_dist[close_group].append(dist)
            group_unc[close_group].append(z_unc)
        elif self._is_close_match(dist, z, half=False):
            group_z[medium_group].append(z)
            group_dist[medium_group].append(dist)
            group_unc[medium_group].append(z_unc)
        else:
            group_z[far_group].append(z)
            group_dist[far_group].append(dist)
            group_unc[far_group].append(z_unc)

    def _normalize_unc(self, unc: None | float) -> None | float:
        """
        Normalize redshift uncertainty for weighting.

        In some catalogs (for example NED2026), zunc=0 can occur for what is effectively
        a negligibly small uncertainty. In this context, treat this as a very small
        uncertainty, i.e. approximately maximal weight.
        """

        if unc is None:
            return None

        if unc == 0:
            return 1e-6

        if unc < 0:
            return None

        return float(unc)

    def _weighted_mean(self, values: Sequence[float], uncs: Sequence[float]) -> float:
        """
        Compute inverse-variance weighted mean using per-measurement uncertainties.
        """

        if len(values) == 0:
            raise ValueError("Cannot compute weighted mean of empty sequence")

        if len(values) != len(uncs):
            raise ValueError("values and uncs must have the same length")

        valid_pairs = [
            (val, norm_unc)
            for val, unc in zip(values, uncs, strict=False)
            if (norm_unc := self._normalize_unc(unc)) is not None
        ]
        if len(valid_pairs) == 0:
            raise ValueError("No usable uncertainties available for weighted mean")

        valid_values = np.array([val for val, _ in valid_pairs], dtype=float)
        weights = np.array(
            [1.0 / (unc * unc) for _, unc in valid_pairs],
            dtype=float,
        )

        return float(np.average(valid_values, weights=weights))

    def _combined_weighted_unc(self, uncs: Sequence[float]) -> float:
        """
        Return uncertainty of the inverse-variance weighted mean.
        """

        valid_uncs = [
            norm_unc
            for unc in uncs
            if (norm_unc := self._normalize_unc(unc)) is not None
        ]
        if len(valid_uncs) == 0:
            raise ValueError("No usable uncertainties available")

        weights = np.array([1.0 / (unc * unc) for unc in valid_uncs], dtype=float)
        return float(np.sqrt(1.0 / np.sum(weights)))

    def _get_lsphotoz_groupz(
        self, t2_res: Mapping[str, Any]
    ) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
        """
        Parse output from T2LSPhotoZTap and investigate whether any matches fulfill group
        redshift criteria.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]
        group_dist: list[list[float]] = [[], [], [], [], [], [], []]
        group_unc: list[list[float]] = [[], [], [], [], [], [], []]
        for lsdata in t2_res.values():
            if lsdata is None:
                continue

            # Warning: all LS checks done with a 10" matching radius,
            # this is thus enforced (in case T2 run with larger radius)
            if lsdata["dist2transient"] > 10:
                self.logger.debug(
                    "No Digest redshift LS estimate.",
                    extra={"dist2transient": lsdata["dist2transient"]},
                )
                continue

            # LS Spec is included in NED, so only photometric redshifts
            if lsdata["z_phot_median"] is not None and lsdata["z_phot_median"] > -1:
                dist = lsdata["dist2transient"]
                z = lsdata["z_phot_median"]
                z_unc = self._normalize_unc(lsdata.get("z_phot_std"))
                if z_unc is not None:
                    self._append_group_with_unc(
                        group_z, group_dist, group_unc, z, dist, z_unc, 2, 5
                    )
            self.logger.debug(f"LS debug phot: {lsdata} yield {group_z}")

        return group_z, group_dist, group_unc

    def _get_catalogmatch_groupz(
        self, t2_res: Mapping[str, Any]
    ) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
        """
        Parse output from T2CatalogMatch.

        Made complicated as returns can be both single and lists.


        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]
        group_dist: list[list[float]] = [[], [], [], [], [], [], []]
        group_unc: list[list[float]] = [[], [], [], [], [], [], []]

        for cat_name, cat_matches in t2_res.items():
            if cat_matches is None or cat_matches is False:
                continue
            # List or dict depending on whether the closest or all matches are returned from.
            if isinstance(cat_matches, list):
                cat_match_list = cat_matches
            elif isinstance(cat_matches, tuple):
                cat_match_list = list(cat_matches)
            else:
                cat_match_list = [cat_matches]

            for cat_match in cat_match_list:
                # All catalogs have different structure, so doing this individually

                if cat_name == "NED2026" and cat_match["z"] is not None:
                    dist = cat_match["dist2transient"]
                    z = cat_match["z"]
                    z_unc = self._normalize_unc(cat_match.get("zunc"))
                    if z_unc is not None:
                        self._append_group_with_unc(
                            group_z, group_dist, group_unc, z, dist, z_unc, 0, 3
                        )

                if cat_name == "LSPhotoZZou":
                    dist = cat_match["dist2transient"]

                    # Spec
                    if cat_match["specz"] is not None and cat_match["specz"] > -0.1:
                        z = float(cat_match["specz"])
                        z_unc = self._normalize_unc(1e-6)
                        if z_unc is not None:
                            self._append_group_with_unc(
                                group_z, group_dist, group_unc, z, dist, z_unc, 1, 4
                            )

                    # Photo-z
                    if cat_match["photoz"] is not None and cat_match["photoz"] > -0.1:
                        z = float(cat_match["photoz"])
                        z_unc = self._normalize_unc(cat_match.get("e_photoz"))
                        if z_unc is not None:
                            self._append_group_with_unc(
                                group_z, group_dist, group_unc, z, dist, z_unc, 2, 5
                            )


                if cat_name == "PS1_photoz" and cat_match["z_phot"] is not None:
                    z = float(cat_match["z_phot"]) / 1000
                    dist = cat_match["dist2transient"]
                    z_unc_raw = cat_match.get("z_photErr")

                    if z > -0.1 and z_unc_raw is not None:
                        z_unc = self._normalize_unc(float(z_unc_raw) / 1000)
                        if z_unc is not None:
                            self._append_group_with_unc(
                                group_z, group_dist, group_unc, z, dist, z_unc, 2, 5
                            )

                # Also check for manual override
                if self.catalogmatch_override:
                    for or_catname, or_catdict in self.catalogmatch_override.items():
                        if or_catname == cat_name:
                            try:
                                cat_z = float(cat_match[or_catdict["z_keyword"]])
                                cat_unc = self._normalize_unc(
                                    float(cat_match[or_catdict["z_unc_keyword"]])
                                )
                                if (
                                    float(cat_match["dist2transient"])
                                    < or_catdict["max_distance"]
                                    and cat_z < or_catdict["max_redshift"]
                                    and cat_unc is not None
                                ):
                                    group_idx = or_catdict["z_group"] - 1
                                    group_z[group_idx].append(cat_z)
                                    group_dist[group_idx].append(
                                        float(cat_match["dist2transient"])
                                    )
                                    group_unc[group_idx].append(cat_unc)
                            except (ValueError, KeyError, TypeError):
                                self.logger.info(
                                    "Cannot parse z", extra={"catdict": cat_match}
                                )

        return group_z, group_dist, group_unc

    def _get_matchbts_groupz(
        self, t2_res: Mapping[str, Any]
    ) -> tuple[list[list[float]], list[list[float]]]:
        """
        Parse output from T2MatchBTS.

        Any transient with a redshift with two decimals (from SN template matching) is put in Group II,
        those with more (from host, high-res spec) are put into Group I.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]
        group_unc: list[list[float]] = [[], [], [], [], [], [], []]

        if (
            isinstance(bts_redshift := t2_res.get("bts_redshift"), str)
            and bts_redshift != "-"
        ):
            # BTS redshifts are stored as strings. Crude way to get to redshift precision for evaluation:
            # Take decimal part, remove initial zeroes and cound digits
            decimals = len(bts_redshift.split(".")[1].lstrip("0"))
            z = float(bts_redshift)

            if decimals > 2:
                group_z[0].append(z)
                group_unc[0].append(10 ** (-decimals))
            else:
                group_z[1].append(z)
                group_unc[1].append(10 ** (-decimals))

        self.logger.debug(f" bts match yield {group_z}")

        return group_z, group_unc

    def get_ampelZ(self, t2_views: Sequence[T2DocView]) -> dict[str, UBson]:
        """

        Parse t2_views from catalogs that were part of the redshift studies.
        Return these together with a "best estimate" - ampel_z

        Main method, separated to be used externally.

        """

        # Loop through all potential T2s with redshift information.
        # Each should return an array of arrays, corresponding to redshift maches
        # found in each category. These will be added to sn redshifts
        group_redshifts: list[list[float]] = [[], [], [], [], [], [], []]
        group_distances: list[list[float]] = [[], [], [], [], [], [], []]
        group_redshift_uncs: list[list[float]] = [[], [], [], [], [], [], []]

        # Loop through t2_views and collect information.
        for t2_view in t2_views:
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            t2_res = get_payload(t2_view)

            if t2_view.unit == "T2LSPhotoZTap":
                new_zs, new_dists, new_uncs = self._get_lsphotoz_groupz(t2_res)
            elif t2_view.unit in {"T2CatalogMatch", "T2CatalogMatchLocal"}:
                new_zs, new_dists, new_uncs = self._get_catalogmatch_groupz(t2_res)
            elif t2_view.unit == "T2MatchBTS":
                new_zs, new_uncs = self._get_matchbts_groupz(t2_res)
                new_dists = [[], [], [], [], [], [], []]
            else:
                self.logger.error(f"No instructions for dealing with {t2_view.unit}")
                return {}

            for k in range(7):
                if len(new_zs[k]) > 0:
                    group_redshifts[k].extend(new_zs[k])
                if len(new_dists[k]) > 0:
                    group_distances[k].extend(new_dists[k])
                if len(new_uncs[k]) > 0:
                    group_redshift_uncs[k].extend(new_uncs[k])
            self.logger.debug(f"group_z after {t2_view.unit}: {group_redshifts}")

        # Check for best match
        t2_output: dict[str, UBson] = {
            "group_zs": group_redshifts,
            "group_dists": group_distances,
            "group_z_uncs": group_redshift_uncs,
        }
        for k in range(7):
            if (k + 1) > self.max_redshift_category:
                # No matches with sufficient precision
                break
            if len(group_redshifts[k]) > 0 and len(group_redshift_uncs[k]) > 0:
                t2_output["ampel_z"] = self._weighted_mean(
                    group_redshifts[k], group_redshift_uncs[k]
                )
                t2_output["group_z_precision"] = self._combined_weighted_unc(
                    group_redshift_uncs[k]
                )
                t2_output["group_z_nbr"] = k + 1
                if len(group_distances[k]) == len(group_redshift_uncs[k]) and len(group_distances[k]) > 0:
                    t2_output["ampel_dist"] = self._weighted_mean(
                        group_distances[k], group_redshift_uncs[k]
                    )
                elif len(group_distances[k]) > 0:
                    t2_output["ampel_dist"] = float(
                        np.mean(group_distances[k])
                    )  # distance transient to matched z sources NOT TO EARTH
                # We then do *not* look for higher group (more uncertain) matches
                break
        if self.catalogmatch_override:
            t2_output["AmpelZ-Warning"] = "Override catalog in use."

        return t2_output

    def get_hostCol(self, t2_views: Sequence[T2DocView]) -> dict:
        """

        Parse t2_views from catalogs and extract host galaxy color information.

        If a Wise catalog match within 3" is found, a W1-W2 color is returned.
        If a PS1 matchin within 10", general host colors are returned from there.

        """

        # Extracted information added here
        host_info: dict[str, UBson] = {}

        # Loop through t2_views and collect information.
        for t2_view in t2_views:
            if t2_view.unit not in {"T2CatalogMatch", "T2CatalogMatchLocal"}:
                continue
            t2_res = get_payload(t2_view)

            for cat_name, cat_matches in t2_res.items():
                if cat_matches is None or cat_matches is False:
                    continue
                # List or dict depending on whether the closest or all matches are returned from.
                if isinstance(cat_matches, list):
                    cat_match_list = cat_matches
                elif isinstance(cat_matches, tuple):
                    cat_match_list = list(cat_matches)
                else:
                    cat_match_list = [cat_matches]

                for cat_match in cat_match_list:
                    # All catalogs have different structure, so doing this individually

                    if cat_name == "WISE" and cat_match["dist2transient"] < 3:
                        # Only W1-W2 color
                        if (
                            cat_match.get("Mag_W1", -1) > 0
                            and cat_match.get("Mag_W2", -1) > 0
                        ):
                            host_info["col_wise_w1w2"] = (
                                cat_match["Mag_W1"] - cat_match["Mag_W2"]
                            )
                    elif cat_name == "PS1" and cat_match["dist2transient"] < 10:
                        # Which color to use? skipping y for now
                        if (
                            cat_match.get("gPSFMag", -1) > 0
                            and cat_match.get("rPSFMag", -1) > 0
                        ):
                            host_info["col_ps1_gr"] = (
                                cat_match["gPSFMag"] - cat_match["rPSFMag"]
                            )
                        if (
                            cat_match.get("iPSFMag", -1) > 0
                            and cat_match.get("rPSFMag", -1) > 0
                        ):
                            host_info["col_ps1_ri"] = (
                                cat_match["rPSFMag"] - cat_match["iPSFMag"]
                            )
                        if (
                            cat_match.get("iPSFMag", -1) > 0
                            and cat_match.get("zPSFMag", -1) > 0
                        ):
                            host_info["col_ps1_iz"] = (
                                cat_match["iPSFMag"] - cat_match["zPSFMag"]
                            )

        return host_info

    def get_redshift(
        self, t2_views
    ) -> tuple[None | list[float], None | str, None | list[float]]:
        """

        Return a single or list of redshifts to be used. Not called in T2DigestRedshift.process
        but provides interface to e.g. fit units.

        """

        # Examine T2s for eventual information
        z: None | list[float] = None
        z_source: None | str = None
        z_weights: None | list[float] = None

        t2_res: UBson | dict[str, Any] = {}
        if self.redshift_kind in [
            "T2DigestRedshifts",
            "AmpelZ",
        ]:
            t2_res = self.get_ampelZ(t2_views)
            if (
                isinstance(t2_res, dict)
                and "ampel_z" in t2_res
                and t2_res["ampel_z"] is not None
                and t2_res["group_z_nbr"] <= self.max_redshift_category
            ):
                z = [float(t2_res["ampel_z"])]
                z_source = "AMPELz_group" + str(t2_res["group_z_nbr"])
        elif self.redshift_kind in [
            "T2MatchBTS",
            "T2DigestRedshifts",
            "T2ElasticcRedshiftSampler",
        ]:
            for t2_view in t2_views:
                if t2_view.unit != self.redshift_kind:
                    continue
                self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                t2_res = get_payload(t2_view)
                # Parse this
                if self.redshift_kind == "T2MatchBTS":
                    if (
                        isinstance(t2_res, dict)
                        and "bts_redshift" in t2_res
                        and t2_res["bts_redshift"] != "-"
                    ):
                        z = [float(t2_res["bts_redshift"])]
                        z_source = "BTS"
                elif self.redshift_kind == "T2ElasticcRedshiftSampler" and isinstance(
                    t2_res, dict
                ):
                    z = t2_res["z_samples"]
                    z_source = t2_res["z_source"]
                    z_weights = t2_res["z_weights"]
        # Check if there is a fixed z set for this run, otherwise keep as free parameter
        elif self.fixed_z is not None:
            if isinstance(self.fixed_z, float):
                z = [self.fixed_z]
            else:
                z = list(self.fixed_z)
            z_source = "Fixed"
        else:
            z = None
            z_source = "Fitted"

        if (z is not None) and (z_source is not None) and self.scale_z:
            z = [onez * self.scale_z for onez in z]
            z_source += f" + scaled {self.scale_z}"

        return z, z_source, z_weights

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """

        Parse t2_views from catalogs that were part of the redshift studies.
        Return these together with a "best estimate" - ampel_z

        """

        if not t2_views:  # Should not happen actually, T2Processor catches that case
            self.logger.error("Missing tied t2 views")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        infocol = self.get_hostCol(t2_views)
        ampelz = self.get_ampelZ(t2_views)
        if ampelz is not None:
            infocol.update(ampelz)
        return infocol