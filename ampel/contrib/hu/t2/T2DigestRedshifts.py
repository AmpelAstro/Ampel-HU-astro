#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2DigestRedshifts.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                06.06.2021
# Last Modified Date:  19.10.2022
# Last Modified By:    atownsend@physik.hu-berlin.de

from collections.abc import Sequence
from typing import Any, Literal, cast

import numpy as np

from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView


class T2DigestRedshifts(AbsTiedLightCurveT2Unit):
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
    redshift_kind: None | Literal[
        "T2MatchBTS", "T2DigestRedshifts", "T2ElasticcRedshiftSampler", "AmpelZ"
    ] = None

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

    def _get_lsphotoz_groupz(
        self, t2_res: dict[str, Any]
    ) -> tuple[list[list[float]], list[list[float]]]:
        """
        Parse output from T2LSPhotoZTap and investigate whether any matches fulfill group
        redshift criteria.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]
        group_dist: list[list[float]] = [[], [], [], [], [], [], []]
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

            # First investigate LS spectroscopic redshift
            if lsdata["z_spec"] is not None and lsdata["z_spec"] > -1:
                if lsdata["z_spec"] < 0.03:
                    # Group I
                    group_z[0].append(lsdata["z_spec"])
                    group_dist[0].append(lsdata["dist2transient"])
                elif lsdata["z_spec"] < 0.15:
                    # Group II
                    group_z[1].append(lsdata["z_spec"])
                    group_dist[1].append(lsdata["dist2transient"])
                elif lsdata["z_spec"] < 0.4:
                    # Group III
                    group_z[2].append(lsdata["z_spec"])
                    group_dist[2].append(lsdata["dist2transient"])
                else:
                    # Group V
                    group_z[4].append(lsdata["z_spec"])
                    group_dist[4].append(lsdata["dist2transient"])
            self.logger.debug(f"LS debug spec: {lsdata} yield {group_z}")

            # Now, photometric redshifts
            if lsdata["z_phot_median"] is not None and lsdata["z_phot_median"] > -1:
                if lsdata["z_phot_median"] < 0.1:
                    # Group IV
                    group_z[3].append(lsdata["z_phot_median"])
                    group_dist[3].append(lsdata["dist2transient"])
                elif lsdata["z_phot_median"] < 0.2:
                    # Group V
                    group_z[4].append(lsdata["z_phot_median"])
                    group_dist[4].append(lsdata["dist2transient"])
                elif lsdata["z_phot_median"] < 0.4:
                    # Group VI
                    group_z[5].append(lsdata["z_phot_median"])
                    group_dist[5].append(lsdata["dist2transient"])
                else:
                    # Group VII
                    group_z[6].append(lsdata["z_phot_median"])
                    group_dist[6].append(lsdata["dist2transient"])
            self.logger.debug(f"LS debug phot: {lsdata} yield {group_z}")

        return group_z, group_dist

    def _get_catalogmatch_groupz(
        self, t2_res: dict[str, Any]
    ) -> tuple[list[list[float]], list[list[float]]]:
        """
        Parse output from T2CatalogMatch.

        Made complicated as returns can be both single and lists.


        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]
        group_dist: list[list[float]] = [[], [], [], [], [], [], []]

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

                if cat_name == "NEDz_extcats":
                    # at some point: verify whether 0.03 was the NEDz_extcats cut.
                    if cat_match["dist2transient"] < 2 and cat_match["z"] < 0.03:
                        group_z[0].append(cat_match["z"])
                        group_dist[0].append(cat_match["dist2transient"])
                    elif cat_match["dist2transient"] < 20 and cat_match["z"] < 0.05:
                        group_z[2].append(cat_match["z"])
                        group_dist[2].append(cat_match["dist2transient"])
                    else:
                        group_z[3].append(cat_match["z"])
                        group_dist[3].append(cat_match["dist2transient"])

                    # Implicit restriction as tests where done with this max matching radius

                if cat_name == "SDSS_spec" and cat_match["dist2transient"] < 10:
                    group_z[1].append(cat_match["z"])
                    group_dist[1].append(cat_match["dist2transient"])
                # Implicit restriction as tests where done with this max matching radius
                if (
                    cat_name == "GLADEv23"
                    and cat_match["dist2transient"] < 10
                    and cat_match["z"] is not None
                ):
                    if cat_match["z"] < 0.05:
                        group_z[2].append(cat_match["z"])
                        group_dist[2].append(cat_match["dist2transient"])
                    else:
                        group_z[3].append(cat_match["z"])
                        group_dist[3].append(cat_match["dist2transient"])

                if cat_name == "LSPhotoZZou":
                    # Spec
                    if cat_match["specz"] is not None and cat_match["specz"] > -0.1:
                        if (
                            cat_match["specz"] < 0.15
                            and cat_match["dist2transient"] < 10
                        ):
                            group_z[1].append(cat_match["specz"])
                            group_dist[1].append(cat_match["dist2transient"])
                        elif cat_match["specz"] < 0.2:
                            group_z[2].append(cat_match["specz"])
                            group_dist[2].append(cat_match["dist2transient"])
                        else:
                            group_z[4].append(cat_match["specz"])
                            group_dist[4].append(cat_match["dist2transient"])

                    # Photo-z
                    if cat_match["photoz"] is not None and cat_match["photoz"] > -0.1:
                        if cat_match["photoz"] < 0.1:
                            group_z[3].append(cat_match["photoz"])
                            group_dist[3].append(cat_match["dist2transient"])
                        elif cat_match["photoz"] < 0.2:
                            group_z[4].append(cat_match["photoz"])
                            group_dist[4].append(cat_match["dist2transient"])
                        elif cat_match["dist2transient"] < 20:
                            group_z[5].append(cat_match["photoz"])
                            group_dist[5].append(cat_match["dist2transient"])
                        else:
                            group_z[6].append(cat_match["photoz"])
                            group_dist[6].append(cat_match["dist2transient"])

                if cat_name == "wiseScosPhotoz" and (
                    cat_match["zPhoto_Corr"] is not None
                    and cat_match["zPhoto_Corr"] > -0.1
                ):
                    if cat_match["zPhoto_Corr"] < 0.2:
                        group_z[4].append(cat_match["zPhoto_Corr"])
                        group_dist[4].append(cat_match["dist2transient"])
                    else:
                        group_z[5].append(cat_match["zPhoto_Corr"])
                        group_dist[5].append(cat_match["dist2transient"])

                if cat_name == "twoMPZ":
                    # Photoz
                    if cat_match["zPhoto"] is not None and cat_match["zPhoto"] > -0.1:
                        if cat_match["zPhoto"] < 0.03:
                            group_z[2].append(cat_match["zPhoto"])
                            group_dist[2].append(cat_match["dist2transient"])
                        else:
                            group_z[3].append(cat_match["zPhoto"])
                            group_dist[3].append(cat_match["dist2transient"])
                    # Specz
                    if cat_match["zSpec"] is not None and cat_match["zSpec"] > -0.1:
                        group_z[1].append(cat_match["zSpec"])
                        group_dist[1].append(cat_match["dist2transient"])

                if cat_name == "PS1_photoz":
                    ps1_z_phot = float(cat_match["z_phot"]) / 1000
                    if cat_match["z_phot"] is not None and ps1_z_phot > -0.1:
                        if ps1_z_phot < 0.2:
                            group_z[3].append(ps1_z_phot)
                            group_dist[3].append(cat_match["dist2transient"])
                        elif ps1_z_phot < 0.4:
                            group_z[4].append(ps1_z_phot)
                            group_dist[4].append(cat_match["dist2transient"])
                        elif cat_match["dist2transient"] < 20:
                            group_z[5].append(ps1_z_phot)
                            group_dist[5].append(cat_match["dist2transient"])
                        else:
                            group_z[6].append(ps1_z_phot)
                            group_dist[6].append(cat_match["dist2transient"])

                # Implicit restriction as tests where done with this max matching radius
                if (
                    cat_name == "NEDz"
                    and cat_match["dist2transient"] < 10
                    and cat_match["z"] < 0.4
                ):
                    group_z[2].append(cat_match["z"])
                    group_dist[2].append(cat_match["dist2transient"])

                # Also check for manual override
                if self.catalogmatch_override:
                    for or_catname, or_catdict in self.catalogmatch_override.items():
                        if or_catname == cat_name:
                            try:
                                cat_z = float(cat_match[or_catdict["z_keyword"]])
                                if (
                                    float(cat_match["dist2transient"])
                                    < or_catdict["max_distance"]
                                    and cat_z < or_catdict["max_redshift"]
                                ):
                                    group_z[or_catdict["z_group"] - 1].append(cat_z)
                            except ValueError:
                                self.logger.info(
                                    "Cannot parse z", extra={"catdict": cat_match}
                                )

        return group_z, group_dist

    def _get_matchbts_groupz(self, t2_res: dict[str, Any]) -> list[list[float]]:
        """
        Parse output from T2MatchBTS.

        Any transient with a redshift with two decimals (from SN template matching) is put in Group II,
        those with more (from host, high-res spec) are put into Group I.

        Return:
        One list for each of the seven redshift cateogries
        """

        group_z: list[list[float]] = [[], [], [], [], [], [], []]

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
            else:
                group_z[1].append(z)

        self.logger.debug(f" bts match yield {group_z}")

        return group_z

    def get_ampelZ(self, t2_views: Sequence[T2DocView]) -> UBson | UnitResult:
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

        # Loop through t2_views and collect information.
        for t2_view in t2_views:
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            t2_res = cast(dict[str, Any], t2_view.get_payload())

            if t2_view.unit == "T2LSPhotoZTap":
                new_zs, new_dists = self._get_lsphotoz_groupz(t2_res)
            elif t2_view.unit in {"T2CatalogMatch", "T2CatalogMatchLocal"}:
                new_zs, new_dists = self._get_catalogmatch_groupz(t2_res)
            elif t2_view.unit == "T2MatchBTS":
                new_zs = self._get_matchbts_groupz(t2_res)
            else:
                self.logger.error(f"No instructions for dealing with {t2_view.unit}")
                return UnitResult(code=DocumentCode.T2_MISSING_INFO)

            for k in range(7):
                if len(new_zs[k]) > 0:
                    group_redshifts[k].extend(new_zs[k])
                if len(new_dists[k]) > 0:
                    group_distances[k].extend(new_dists[k])
            self.logger.debug(f"group_z after {t2_view.unit}: {group_redshifts}")

        # Check for best match
        t2_output: dict[str, UBson] = {
            "group_zs": group_redshifts,
            "group_dists": group_distances,
        }
        for k in range(7):
            if (k + 1) > self.max_redshift_category:
                # No matches with sufficient precision
                break
            if len(group_redshifts[k]) > 0:
                t2_output["ampel_z"] = float(np.mean(group_redshifts[k]))
                t2_output["group_z_precision"] = self.category_precision[k]
                t2_output["group_z_nbr"] = k + 1
                t2_output["ampel_dist"] = float(
                    np.mean(group_distances[k])
                )  # distance transient to matched z sources NOT TO EARTH
                # We then do *not* look for higher group (more uncertain) matches
                break
        if self.catalogmatch_override:
            t2_output["AmpelZ-Warning"] = "Override catalog in use."

        return t2_output

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
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )
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
        light_curve: LightCurve,
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """

        Parse t2_views from catalogs that were part of the redshift studies.
        Return these together with a "best estimate" - ampel_z

        """

        if not t2_views:  # Should not happen actually, T2Processor catches that case
            self.logger.error("Missing tied t2 views")
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        return self.get_ampelZ(t2_views)
