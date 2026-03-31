#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t1/T1PositionalStreamCombine.py
# License:             BSD-3-Clause
# Author:              Jannis Necker <jannis.necker@gmail.com>
# Date:                24.03.2026
# Last Modified Date:  24.03.2026
# Last Modified By:    Jannis Necker <jannis.necker@gmail.com>
from collections.abc import Callable

import numpy as np

from ampel.abstract.AbsT1CombineUnit import AbsT1CombineUnit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.util.meanpos import mean_position
from ampel.model.operator.AllOf import AllOf
from ampel.model.operator.AnyOf import AnyOf
from ampel.struct.T1CombineResult import T1CombineResult
from ampel.types import Tag

ARCSEC_IN_RAD = np.pi / 180 / 3600
SQDEG_IN_SR = (np.pi / 180) ** 2
FWHM_TO_STD = 1 / (2 * np.sqrt(2 * np.log(2)))


class T1PositionalStreamCombine(AbsT1CombineUnit):
    """
    Combine based on the position of alerts.
    Based on https://www.overleaf.com/read/hpbjsjrrxpym#b8d122
    """

    # the minimum posterior probability for any match
    min_posterior: float = 0.9

    # The density of sources in the primary alert stream in
    # 1/deg^2 used as the prior.
    # In the future this could be replaced with a
    # general PriorModel that can be evaluated per sky position.
    nu1: float

    # the resolutions of the two alert streams in arcseconds
    sigma1: float
    sigma2: float

    # tag(s) of the primary stream
    primary_tag: AnyOf[Tag] | AllOf[Tag] | Tag

    # tag(s) of the secondary streams to be included
    secondary_tag: AnyOf[Tag] | AllOf[Tag] | Tag

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # compile primary tag and operator
        self._primary_operator, self._primary_tag = self._compile_operator(
            self.primary_tag
        )
        self._secondary_operator, self._secondary_tag = self._compile_operator(
            self.secondary_tag
        )
        self._rho1 = 4 * np.pi * self.nu1 / SQDEG_IN_SR

        # identify the prior in the database
        # for now it can just be the prior value
        self._prior_hash = self.nu1

    @staticmethod
    def _compile_operator(
        tag_spec: AnyOf[Tag] | AllOf[Tag] | Tag,
    ) -> tuple[Callable, set[Tag]]:
        if isinstance(tag_spec, AnyOf):
            operator = any
            tags = set(tag_spec.any_of)
        elif isinstance(tag_spec, AllOf):
            operator = all
            tags = set(tag_spec.all_of)
        else:
            operator = any
            tags = {tag_spec}
        return operator, tags

    def combine(self, datapoints_list: list[DataPoint]) -> T1CombineResult:
        """
        :param datapoints_list: list of datapoints to combine
        :return: tuple of UBson or UnitResult and StockId
        """

        primary_dps = [
            dp
            for dp in datapoints_list
            if self._primary_operator([t in dp["tag"] for t in self._primary_tag])
        ]
        primary_stock = primary_dps[0]["stock"]
        secondary_dps = [
            dp
            for dp in datapoints_list
            if self._secondary_operator([t in dp["tag"] for t in self._secondary_tag])
        ]

        # calculate position of primary stream
        ra, dec, dra, ddec = mean_position(
            [dp["body"]["ra"] for dp in primary_dps if "ra" in dp["body"]],
            [dp["body"]["dec"] for dp in primary_dps if "dec" in dp["body"]],
        )

        # sort the alerts per source name
        source_names = np.unique([dp["stock"] for dp in secondary_dps])
        sorted_dps_dict = {
            zn: [dp for dp in secondary_dps if dp["stock"] == zn] for zn in source_names
        }

        # loop over all sources and find the closest one
        matching_result = np.zeros(
            len(sorted_dps_dict), dtype=[("zn", "int"), ("p", "<f4"), ("dist", "<f4")]
        )  # zn, posterior, distance
        for i, (zn, ads) in enumerate(sorted_dps_dict.items()):
            mp = mean_position(
                [dp["body"]["ra"] for dp in ads],
                [dp["body"]["dec"] for dp in ads],
            )

            # fallback for a single datapoint
            if len(ads) == 1:
                std_from_fwhm = ads[0]["body"]["fwhm"] * FWHM_TO_STD
                mp = mp[0], mp[1], std_from_fwhm, std_from_fwhm

            distance = np.sqrt(
                (np.radians(mp[0]) - np.radians(ra)) ** 2 * np.cos(np.radians(dec)) ** 2
                + (np.radians(mp[1]) - np.radians(dec)) ** 2
            )

            # calculate posterior association probability
            sigma1_sq = dra**2 + ddec**2
            sigma2_sq = mp[2] ** 2 + mp[3] ** 2
            sigma_sq_rad = (sigma1_sq + sigma2_sq) * SQDEG_IN_SR
            posterior = 1 / (
                self._rho1 * sigma_sq_rad / 2 * np.exp(distance**2 / (2 * sigma_sq_rad))
                + 1
            )

            matching_result[i] = (zn, posterior, distance / ARCSEC_IN_RAD)

        posterior_beyond_threshold = matching_result["p"] > self.min_posterior

        # if there is no good match combine only datapoints from primary stream
        if sum(posterior_beyond_threshold) == 0:
            return T1CombineResult(dps=[dp["id"] for dp in primary_dps])

        # select best matching source from secondary stream
        best_match_index = np.argmax(matching_result["p"])
        best_match = matching_result[best_match_index]

        # note the association in the database
        selected_dps = [dp["id"] for dp in primary_dps + sorted_dps_dict[best_match[0]]]
        body = {
            "primary_stock": primary_stock,
            "best_match": best_match.tolist(),
            "other_matches": np.delete(
                matching_result, best_match_index, axis=0
            ).tolist(),
        }

        # combine data from primary stream and secondary source
        return T1CombineResult(
            dps=selected_dps,
            body=body,
        )
