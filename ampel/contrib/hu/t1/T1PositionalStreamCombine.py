#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t1/T1PositionalStreamCombine.py
# License:             BSD-3-Clause
# Author:              Jannis Necker <jannis.necker@gmail.com>
# Date:                24.03.2026
# Last Modified Date:  24.03.2026
# Last Modified By:    Jannis Necker <jannis.necker@gmail.com>
import json
from collections.abc import Callable
from hashlib import md5

import numpy as np
import pymongo
from pymongo.errors import DuplicateKeyError

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

    # MongoDB to save matches
    mongo_uri: str
    database_name: str
    collection_name: str = "associations"
    reset_db: bool = False

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

        # set up Mongo collection to store matches
        self._client = pymongo.MongoClient(self.mongo_uri)
        if (self.database_name in self._client.list_database_names()) and self.reset_db:
            self._client.drop_database(self.database_name)

        # if database does not exist or was just dropped because self.reset_db
        if self.database_name not in self._client.list_database_names():
            col = self._client[self.database_name][self.collection_name]
            col.create_index(
                "secondary_stock", partialFilterExpression={"superseded": False}
            )
            col.create_index(
                "primary_stock", partialFilterExpression={"superseded": False}
            )
            col.create_index("match_id", unique=True)

        self._col = self._client[self.database_name][self.collection_name]

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
        posteriors = []
        for zn, ads in sorted_dps_dict.items():
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

            posteriors.append((zn, posterior, distance))

        posteriors = np.array(posteriors)
        posterior_beyond_threshold = posteriors[:, 1] > self.min_posterior
        no_match_res = T1CombineResult(dps=[dp["id"] for dp in primary_dps])

        # if there is no good match combine only datapoints from primary stream
        if sum(posterior_beyond_threshold) == 0:
            return no_match_res

        # select best matching source from secondary stream
        best_match = posteriors[np.argmax(posteriors[:, 1])]
        secondary_stock = best_match[0]

        # check if secondary stock is already better associated to another source
        prv_match_secondary = None
        col = self._col
        for d in col.find({"secondary_stock": secondary_stock, "superseded": False}):
            if prv_match_secondary:
                raise RuntimeError(
                    f"Corrupt match database: more than one match found for {secondary_stock}!"
                )
            prv_match_secondary = d

        # if the association is better, the primary source has no match
        if (
            prv_match_secondary
            and (prv_match_secondary["primary_stock"] != primary_stock)
            and (prv_match_secondary["p_association"] > best_match[1])
        ):
            return no_match_res

        # note the association in the database
        selected_secondary_dps = sorted_dps_dict[best_match[0]]
        selected_dps = [dp["id"] for dp in primary_dps + selected_secondary_dps]
        match_id = md5(
            json.dumps(
                [self._prior_hash, *selected_dps], separators=(",", ":")
            ).encode()
        ).hexdigest()
        body = {
            "primary_stock": primary_stock,
            "secondary_stock": secondary_stock,
            "dist": best_match[2],
            "posterior": best_match[1],
            "superseded": False,
            "superseded_by": None,
            "match_id": match_id,
        }
        try:
            col.insert_one(body)

        # case the match already exists in the database, e.g. for a re-run
        except DuplicateKeyError:
            d = col.find_one({"match_id": match_id})
            for k, v in d.items():
                assert body[k] == v, (
                    f"Match {match_id} already exists in database but with different values! {k}: {body[k]} vs {v}"
                )

        # supersede previous association
        col.update_many(
            {"primary_stock": primary_stock, "superseded": False},
            {"$set": {"superseded": True, "superseded_by": match_id}},
        )

        # combine data from primary stream and secondary source
        return T1CombineResult(
            dps=selected_dps,
            # this is no meta info and should better be stored in a body!
            meta={
                "stock": selected_secondary_dps[0]["stock"],
                "p_association": best_match[1],
            },
        )
