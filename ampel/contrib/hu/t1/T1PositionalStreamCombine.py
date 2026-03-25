#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t1/T1PositionalStreamCombine.py
# License:             BSD-3-Clause
# Author:              Jannis Necker <jannis.necker@gmail.com>
# Date:                24.03.2026
# Last Modified Date:  24.03.2026
# Last Modified By:    Jannis Necker <jannis.necker@gmail.com>


import numpy as np
import pymongo

from ampel.abstract.AbsT1CombineUnit import AbsT1CombineUnit
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.util.meanpos import mean_position
from ampel.core.ContextUnit import ContextUnit
from ampel.model.operator.AllOf import AllOf
from ampel.model.operator.AnyOf import AnyOf
from ampel.struct.T1CombineResult import T1CombineResult
from ampel.types import Tag

ARCSEC_IN_RAD = np.pi / 180 / 3600
SQDEG_IN_SR = (np.pi / 180) ** 2


class T1PositionalStreamCombine(AbsT1CombineUnit, ContextUnit):
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
    primary_tag: AnyOf[Tag] | AllOf[Tag]

    # tag(s) of the secondary streams to be included
    secondary_tag: AnyOf[Tag] | AllOf[Tag]

    # MongoDB to save matches
    mongo_uri: str
    database_name: str
    collection_name: str = "associations"
    reset_db: bool = False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._primary_operator = all if isinstance(self.primary_tag, AllOf) else any
        self._secondary_operator = all if isinstance(self.secondary_tag, AllOf) else any
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
            col.create_index("primary_stock")
            col.create_index("secondary_stock")

        self._col = self._client[self.database_name][self.collection_name]

    def select_closest_source(
        self, datapoints: list[DataPoint], ra_ref: float, dec_ref: float
    ) -> tuple[list[DataPoint], tuple[float, float, float, float], float]:
        # sort the alerts per source name
        source_names = np.unique([dp["stock"] for dp in datapoints])
        sorted_dps_dict = {
            zn: [dp for dp in datapoints if dp["stock"] == zn] for zn in source_names
        }

        # loop over all sources and find the closest one
        min_distance = np.inf
        closest_source_name = None
        associated_mp = None
        for zn, ads in sorted_dps_dict.items():
            mp = mean_position(
                [dp["body"]["ra"] for dp in ads],
                [dp["body"]["dec"] for dp in ads],
            )
            distance = np.sqrt(
                (np.radians(mp[0]) - np.radians(ra_ref)) ** 2
                * np.cos(np.radians(dec_ref)) ** 2
                + (np.radians(mp[1]) - np.radians(dec_ref)) ** 2
            )
            if distance < min_distance:
                min_distance = distance
                closest_source_name = zn
                associated_mp = mp

        if closest_source_name is None:
            # We checked earlier already that there are some alerts so this should not happen
            raise ValueError("No source found!")

        # selecting the alerts for the closest source
        return (
            sorted(sorted_dps_dict[closest_source_name], key=lambda x: x["body"]["jd"]),
            associated_mp,
            min_distance,
        )

    def combine(self, datapoints_list: list[DataPoint]) -> T1CombineResult:
        """
        :param datapoints_list: list of datapoints to combine
        :return: tuple of UBson or UnitResult and StockId
        """

        primary_dps = [
            dp
            for dp in datapoints_list
            if self._primary_operator([t in self.primary_tag for t in dp["tag"]])
        ]
        primary_stock = primary_dps[0]["stock"]
        secondary_dps = [
            dp
            for dp in datapoints_list
            if self._primary_operator([t in self.primary_tag for t in dp["tag"]])
        ]
        secondary_stocks = np.unique([dp["stock"] for dp in secondary_dps])

        # check if primary stock has already a match
        prv_match = None
        col = self._col
        for d in col.find({"primary_stock": primary_stock}):
            if prv_match:
                raise RuntimeError(
                    f"Corrupt match database: more than one match found for {primary_stock}!"
                )
            prv_match = d

        # if there was a previous match, get associated data
        if prv_match:
            # make sure the previous match is in the selected data
            if prv_match["secondary_stock"] in secondary_stocks:
                posterior = prv_match["p_association"]
                secondary_stock = prv_match["secondary_stock"]
                dist = prv_match["dist"]

        # if no match proceed with calculation
        else:
            # calculate position of primary stream
            ra, dec, dra, ddec = mean_position(
                [dp["body"]["ra"] for dp in primary_dps],
                [dp["body"]["dec"] for dp in primary_dps],
            )

            # select closest source from secondary stream
            closest_dps, pos, dist = self.select_closest_source(secondary_dps, ra, dec)
            secondary_stock = closest_dps[0]["stock"]

            # calculate posterior association probability
            sigma1_sq = dra**2 + ddec**2
            sigma2_sq = pos[2] ** 2 + pos[3] ** 2
            sigma_sq_rad = (sigma1_sq + sigma2_sq) * SQDEG_IN_SR
            posterior = 1 / (
                self._rho1 * sigma_sq_rad / 2 * np.exp(dist**2 / (2 * sigma_sq_rad)) + 1
            )

        # check if secondary stock is already better associated to another source
        prv_match_secondary = None
        for d in col.find({"secondary_stock": secondary_stock}):
            if prv_match_secondary:
                raise RuntimeError(
                    f"Corrupt match database: more than one match found for {secondary_stock}!"
                )
            prv_match_secondary = d

        # if the association is better, the primary source has no match
        if prv_match_secondary and (prv_match_secondary["p_association"] > posterior):
            # update database accordingly
            if prv_match:
                col.update_one(
                    {"primary_stock": primary_stock},
                    {"$set": {"superseded": prv_match_secondary["secondary_stock"]}},
                )

            # combine only data from primary stream
            return T1CombineResult(dps=[dp["id"] for dp in primary_dps])

        # note the association in the database
        col.insert_one(
            {
                "primary_stock": primary_stock,
                "secondary_stock": secondary_stock,
                "dist": dist,
                "posterior": posterior,
            }
        )

        # combine data from primary stream and secondary source
        return T1CombineResult(
            dps=[dp["id"] for dp in primary_dps + closest_dps],
            meta={"stock": closest_dps[0]["stock"], "p_association": posterior},
        )
