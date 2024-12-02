#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ScoreTNSObjects
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                8.1.2024
# Last Modified Date:  8.1.2024
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from typing import Any

from astropy.time import Time

from ampel.contrib.hu.t3.AbsScoreCalculator import AbsScoreCalculator
from ampel.contrib.hu.t3.tns.TNSClient import TNSClient
from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.secret.NamedSecret import NamedSecret


class ScoreTNSObjects(AbsScoreCalculator):
    """

    Calculate score based on detection time reported to TNS, if any.

    """

    # Bot api key frm TNS
    tns_api_key: NamedSecret[str] = NamedSecret[str](label="tns/api/token")
    tns_id: int = 49023
    tns_name: str = "ZTF_AMPEL_NEW"

    # Identifying the SN, done through coordinate matching
    maxdist: float = 2.0  # max squared dist, in arcsec.
    # TNS prefix has to be amogn these to score
    tns_prefix: list[str] = [
        "SN"
    ]  # use ['SN,'AT'] to also include unclassified transients

    # Calculating the score
    powerscale: float = -2.0  # timediff to tzero to this power gives scaling.
    mindiff: float = 0.1  # lower limit for timediffs, to avoid zeros for ZTF detections + some TNS report incosistencies

    t2scores_from: list[str] = ["T2InfantCatalogEval"]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.client = TNSClient(
            TNSToken(
                id=self.tns_id, name=self.tns_name, api_key=self.tns_api_key.get()
            ),
            5,  # timeout
            1,  # parallel
            self.logger,
        )

    def get_tns_discovery(self, ra: float, dec: float) -> None | float:
        tdisc = None
        for doc in self.client.search(
            ra=ra,
            dec=dec,
            radius=self.maxdist,
            units="arcsec",
        ):
            self.logger.debug("got from tns", extra={"doc": doc})
            if doc["name_prefix"] not in self.tns_prefix:
                continue
            if (discdate := doc.get("discoverydate", None)) is None:
                continue
            tdisc = Time(discdate, format="iso", scale="utc").jd
            self.logger.debug("found a discovery date", extra={"tdisc": tdisc})
        return tdisc

    def evaluate(self, t2unit: str, t2_result: dict[str, Any]) -> float:
        """
        Check whether this transient was submitted to TNS and, if so, add a
        score depending on how close to the reported detection time we got.
        """

        # Assuming that we have info from T2InfantCatalogEval, could of course grab coordinates and time from other sources.
        if t2unit != "T2InfantCatalogEval":
            return 0

        if (ra := t2_result.get("ra")) and (dec := t2_result.get("dec")):
            tjd = self.get_tns_discovery(ra, dec)
        else:
            tjd = None

        # Object not reported to TNS - assuming this was not an interesting transient!
        if tjd is None:
            return 0

        # calculate score
        tdiff = max(t2_result["t_max"] - tjd, self.mindiff)
        return (tdiff) ** self.powerscale
