#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/t3/HealpixTokenGenerator.py
# License:             BSD-3-Clause
# Author:              jnordin
# Date:                24.03.2023
# Last Modified Date:  24.03.2023
# Last Modified By:    jnordin <jnordin@physik.hu-berlin.de>

import random
import time

from astropy.time import Time  # type: ignore
from requests_toolbelt.sessions import BaseUrlSession

from ampel.abstract.AbsT3PlainUnit import AbsT3PlainUnit
from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix, deres
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.Resource import Resource
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


class HealpixTokenGenerator(AbsT3PlainUnit):
    """
    Based on a URL to a Healpix map:
    - find pixels given requested prob contour.
    - request archive token for this stream.
    """

    # Process pixels with p-values lower than this limit
    pvalue_limit: float = 0.9

    # Name (signifier)
    map_name: str

    # URL for healpix retrieval
    map_url: str
    map_dir: str  # Local dir where map is saved. File with this name del

    archive_token: NamedSecret[str] = NamedSecret(label="ztf/archive/token")
    #: Base URL of archive service
    archive: str = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/"

    date_str: None | str = (
        None  # Start of time window we are interested in (default: event trigger time)
    )
    date_format: str = "iso"  # "%Y-%m-%d"
    delta_time: None | float = (
        None  # Length of time window in days (default: until now)
    )

    # overwrite candidate section in query
    candidate: dict = {}

    chunk_size: int | None = 500

    #: seconds to wait for query to complete
    timeout: float = 60

    debug: bool = False

    def process(self, t3s: T3Store) -> UBson | UnitResult:
        # Retrieve and process map
        ah = AmpelHealpix(
            map_name=self.map_name, map_url=self.map_url, save_dir=self.map_dir
        )
        map_hash = ah.process_map()

        # Get list of pixels within requested significance contour
        pixels = ah.get_pixelmask(self.pvalue_limit)
        self.logger.info(
            "",
            extra={
                "map": self.map_name,
                "hash": map_hash,
                "size": len(pixels),
                "nside": ah.nside,
            },
        )

        # JD time range
        if self.delta_time:
            if self.date_str:
                start_jd = Time(
                    self.date_str,
                    format=self.date_format,
                    scale="utc",
                ).jd
            else:
                start_jd = ah.trigger_time

            end_jd = start_jd + self.delta_time

        else:
            if self.date_str:
                start_jd = Time(
                    self.date_str,
                    format=self.date_format,
                    scale="utc",
                ).jd
            else:
                start_jd = ah.trigger_time
            end_jd = Time.now().jd

        session = BaseUrlSession(
            self.archive if self.archive.endswith("/") else self.archive + "/"
        )
        session.headers["authorization"] = f"bearer {self.archive_token.get()}"

        # Combine pixels when possible
        deresdict = deres(ah.nside, pixels)
        healpix_regions = [
            {"nside": nside, "pixels": members} for nside, members in deresdict.items()
        ]
        count = sum([len(region["pixels"]) for region in healpix_regions])

        hp_area = ah.get_maparea(self.pvalue_limit)

        candidate_dict = {
            "rb": {"$gt": 0.3},
            "magpsf": {"$gt": 15},
            "ndethist": {"$gt": 0, "$lte": 10},
            "jdstarthist": {"$gt": start_jd},
        }

        candidate_dict.update(self.candidate)

        # TODO: candidate optional input, jdstarthis = start_jd + epsilon
        query_dict = {
            "jd": {"$gt": start_jd, "$lt": end_jd},
            "regions": healpix_regions,
            "candidate": candidate_dict,
            "chunk_size": self.chunk_size,
        }

        # count alerts before trigger time without querying for them
        count_query_dict = {
            "jd": {"$gt": start_jd, "$lt": end_jd},
            "regions": healpix_regions,
        }
        endpoint_count = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/alerts/healpix/skymap/count"
        response_count = session.post(endpoint_count, json=count_query_dict)
        alert_count_nofilter = response_count.json()["count"]
        # print("ALERT COUNT NO FILTER", alert_count_nofilter)

        # print("HEALPIXTOKENGENERATOR::", query_dict)

        response = session.post(
            "streams/from_query",
            json=query_dict,
        )
        response.raise_for_status()

        rd = response.json()
        try:
            token = rd.pop("resume_token")
        except KeyError as exc:
            raise ValueError(f"Unexpected response: {rd}") from exc

        # wait for query to finish - is this needed, or handled by alert consumer?
        t0 = time.time()
        delay = 1
        while time.time() - t0 < self.timeout:
            response = session.get(f"stream/{token}")
            if response.status_code != 423:
                break
            time.sleep(random.uniform(0, delay))
            delay *= 2
        else:
            raise RuntimeError(
                f"{session.base_url}stream/{token} still locked after {time.time() - t0:.0f} s"
            )
        self.logger.info("Stream created", extra=response.json())

        # print("TOKENGENERATOR AAAAAAAAAAAAAAAAA::", response.json()["remaining"]["items"])
        if response.json().get("remaining"):
            queried_alerts = response.json()["remaining"]["items"]
        else:
            queried_alerts = 0
        # Package resource needed
        resource = {
            "map_name": self.map_name,
            "map_dir": self.map_dir,
            "map_url": self.map_url,
            "hash": map_hash,
            "token": token,
            "jd": ah.trigger_time,
            "map_area": hp_area,
            "alert_count_nofilter": alert_count_nofilter,
            "alert_count_query": queried_alerts,
        }

        r = Resource(name=self.map_name, value=resource)
        t3s.add_resource(r)
        r = Resource(name=self.map_name + "_token", value=token)
        t3s.add_resource(r)
        r = Resource(name="healpix_map_dir", value=self.map_dir)
        t3s.add_resource(r)
        r = Resource(name="healpix_map_hash", value=map_hash)
        t3s.add_resource(r)
        r = Resource(name="healpix_map_name", value=self.map_name)
        t3s.add_resource(r)
        r = Resource(name="map_area", value=hp_area)
        t3s.add_resource(r)
        r = Resource(name="alert_count_nofilter", value=alert_count_nofilter)
        t3s.add_resource(r)

        if self.debug:
            return r.dict()

        return None
