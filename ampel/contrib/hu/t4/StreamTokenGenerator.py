#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-ZTF/ampel/ztf/t4/resource/T3ZTFArchiveTokenGenerator.py
# License:             BSD-3-Clause
# Author:              Akshay Eranhalodi
# Date:                16.08.2024
# Last Modified Date:  16.08.2024
# Last Modified By:    akshay eranhalodi <firstname.lastname@desy.de>

import random
import time
from typing import Any, Literal

from astropy.time import Time
from requests_toolbelt.sessions import BaseUrlSession

from ampel.abstract.AbsT4Unit import AbsT4Unit
from ampel.secret.NamedSecret import NamedSecret


class StreamTokenGenerator(AbsT4Unit):
    """
    Stream based token generator for:
    1. Cone search
    2. Healpix search
    3. Epoch search
    """

    archive_token: NamedSecret[str] = NamedSecret(label="ztf/archive/token")

    #: Base URL of archive service
    archive: str = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/"
    resource_name: str = "ztf_stream_token"

    mode: Literal["cone", "healpix", "epoch"] = "cone"

    # supply ra, dec and radius if mode is "cone"
    cone: None | dict[str, float] = None
    candidate: None | dict[str, Any] = None

    nside: None | int = None  # use when mode is "healpix"
    pixels: None | list = None  # use when mode is "healpix"

    # default jd range: [ztf_start, jd_now]
    jd_range: tuple[float, float] = (2458195.0, float(Time.now().jd))

    #: seconds to wait for query to complete
    timeout: float = 60

    debug: bool = False
    query: dict[str, Any] | None = None

    def do(self) -> dict[str, str]:
        if self.candidate:
            candidate = self.candidate
        else:  # default candidate filter
            candidate = {
                "rb": {"$gt": 0.3},
                "ndethist": {"$gte": 2},
                "isdiffpos": {"$in": ["t", "1"]},
            }

        if self.mode == "cone":
            if not self.cone:
                raise TypeError(
                    "Missing required argument==> cone={ra:value, dec:value, radius:value}"
                )
            self.query = {"cone": self.cone, "candidate": candidate}

        elif self.mode == "healpix":
            if not (self.nside and self.pixels):
                raise TypeError('Missing required argument: "nside" and/or "pixels"')

            self.query = {
                "regions": [{"nside": self.nside, "pixels": self.pixels}],
                "jd": {"$gt": self.jd_range[0], "$lt": self.jd_range[1]},
                "candidate": candidate,
            }

        elif self.mode == "epoch":
            self.query = {
                "jd": {"$gt": self.jd_range[0], "$lt": self.jd_range[1]},
                "candidate": candidate,
            }

        session = BaseUrlSession(
            self.archive if self.archive.endswith("/") else self.archive + "/"
        )
        session.headers["authorization"] = f"bearer {self.archive_token.get()}"

        response = session.post("streams/from_query?programid=1", json=self.query)

        rd = response.json()

        try:
            token = rd.pop("resume_token")
        except KeyError as exc:
            raise ValueError(f"Unexpected response: {rd}") from exc

        # wait for query to finish
        t0 = time.time()
        delay = 1
        while time.time() - t0 < self.timeout:
            response = session.get(f"stream/{token}")
            if response.status_code not in (423, 404):
                break
            time.sleep(random.uniform(0, delay))
            delay *= 2
        else:
            raise RuntimeError(
                f"{response.url} still locked after {time.time() - t0:.0f} s"
            )
        response.raise_for_status()
        self.logger.info("Stream created", extra=response.json())

        return {self.resource_name: token}
