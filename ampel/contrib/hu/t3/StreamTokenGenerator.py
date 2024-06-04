#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-ZTF/ampel/ztf/t3/resource/T3ZTFArchiveTokenGenerator.py
# License:             BSD-3-Clause
# Author:              Akshay Eranhalodi
# Date:                19.12.2023
# Last Modified Date:  19.12.2023
# Last Modified By:    akshay eranhalodi <firstname.lastname@desy.de>

import random
import time
from typing import Any, Literal

from requests_toolbelt.sessions import BaseUrlSession

from ampel.abstract.AbsT3PlainUnit import AbsT3PlainUnit
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.Resource import Resource
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


class StreamTokenGenerator(AbsT3PlainUnit):
    archive_token: NamedSecret[str] = NamedSecret(label="ztf/archive/token")

    #: Base URL of archive service
    archive: str = "https://ampel.zeuthen.desy.de/api/ztf/archive/v3/"
    resource_name: str = "ztf_stream_token"
    mode: Literal["cone", "healpix"] = "healpix"

    cone: None | dict[str, float] = None
    candidate: None | dict[str, Any] = None

    nside: None | int = None
    pixels: None | list = None
    jd_range: list = [0, 9000000]

    #: seconds to wait for query to complete
    timeout: float = 60

    debug: bool = False

    def process(self, t3s: T3Store) -> UBson | UnitResult:
        if self.candidate:
            candidate = self.candidate
        else:
            candidate = {
                "rb": {"$gt": 0.3},
                "ndethist": {"$gte": 2},
                "isdiffpos": {"$in": ["t", "1"]},
            }

        if self.mode == "cone":
            if not (self.cone):
                raise TypeError("Missing required argument: cone")
            query = {"cone": self.cone, "candidate": candidate}

        elif self.mode == "healpix":
            if not (self.nside and self.pixels):
                raise TypeError("Missing required argument: nside and/or pixels")

            query = {
                "regions": [{"nside": self.nside, "pixels": self.pixels}],
                "jd": {"$gt": self.jd_range[0], "$lt": self.jd_range[1]},
                "candidate": candidate,
            }
        else:
            raise ValueError("Invalid mode; must be cone OR healpix")

        session = BaseUrlSession(
            self.archive if self.archive.endswith("/") else self.archive + "/"
        )
        session.headers["authorization"] = f"bearer {self.archive_token.get()}"

        response = session.post("streams/from_query?programid=1", json=query)

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
            if response.status_code != 423:
                break
            time.sleep(random.uniform(0, delay))
            delay *= 2
        else:
            raise RuntimeError(
                f"{session.base_url}stream/{token} still locked after {time.time() - t0:.0f} s"
            )
        response.raise_for_status()
        self.logger.info("Stream created", extra=response.json())

        r = Resource(name=self.resource_name, value=token)
        t3s.add_resource(r)

        if self.debug:
            return r.dict()

        return None
