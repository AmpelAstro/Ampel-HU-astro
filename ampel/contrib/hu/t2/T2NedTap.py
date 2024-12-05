#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2NedTap.py
# License:             BSD-3-Clause
# Author:              valery brinnel <firstname.lastname@gmail.com>
# Date:                09.03.2021
# Last Modified Date:  24.11.2021
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

import json
from math import acos, cos, pi, sin
from typing import Any

import requests

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson

RAD = pi / 180


class T2NedTap(AbsPointT2Unit):
    """
    See also:
    https://ned.ipac.caltech.edu/tap/sync?QUERY=SELECT+*+FROM+TAP_SCHEMA.tables&REQUEST=doQuery&LANG=ADQL&FORMAT=text
    Export all NED:
    https://ned.ipac.caltech.edu/tap/sync?QUERY=SELECT+*+FROM+NEDTAP.objdir&REQUEST=doQuery&LANG=ADQL&FORMAT=text
    """

    max_res: int = 5
    radius: float = 1 / 3600 * 20
    z_cut: float = 0.15

    request_timeout: int = 5
    verbose: bool = False

    # Example:
    # https://ned.ipac.caltech.edu/tap/sync?query=SELECT+TOP+5+prefname,pretype,ra,dec,z,zunc,zflag,n_spectra+FROM+objdir+WHERE+CONTAINS(POINT(%27J2000%27,ra,dec),CIRCLE(%27J2000%27,141.0678871,49.2484661,0.01))=1&LANG=ADQL&REQUEST=doQuery&FORMAT=json
    query: str = "https://ned.ipac.caltech.edu/tap/sync?query=SELECT+TOP+%i+prefname,pretype,ra,dec,z,zunc,zflag,n_spectra+FROM+objdir+WHERE+CONTAINS(POINT('J2000',ra,dec),CIRCLE('J2000',%s,%s,%s))=1+and+z<%s&LANG=ADQL&REQUEST=doQuery&FORMAT=json"

    #: Path to file created by mongoexport. No query will be performed:
    #: ned results will be imported using the mongo export file.
    #: Export command example: mongoexport -q '{"unit": "T2NedTap"}' -d Dipole_data -c t2 > /home/ampel/ned.export
    mongo_data: None | str

    #: Path to file created by mongoexport. No query will be performed:
    #: mongoexport -q '{"query": {"$exists": true}, "radius": {"$exists": true}}' -d Dipole_ext -c confid > /home/ampel/ned.confid
    mongo_confid: None | str

    do_query_if_missing: bool = False
    do_query_if_no_match: bool = False

    def process(self, datapoint: DataPoint) -> UBson | UnitResult:
        try:
            ra = datapoint["body"]["ra"]
            dec = datapoint["body"]["dec"]
        except KeyError:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        if getattr(self.logger, "verbose", 0) > 1:
            self.logger.debug(f"Transient position (ra, dec): {ra:.4f}, {dec:.4f} deg")

        if self.mongo_data and self.mongo_confid:
            e = self._load_mongo_export(
                self.mongo_data, f"stock:\"{datapoint['stock']}"
            )
            if e is None:
                if not self.do_query_if_missing:
                    return UnitResult(tag=["NO_QUERY", "EXT_NOT_FOUND"], code=2)
            else:
                if (
                    (c := self._load_mongo_export(self.mongo_confid, e["config"]))
                    and c.get("radius") == self.radius
                    and c.get("z_cut") == self.z_cut
                    and e.get("code") == 0
                    and (data := e.get("body", [{}])[0].get("data"))
                ):
                    return UnitResult(code=0, body=e["body"][0], tag=e["tag"])

                if not self.do_query_if_no_match:
                    return UnitResult(tag=["NO_QUERY", "EXT_NO_MATCH"], code=3)

        try:
            if getattr(self.logger, "verbose", 0) > 1:
                self.logger.info(
                    self.query % (self.max_res, ra, dec, self.radius, self.z_cut)
                )

            resp = requests.get(
                self.query % (self.max_res, ra, dec, self.radius, self.z_cut),
                timeout=self.request_timeout,
            )

            if resp.status_code != 200:
                self.logger.error(f"NED request response code: {resp.status_code}")
                return UnitResult(code=DocumentCode.RERUN_REQUESTED)

            r: dict[str, Any] = resp.json()

        except Exception as e:
            self.logger.error(
                "Connection error while sending request to NED", exc_info=e
            )
            # Service temporarily? unavail
            return UnitResult(code=DocumentCode.RERUN_REQUESTED)

        # ex: {0: 'prefname', 1: 'pretype', 2: 'ra', ...
        d: dict[int, str] = {i: v["name"] for i, v in enumerate(r["metadata"])}

        # Example:
        # [{
        #  'prefname': 'NGC 2856',
        #  'pretype': 'G',
        #  'ra': 141.0667033142,
        #  'dec': 49.249191788100006,
        #  'z': 0.00879899971,
        #  'zunc': 6.29999995e-05,
        #  'zflag': '',
        #  'n_spectra': 2
        # },
        # {'prefname': 'SDSS J092415.46+491506.9',
        #  'pretype': 'PofG',
        #  'ra': 141.06445138549998,
        #  'dec': 49.251921356000004,
        #  'z': 0.009257,
        #  'zunc': 5.7e-05,
        #  'zflag': None,
        #  'n_spectra': 1
        # }]
        reshaped: list[dict[str, Any]] = [
            {d[i]: v.strip() if isinstance(v, str) else v for i, v in enumerate(el)}
            for el in r["data"]
        ]

        unsorted = [el for el in reshaped if el["z"]]

        if not unsorted:
            self.logger.info("No catalog match")
            return UnitResult(tag="NED_NO_MATCH", code=1)

        # Compute separation
        for el in unsorted:
            el["sep"] = (
                acos(
                    sin(dec * RAD) * sin(el["dec"] * RAD)
                    + cos(dec * RAD) * cos(el["dec"] * RAD) * cos((ra - el["ra"]) * RAD)
                )
                * 206264.8062
            )  # to arcsecs

        data = list(sorted(unsorted, key=lambda k: k["sep"]))
        tags = ["NED_MATCH"]

        if len(data) > 1:
            tags.append("NED_MULTI_MATCH")

        if self.is_spec(data[0]):
            tags.append("NED_NEAREST_IS_SPEC")
        elif self.is_unsafe_spec(data[0]):
            tags.append("NED_NEAREST_IS_UNSAFE_SPEC")
        else:
            tags.append("NED_NEAREST_NOT_SPEC")

        if [el for el in data if self.is_spec(el)]:
            tags.append("NED_HAS_SPEC")

        if [el for el in data if self.is_unsafe_spec(el)]:
            tags.append("NED_HAS_UNSAFE_SPEC")

        if "NED_HAS_SPEC" not in tags and "NED_HAS_UNSAFE_SPEC" not in tags:
            tags.append("NED_NO_SPEC")

        return UnitResult(tag=tags, body={"data": data})

    def _load_mongo_export(self, fpath, match) -> None | dict:
        with open(fpath) as f:
            for l in f:
                if match in l:
                    return json.loads(l)
        return None

    def is_spec(self, d: dict) -> bool:
        return d.get("n_spectra", 0) > 0 or d["zflag"] == "SPEC"

    def is_unsafe_spec(self, d: dict) -> bool:
        return d["zflag"] is None and d.get("zunc", 0) < 0.0001 and len(str(d["z"])) > 9
