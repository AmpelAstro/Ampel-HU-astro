#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/AstroColibriClient.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                13.11.2022
# Last Modified Date:  13.11.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

# FIXME: restore mypy when this is actually ready
# type: ignore

import base64
import json
from typing import Any

import backoff
import requests
from requests import HTTPError
from requests.auth import HTTPBasicAuth

# from datetime import datetime
# from astropy.time import Time

# url = 'https://astro-colibri.herokuapp.com/add_ampel_transient'


class AstroColibriPhot:
    mag: float
    band: str
    datetime: str  # e.g. datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")


class AstroColibriPost:
    timestamp: str  # Timestamp of first detection
    trigger_id: str  # AstroColibri identifier -> should be TNS ID
    source_name: str  #  SN|AT[space]2022ajn
    ra: float
    dec: float
    err: float
    type: str = "ot_sn"  # 'ot_sn'
    observatory: str | None = "ztf"
    ampel_attributes: list[str] = []
    photometry: dict[str, AstroColibriPhot]
    discovery_name: str  #  ZTFxxx
    # dev
    # ampel_id: str             # Internal AMPEL ID. To be added
    # ivorn: str | None = None  # If existing, this can be added.


class AstroColibriClient:
    """
    Initiate a session for communicating with an AstroColibri instance.

    """

    # url (currently test version)
    # colibri_url = 'https://astro-colibri.herokuapp.com/add_ampel_transient'
    api_url = "https://astro-colibri.science"

    def __init__(self, colibri_username: str, colibri_password: str, logger):
        self.logger = logger

        self.session = requests.session()
        self.session.auth = HTTPBasicAuth(colibri_username, colibri_password)
        # or ?
        # self.auth=HTTPBasicAuth(colibri_username, colibri_password)

    def store_image(self, image_file: str):
        with open(image_file, "rb") as f:
            im_bytes = f.read()
        im_b64 = base64.b64encode(im_bytes).decode("utf8")
        headers = {"Content-type": "application/json", "Accept": "text/plain"}
        payload = json.dumps(
            {"image": im_b64, "file_name": image_file, "other_key": "value"}
        )
        #        print('image payload', payload)
        #        response = requests.post(self.api_url+'/add_image_to_storage',
        #                            data=payload, headers=headers, auth=self.session.auth)
        response = self.session.post(
            self.api_url + "/add_image_to_storage", data=payload, headers=headers
        )
        return response.json()["url_lc"]

    # robustify post
    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: not isinstance(e, HTTPError)
        or e.response.status_code not in {503, 504, 429, 408},
        max_time=60,
    )
    def firestore_post(
        self, data: AstroColibriPost, image_path: str | None = None
    ) -> dict[Any, Any]:
        # Upload image file if provided
        if image_path is not None:
            lc_url = self.store_image(image_path)
            self.logger.info(
                "AstroClibri image upload", extra={"local": image_path, "url": lc_url}
            )
            data["lc_url"] = lc_url

        # TODO: not working
        response = self.session.post(self.api_url + "/add_ampel_transient", json=data)

        if response.ok:
            self.logger.debug(
                "AstroColibriClient submit success", extra={"payload": data}
            )
            return {"success": True, **response.json()}

        self.logger.info("AstroColibriClient submit fail", extra={"payload": data})
        return {"success": False, "response": response.status_code}
