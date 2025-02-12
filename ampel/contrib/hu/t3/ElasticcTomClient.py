#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ElasticcTomClient.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                11.04.2022
# Last Modified Date:  11.04.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from collections.abc import Sequence
from typing import Any

import backoff
import requests
from requests import HTTPError
from typing_extensions import TypedDict


class ClassificationDict(TypedDict):
    classifierName: str
    classifierParams: str
    classId: int
    probability: float


class ElasticcClassification(TypedDict):
    alertId: int
    diaSourceId: int
    elasticcPublishTimestamp: int
    brokerIngestTimestamp: int
    brokerName: str
    brokerVersion: str
    classifications: Sequence[ClassificationDict]


class Elasticc2ClassificationDict(TypedDict):
    classId: int
    probability: float


class Elasticc2Report(TypedDict):
    alertId: int
    diaSourceId: int
    elasticcPublishTimestamp: int
    brokerIngestTimestamp: int
    brokerName: str
    brokerVersion: str
    classifierName: str
    classifierParams: str
    classifications: Sequence[Elasticc2ClassificationDict]


class ElasticcTomClient:
    """
    Initiate a session for reporting ELEAsTICC classification records to the
    DESC TOM system.

    Requires a desc username and password. Each call to tom_post
    attempts to put a ElasticcClassification dict to the tom.

    A successful put returns a 'dbMessageIndex' value.

    todo: do we need to robustify also the step where the session is created?
    (is the session best started in a post_init()?)

    """

    def __init__(
        self,
        desc_username: str,
        desc_password: str,
        logger,
        tom_url: str = "https://desc-tom.lbl.gov",
        endpoint: str = "/elasticc/brokermessage/",
    ):
        self.logger = logger

        self.tom_url = tom_url
        self.endpoint = endpoint

        # Setup django connection. From Rob Knop:
        # There's a bit of a dance to log in since django
        # requires a csrftoken not only in its headers, but
        # also in the POST data to log in; do a quick GET
        # to the login URI to get that token.  (There must
        # be a cleaner way.)
        self.session = requests.session()
        self.session.get(f"{self.tom_url}/accounts/login/")
        self.session.post(
            f"{self.tom_url}/accounts/login/",
            data={
                "username": desc_username,
                "password": desc_password,
                "csrfmiddlewaretoken": self.session.cookies["csrftoken"],
            },
        )
        self.csrfheader = {"X-CSRFToken": self.session.cookies["csrftoken"]}

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
    def tom_post(
        self,
        classification: Any,
    ) -> dict[Any, Any]:
        response = self.session.put(
            f"{self.tom_url}{self.endpoint}",
            json=classification,
            headers=self.csrfheader,
        )

        if response.ok:
            self.logger.debug(
                "ElasticcTomClient submit done", extra={"payload": classification}
            )
            return {"success": True, **response.json()}

        self.logger.info(
            "ElasticcTomClient submit fail", extra={"payload": classification}
        )
        return {
            "success": False,
            "response": response.status_code,
            "response_body": response.text,
        }
