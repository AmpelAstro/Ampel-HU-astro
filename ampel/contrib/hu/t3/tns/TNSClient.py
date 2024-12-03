#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/aiotns
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                13.12.2018
# Last Modified Date:  04.11.2019
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import json
import sys
import time
import traceback
from collections.abc import Generator
from typing import Any

import backoff
import requests
from requests_toolbelt import MultipartEncoder
from requests_toolbelt.sessions import BaseUrlSession

from ampel.protocol.LoggerProtocol import LoggerProtocol

from .TNSToken import TNSToken


class TNSSession:
    def __init__(self, base_url: str, token: TNSToken) -> None:
        self._session = BaseUrlSession(base_url)
        self._session.headers["User-Agent"] = "tns_marker" + json.dumps(
            {"tns_id": token.id, "name": token.name, "type": "bot"}
        )
        self._token = token

    def __call__(
        self,
        method: str,
        data: dict | int,
        payload_label: str = "data",
        max_retries: int = 10,
    ) -> dict[str, Any]:
        """
        post to TNS
        """

        m = MultipartEncoder(
            fields={
                "api_key": (
                    None,
                    self._token.api_key.encode(),
                    "text/plain; charset=utf-8",
                ),
                payload_label: (None, json.dumps(data), "application/json"),
            }
        )

        for _ in range(max_retries):
            resp = self._session.post(
                f"https://www.wis-tns.org/api/{method}",
                data=m,
                headers={"content-type": m.content_type},
            )
            if resp.status_code == 429:
                wait = max(
                    (
                        int(resp.headers.get("x-cone-rate-limit-reset", 0))
                        if method.endswith("search")
                        else 0,
                        int(resp.headers.get("x-rate-limit-reset", 0)),
                        1,
                    )
                )
                time.sleep(wait)
                continue
            else:  # noqa: RET507
                break
        resp.raise_for_status()
        return resp.json()


class TNSClient:
    def __init__(
        self,
        token: TNSToken,
        timeout: float,
        maxParallelRequests: int,
        logger: LoggerProtocol,
    ) -> None:
        self.token = token
        self.logger = logger

        session = TNSSession("https://www.wis-tns.org/api/", self.token)
        # robustify tns_post
        self.request = backoff.on_exception(
            backoff.expo,
            requests.exceptions.RequestException,
            logger=None,
            giveup=self.is_permanent_error,
            on_giveup=self._on_giveup,
            on_backoff=self._on_backoff,
            max_time=timeout,
        )(session.__call__)

    def _on_backoff(self, details):
        exc_typ, exc, _ = sys.exc_info()
        err = (
            traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
            if exc
            else None
        )
        self.logger.warn(
            "backoff",
            extra={"exc": err, "wait": details["wait"], "tries": details["tries"]},
        )

    def _on_giveup(self, details):
        exc_typ, exc, _ = sys.exc_info()
        err = (
            traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
            if exc
            else None
        )
        self.logger.warn("gave up", extra={"exc": err, "tries": details["tries"]})

    @staticmethod
    def is_permanent_error(exc: Exception) -> bool:
        if isinstance(exc, requests.exceptions.HTTPError):
            return exc.response.status_code not in {500, 429, 404}
        return False

    def search(
        self,
        exclude: None | set[str] = None,
        **params: Any,
    ) -> Generator[dict[str, Any], None, None]:
        response = self.request("get/search", params)
        names = [item["objname"] for item in response["data"]["reply"]]
        self.logger.info("TNS search", extra={"params": params, "results": len(names)})
        yield from (
            self.request("get/object", {"objname": name})
            for name in names
            if not (exclude and name in exclude)
        )
        return

    def get(self, objname: str) -> dict[str, Any]:
        return self.request("get/object", {"objname": objname})

    def sendReport(self, report) -> None | int:
        response = self.request("bulk-report", report)
        if response["id_code"] == 200:
            return response["data"]["report_id"]
        return None

    def reportReply(self, report_id: int) -> Any:
        response = self.request("bulk-report-reply", report_id, "report_id")
        return response["data"]["feedback"]
