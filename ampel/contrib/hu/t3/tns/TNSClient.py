#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/aiotns
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                13.12.2018
# Last Modified Date:  04.11.2019
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import asyncio
import json
import sys
import traceback
from functools import partial

import aiohttp
import backoff
from aiohttp import ClientSession
from aiohttp.client_exceptions import (
    ClientConnectionError,
    ClientConnectorError,
    ClientResponseError,
    ServerDisconnectedError,
)

from .TNSToken import TNSToken


async def tns_post(
    session: ClientSession,
    semaphore: asyncio.Semaphore,
    method: str,
    token: TNSToken,
    data: dict,
    max_retries: int = 10,
) -> dict:
    """
    post to TNS, asynchronously
    """
    async with semaphore:
        for _ in range(max_retries):
            with aiohttp.MultipartWriter("form-data") as mpwriter:
                p: aiohttp.Payload = aiohttp.StringPayload(token.api_key)
                p.set_content_disposition("form-data", name="api_key")
                mpwriter.append(p)
                p = aiohttp.JsonPayload(data)
                p.set_content_disposition("form-data", name="data")
                mpwriter.append(p)
                resp = await session.post(
                    "https://www.wis-tns.org/api/" + method, data=mpwriter
                )
            if resp.status == 429:
                wait = max(
                    (
                        int(resp.headers.get("x-cone-rate-limit-reset", 0))
                        if method.endswith("search")
                        else 0,
                        int(resp.headers.get("x-rate-limit-reset", 0)),
                        1,
                    )
                )
                await asyncio.sleep(wait)
                continue
            else:  # noqa: RET507
                break
        resp.raise_for_status()
        return await resp.json()


class TNSClient:
    def __init__(self, token: TNSToken, timeout, maxParallelRequests, logger):
        self.token = token
        self.maxParallelRequests = maxParallelRequests
        self.logger = logger
        # robustify tns_post
        self.tns_post = backoff.on_exception(
            backoff.expo,
            (
                TimeoutError,
                ClientResponseError,
                ClientConnectorError,
                ClientConnectionError,
                ServerDisconnectedError,
            ),
            logger=None,
            giveup=self.is_permanent_error,
            on_giveup=self._on_giveup,
            on_backoff=self._on_backoff,
            max_time=timeout,
        )(tns_post)

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
    def is_permanent_error(exc):
        if isinstance(exc, ClientResponseError):
            return exc.code not in {500, 429}
        return False

    async def search(self, exclude=set(), **params):  # noqa: B006
        semaphore = asyncio.Semaphore(self.maxParallelRequests)
        async with ClientSession(
            headers={
                "User-Agent": "tns_marker"
                + json.dumps(
                    {"tns_id": self.token.id, "name": self.token.name, "type": "bot"}
                )
            },
        ) as session:
            search = partial(
                self.tns_post, session, semaphore, "get/search", self.token
            )
            get = partial(self.tns_post, session, semaphore, "get/object", self.token)
            response = await search(params)
            names = [item["objname"] for item in response["data"]["reply"]]
            tasks = [
                get({"objname": name})
                for name in names
                if not (exclude and name in exclude)
            ]
            self.logger.info(
                "TNS search", extra={"params": params, "results": len(names)}
            )
            for fut in asyncio.as_completed(tasks):
                item = await fut
                yield item["data"]["reply"]

    async def get(self, session, semaphore, objname):
        return await self.tns_post(
            session, semaphore, "get/object", self.token, {"objname": objname}
        )

    async def sendReport(self, report):
        semaphore = asyncio.Semaphore(self.maxParallelRequests)
        async with ClientSession(
            headers={
                "User-Agent": "tns_marker"
                + json.dumps(
                    {"tns_id": self.token.id, "name": self.token.name, "type": "bot"}
                )
            },
        ) as session:
            postreport = partial(
                self.tns_post, session, semaphore, "bulk-report", self.token
            )
            response = await postreport(report)
            print(response)
            if response["id_code"] == 200:
                return response["data"]["report_id"]
            return False

    async def reportReply(self, report_id):
        semaphore = asyncio.Semaphore(self.maxParallelRequests)
        async with ClientSession(
            headers={
                "User-Agent": "tns_marker"
                + json.dumps(
                    {"tns_id": self.token.id, "name": self.token.name, "type": "bot"}
                )
            },
        ) as session:
            reply_data = {"api_key": self.token.api_key, "report_id": report_id}
            print("replydata", reply_data)
            postreport = partial(
                self.tns_post, session, semaphore, "bulk-report-reply", self.token
            )
            response = await postreport(reply_data)
            return response["data"]["feedback"]
