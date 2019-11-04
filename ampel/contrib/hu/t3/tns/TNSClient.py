#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/aiotns
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 13.12.2018
# Last Modified Date: 04.11.2019
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

import aiohttp
import asyncio
import backoff

from functools import partial
from aiohttp import ClientSession, ClientTimeout
from aiohttp.client_exceptions import ServerDisconnectedError, ClientConnectorError, ClientConnectionError, ClientResponseError

async def tns_post(session, semaphore, method, api_key, data):
    """
    post to TNS, asynchronously
    """
    async with semaphore:
        with aiohttp.MultipartWriter('form-data') as mpwriter:
            p = aiohttp.StringPayload(api_key)
            p.set_content_disposition('form-data', name='api_key')
            mpwriter.append(p)
            p = aiohttp.JsonPayload(data)
            p.set_content_disposition('form-data', name='data')
            mpwriter.append(p)
            resp = await session.post('https://wis-tns.weizmann.ac.il/api/'+method, data=mpwriter)
        return await resp.json()

class TNSClient:
    def __init__(self, apiKey, timeout, maxParallelRequests, logger):
        self.apiKey = apiKey
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
            max_time=timeout
        )(tns_post)

    def _on_backoff(self, details):
        exc_typ, exc, _ = sys.exc_info()
        if exc is not None:
            exc = traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
        self.logger.warn('backoff', extra={'exc': exc, 'wait': details['wait'], 'tries': details['tries']})

    def _on_giveup(self, details):
        exc_typ, exc, _ = sys.exc_info()
        if exc is not None:
            exc = traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
        self.logger.warn('gave up', extra={'exc': exc, 'tries': details['tries']})

    @staticmethod
    def is_permanent_error(exc):
        if isinstance(exc, ClientResponseError):
            return exc.code not in {500}
        else:
            return False

    async def search(self, exclude=set(), **params):
        semaphore = asyncio.Semaphore(self.maxParallelRequests)
        async with ClientSession(raise_for_status=True) as session:
            search = partial(self.tns_post, session, semaphore, 'get/search', self.apiKey)
            get = partial(self.tns_post, session, semaphore, 'get/object', self.apiKey)
            response = await search(params)
            names = [item['objname'] for item in response['data']['reply']]
            tasks = [get({'objname': name}) for name in names if not (exclude and name in exclude)]
            self.logger.info("TNS search", extra={'params': params, 'results': len(names)})
            for fut in asyncio.as_completed(tasks):
                item = await fut
                yield item['data']['reply']

    async def get(self, session, semaphore, objname):
        return await self.tns_post(session, semaphore, 'get/object', self.apiKey, {'objname': objname})
