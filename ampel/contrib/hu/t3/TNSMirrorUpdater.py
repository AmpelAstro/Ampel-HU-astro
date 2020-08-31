#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TNSMirrorUpdater.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 13.12.2018
# Last Modified Date: 13.08.2020
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

import asyncio
import datetime

from ampel.abstract.AbsT3Unit import AbsT3Unit

from ampel.contrib.hu.t3.tns.TNSClient import TNSClient
from ampel.contrib.hu.t3.tns.TNSMirrorDB import TNSMirrorDB


class TNSMirrorUpdater(AbsT3Unit):
    """
    Sync a local mirror of the TNS database
    """

    require = ("extcats.writer",)

    api_key: str
    timeout: float = 60.0
    max_parallel_requests: int = 8
    dry_run: bool = False

    def add(self, transients):
        ...

    def done(self) -> None:
        if self.context and "last_run" in self.context:
            last_run = datetime.datetime.fromtimestamp(self.context["last_run"])
        else:
            last_run = datetime.datetime(2020,9,1)

        async def fetch():
            tns = TNSClient(
                self.api_key, self.timeout, self.max_parallel_requests, self.logger
            )
            return [item async for item in tns.search(public_timestamp=str(last_run))]

        new_reports = asyncio.get_event_loop().run_until_complete(fetch())

        assert self.resource
        db = TNSMirrorDB(self.resource["extcats.writer"], logger=self.logger)
        db.add_sources(new_reports)
