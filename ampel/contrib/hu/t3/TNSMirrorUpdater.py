#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/TNSMirrorUpdater.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                13.12.2018
# Last Modified Date:  29.12.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import datetime
from typing import Any

from pymongo import MongoClient

from ampel.abstract.AbsOpsUnit import AbsOpsUnit
from ampel.contrib.hu.t3.tns.TNSClient import TNSClient
from ampel.contrib.hu.t3.tns.TNSMirrorDB import TNSMirrorDB
from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.secret.NamedSecret import NamedSecret


class TNSMirrorUpdater(AbsOpsUnit):
    """
    Sync a local mirror of the TNS database
    """

    extcats_auth: NamedSecret[dict] = NamedSecret[dict](label="extcats/writer")
    api_key: NamedSecret[dict]
    timeout: float = 60.0
    max_parallel_requests: int = 8
    dry_run: bool = False

    def run(self, beacon: None | dict[str, Any] = None) -> None | dict[str, Any]:
        now = datetime.datetime.now(tz=datetime.timezone.utc)
        last_run = beacon["updated"] if beacon else now - datetime.timedelta(days=7)

        def fetch():
            tns = TNSClient(
                TNSToken(**self.api_key.get()),
                self.timeout,
                self.max_parallel_requests,
                self.logger,
            )
            return tns.search(public_timestamp=str(last_run))

        if not self.dry_run:
            new_reports = fetch()
            TNSMirrorDB(
                MongoClient(
                    self.context.config.get("resource.extcats", str, raise_exc=True),
                    **self.extcats_auth.get(),
                ),
                logger=self.logger,
            ).add_sources(new_reports)

        return {"updated": now}
