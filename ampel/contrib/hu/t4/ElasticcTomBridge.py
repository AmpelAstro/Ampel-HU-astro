#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t4/ElasticcTomBridge.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                12.02.2025
# Last Modified Date:  12.02.2025
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

from collections.abc import Generator
from math import isfinite

import fastavro
from pydantic import TypeAdapter

from ampel.abstract.AbsT4Unit import AbsT4Unit
from ampel.base.AuxUnitRegister import AuxUnitRegister
from ampel.lsst.alert.load.KafkaAlertLoader import KafkaAlertLoader
from ampel.model.UnitModel import UnitModel
from ampel.util.collections import get_chunks

from ..t3.ElasticcTomClient import ElasticcClassification, ElasticcTomClient


class ElasticcTomBridge(AbsT4Unit):
    """
    Bridge between Elasticc classifications and the DESC TOM system.

    This T4 unit reads Elasticc classifications from a Kafka topic and
    posts them to the DESC TOM system.

    """

    desc_username: str
    desc_password: str

    tom_url: str = "https://desc-tom.lbl.gov"
    endpoint: str = "/elasticc/brokermessage/"

    timestamp_field: str = "brokerIngestTimestamp"
    chunk_size: int = 1000

    dry_run: bool = False

    loader: UnitModel

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tom_client = ElasticcTomClient(
            tom_url=self.tom_url,
            desc_username=self.desc_username,
            desc_password=self.desc_password,
            endpoint=self.endpoint,
            logger=self.logger,
        )
        self.consumer = AuxUnitRegister.new_unit(self.loader, sub_type=KafkaAlertLoader)

    def chunks(self) -> Generator[list[ElasticcClassification], None, None]:
        """
        Yield chunks of messages, forever
        """
        validator = TypeAdapter(ElasticcClassification).validate_python
        # disable logical type conversion, in particular int -> datetime for timestamp-millis
        fastavro.read.LOGICAL_READERS.clear()
        for chunk in get_chunks(self.consumer, self.chunk_size):
            meta_records = [message.pop("__kafka") for message in chunk]
            yield [
                validator(
                    {self.timestamp_field: meta["timestamp"]["created"], **message}
                )
                for message, meta in zip(chunk, meta_records, strict=True)
            ]
            self.consumer.acknowledge(meta_records)
        self.logger.debug("no more chunks")

    def do(self) -> None:
        try:
            for chunk in self.chunks():
                for report in chunk:
                    for classification in report["classifications"]:
                        if not isfinite(classification["probability"]):
                            classification["probability"] = 0
                if self.dry_run:
                    self.logger.info(f"would post {len(chunk)} classifications")
                    continue
                self.logger.info(f"posting {len(chunk)} classifications")
                response = self.tom_client.tom_post(chunk)
                if not response["success"]:
                    self.logger.error(response["response_body"])
                    raise RuntimeError(
                        f"POST failed with status {response['response']}"
                    )
        except KeyboardInterrupt:
            self.logger.info("exiting")
