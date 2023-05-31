#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/ElasticcTomClient.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                11.04.2022
# Last Modified Date:  11.04.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

from typing import Generator, Sequence, Dict, Any, Union

import logging
import requests
import json
import backoff
from requests import HTTPError
from requests.models import InvalidJSONError
from io import BytesIO
from math import isfinite

from ampel.lsst.alert.load.HttpSchemaRepository import parse_schema
from ampel.ztf.t0.load.AllConsumingConsumer import AllConsumingConsumer
import fastavro
from ampel.util.collections import get_chunks

class ClassificationDict():
    classifierName: str
    classifierParams: str
    classId: int
    probability: float

class ElasticcClassification():
    alertId: int
    diaSourceId: int
    elasticcPublishTimestamp: int
    brokerIngestTimestamp: int
    brokerName: str
    brokerVersion: str
    classifications: Sequence[ClassificationDict]


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
        endpoint: str = "/elasticc/brokermessage/"
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
        self.session.get( f'{self.tom_url}/accounts/login/' )
        self.session.post( f'{self.tom_url}/accounts/login/',
              data={ "username": desc_username,
                     "password": desc_password,
                     "csrfmiddlewaretoken": self.session.cookies['csrftoken'] } )
        self.csrfheader = { 'X-CSRFToken': self.session.cookies['csrftoken'] }



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
        giveup=lambda e: not isinstance(e, HTTPError) or e.response.status_code not in {503, 504, 429, 408},
        max_time=60,
        )
    def tom_post(self, classification: Union[ElasticcClassification,list[ElasticcClassification]])->Dict[Any,Any]:
        try:
            response = self.session.put(f'{self.tom_url}{self.endpoint}',
                                    json=classification, headers=self.csrfheader)
        except InvalidJSONError:
            for report in (classification if isinstance(classification, list) else [classification]):
                try:
                    json.dumps(report, allow_nan=False)
                except ValueError:
                    self.logger.error(f"invalid report: {report}")
            raise

        if response.ok:
            self.logger.debug('ElasticcTomClient submit done', extra={"payload": classification})
            return {'success':True, **response.json()}

        self.logger.info('ElasticcTomClient submit fail', extra={"payload": classification})
        return {'success':False, 'response':response.status_code, 'response_body': response.text}

logger = logging.getLogger()

class AvroDeserializer:
    def __init__(self, schema: str | dict):
        self._schema = parse_schema(schema)
    def __call__(self, message) -> dict:
        return fastavro.schemaless_reader(
            BytesIO(message.value()), self._schema
        )

def chunks_from_kafka(broker: str, topic: str, group_id: str, avro_schema: str | dict, timeout: float=30, chunk_size: int=1000, **consumer_config) -> Generator[list[dict], None, None]:
    """
    Yield chunks of messages a topic, forever
    """
    consumer = AllConsumingConsumer(broker, timeout=timeout, topics=[topic], auto_commit=False, logger=logger, **{"group.id": group_id}, )
    deserializer = AvroDeserializer(avro_schema)
    while True:
        for chunk in get_chunks(map(deserializer, consumer), chunk_size):
            yield chunk
            consumer.commit()
        logger.debug("no more chunks")

if __name__ == "__main__":
    from argparse import ArgumentParser

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(name)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    parser = ArgumentParser()
    parser.add_argument("--bootstrap")
    parser.add_argument("--topic")
    parser.add_argument("--group")
    parser.add_argument("--schema")
    parser.add_argument("--tom")
    parser.add_argument("-u", "--user")
    parser.add_argument("-p", "--password")
    parser.add_argument("--endpoint", default="/elasticc2/brokermessage/")
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--chunk-size", type=int, default=1000)

    args = parser.parse_args()

    tom_client = ElasticcTomClient(
        tom_url=args.tom,
        desc_username=args.user,
        desc_password=args.password,
        endpoint=args.endpoint,
        logger=logger
    )

    # disable logical type conversion, in particular int -> datetime for timestamp-millis
    fastavro.read.LOGICAL_READERS.clear()
    try:
        for chunk in chunks_from_kafka(broker=args.bootstrap, topic=args.topic, group_id=args.group, avro_schema=args.schema, timeout=args.timeout, chunk_size=args.chunk_size):
            for report in chunk:
                for classification in report["classifications"]:
                    if not isfinite(classification["probability"]):
                        classification["probability"] = 0
            response = tom_client.tom_post(chunk)
            if not response["success"]:
                logger.error(response["response_body"])
                raise RuntimeError(f"POST failed with status {response['response']}")
            logger.info(f"posted {len(chunk)} classifications")
    except KeyboardInterrupt:
        logger.info("exiting")

    