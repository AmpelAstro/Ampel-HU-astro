import importlib
import json
from functools import cached_property
from typing import Annotated, Any

import py_avro_schema
from fastavro.schema import parse_schema
from hop.auth import Auth
from hop.io import Producer
from hop.models import AvroBlob
from pydantic import BaseModel, BeforeValidator

from ampel.abstract.AbsUnitResultAdapter import AbsUnitResultAdapter
from ampel.base.AmpelUnit import AmpelUnit
from ampel.lsst.kafka.KafkaAuthentication import SASLAuthentication
from ampel.struct.UnitResult import UnitResult
from ampel.util.mappings import get_by_path


def _get_model(value: Any) -> type[BaseModel]:
    if isinstance(value, type) and issubclass(value, BaseModel):
        return value
    if isinstance(value, str):
        mod, attr = value.rsplit(".", 1)
        if not mod.startswith("ampel."):
            raise TypeError("Only models from the ampel package are supported")
        model = getattr(importlib.import_module(mod), attr)
        if isinstance(model, type) and issubclass(model, BaseModel):
            return model
        raise TypeError(f"{value} is not a BaseModel subclass")
    raise TypeError("model must be a BaseModel subclass a fully-qualified name")


class HopskotchProducer(AmpelUnit):
    broker: str = "kafka.scimma.org"
    auth: SASLAuthentication
    topic: str

    model: Annotated[type[BaseModel], BeforeValidator(_get_model)]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._producer = Producer(
            broker_addresses=[self.broker],
            topics=[self.topic],
            auth=Auth(self.auth.username.get(), self.auth.password.get()),
            automatic_offload=False,
        )

    @cached_property
    def schema(self) -> dict:
        return parse_schema(  # type: ignore[return-value]
            json.loads(
                py_avro_schema.generate(
                    self.model, options=py_avro_schema.Option.AUTO_NAMESPACE_MODULE
                ).decode()
            )
        )

    def send(self, record: dict) -> None:
        self._producer.write(
            message=AvroBlob(content=record, schema=self.schema, single_record=True),
            topic=self.topic,
        )

    def flush(self):
        if (in_queue := self._producer.flush()) > 0:
            raise TimeoutError(f"{in_queue} messages still in queue")


class HopskotchAdapter(AbsUnitResultAdapter, HopskotchProducer):
    #: Where to find message in UnitResult.body
    # NB: use list[int | str] to prevent coercion to str
    message_path: None | int | str | list[int | str] = None
    raise_exc: bool = False
    #: Clear the body before returning
    drop_body: bool = False

    def _get_by_path(self, body: dict) -> Any:
        if self.message_path is None:
            return body
        return get_by_path(body, self.message_path)

    def handle(self, ur: UnitResult) -> UnitResult:
        assert isinstance(ur.body, dict)
        if isinstance(message := self._get_by_path(ur.body), dict):
            self.send(message)
            self.flush()
        elif isinstance(message, list):
            for part in message:
                self.send(part)
            self.flush()
        elif self.raise_exc:
            raise KeyError(f"{self.message_path} not found in {ur.body!r}")
        if self.drop_body:
            ur.body = None

        return ur
