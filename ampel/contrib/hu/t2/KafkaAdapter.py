from io import BytesIO
from typing import Any

from confluent_kafka import Producer
from fastavro import schemaless_writer

from ampel.abstract.AbsUnitResultAdapter import AbsUnitResultAdapter
from ampel.base.AmpelUnit import AmpelUnit
from ampel.lsst.alert.load.HttpSchemaRepository import parse_schema
from ampel.lsst.alert.load.KafkaAlertLoader import SASLAuthentication
from ampel.struct.UnitResult import UnitResult
from ampel.util.mappings import get_by_path


class KafkaReporter(AmpelUnit):
    broker: str
    topic: str
    avro_schema: dict | str

    auth: None | SASLAuthentication = None

    producer_config: dict[str, Any] = {}
    delivery_timeout: float = 10.0

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._schema = parse_schema(self.avro_schema)
        self._producer = Producer(
            **{
                "bootstrap.servers": self.broker,
            }
            | (self.auth.librdkafka_config() if self.auth else {})
            | self.producer_config
        )

    def serialize(self, record: dict) -> bytes:
        buf = BytesIO()
        schemaless_writer(buf, self._schema, record)
        return buf.getvalue()

    def send(self, record: dict) -> None:
        self._producer.poll(0)
        self._producer.produce(self.topic, self.serialize(record))

    def flush(self):
        if (in_queue := self._producer.flush(self.delivery_timeout)) > 0:
            raise TimeoutError(
                f"{in_queue} messages still in queue after {self.delivery_timeout}s"
            )
        self._producer.poll(0)


class KafkaAdapter(AbsUnitResultAdapter, KafkaReporter):
    #: Where to find message in UnitResult.body
    # NB: use list[int | str] to prevent coercion to str
    message_path: int | str | list[int | str]
    raise_exc: bool = False
    #: Clear the body before returning
    drop_body: bool = False

    def handle(self, ur: UnitResult) -> UnitResult:
        assert isinstance(ur.body, dict)
        if isinstance(message := get_by_path(ur.body, self.message_path), dict):
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
