from confluent_kafka import Producer
from fastavro import schemaless_writer
from typing import Any
from io import BytesIO

from ampel.abstract.AbsUnitResultAdapter import AbsUnitResultAdapter
from ampel.base.AmpelUnit import AmpelUnit
from ampel.struct.UnitResult import UnitResult
from ampel.util.mappings import get_by_path

from ampel.lsst.alert.load.HttpSchemaRepository import parse_schema


class KafkaReporter(AmpelUnit):

    broker: str
    topic: str
    avro_schema: dict | str

    producer_config: dict[str, Any] = {}

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._schema = parse_schema(self.avro_schema)
        self._producer = Producer(
            **{
                "bootstrap.servers": self.broker,
            }
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
        self._producer.flush()
        self._producer.poll(0)


class KafkaAdapter(AbsUnitResultAdapter, KafkaReporter):

    #: Where to find message in UnitResult.body
    # NB: use list[int | str] to prevent coercion to str
    message_path: int | str | list[int | str]
    raise_exc: bool = False

    def handle(self, ur: UnitResult) -> UnitResult:
        if isinstance(message := get_by_path(ur.body, self.message_path), dict):
            self.send(message)
            self.flush()
        elif self.raise_exc:
            raise KeyError(f"{self.message_path} not found in {ur.body!r}")

        return ur
