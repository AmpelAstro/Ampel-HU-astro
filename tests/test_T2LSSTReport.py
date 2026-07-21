from collections.abc import Mapping
from pathlib import Path

import bson
import numpy as np
import pytest

from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2LSSTReport import Classification, T2LSSTReport
from ampel.kafka.HopskotchAdapter import MessageSerializer


@pytest.fixture
def bundle():
    return bson.decode(
        (Path(__file__).parent / "data" / "25409136044802058.bson").read_bytes()
    )


@pytest.fixture
def datapoints(bundle):
    dps = bundle["datapoints"]
    for dp in reversed(dps):
        if "LSST_OBJ" in dp["tag"]:
            # ensure object record fields are not NaN
            dp["body"] = {
                k: v if not isinstance(v, float) or not np.isnan(v) else -1
                for k, v in dp["body"].items()
            }
            break
    return dps


@pytest.fixture
def compound(bundle):
    return T1Document(
        dps=[dp["id"] for dp in bundle["datapoints"]],
        stock=bundle["datapoints"][0]["stock"][0],
        link=1,
    )


class DummyT2LSSTReport(T2LSSTReport):
    def process_t2s(self, payloads: dict[str, Mapping[str, any]]):
        return [
            Classification(
                name="test_classifier",
                version="1.0",
                info="Test classification",
                models=[
                    {
                        "model": "test_model",
                        "probabilities": {"class_A": 0.8, "class_B": 0.2},
                    }
                ],
            )
        ]


def test_process(compound, datapoints, ampel_logger, mock_context):

    t2 = DummyT2LSSTReport(logger=ampel_logger, t2_dependency=[])
    result = t2.process(compound, datapoints, [])
    assert result.code is None
    assert result.journal is None
    assert isinstance(result.body, dict)

    # message serialization test
    serializer = MessageSerializer(model="ampel.contrib.hu.model.LSSTReport.LSSTReport")
    message = serializer.serialize(result.body)
    blob = message.serialize()["content"]
    assert isinstance(blob, bytes)
    assert message.deserialize(blob).content == message.content

    # ensure that ids were not corrupted by secret rounding
    dp = {
        doc["body"]["diaSourceId"]: doc for doc in datapoints if "LSST_DP" in doc["tag"]
    }
    fp = {
        doc["body"]["diaForcedSourceId"]: doc
        for doc in datapoints
        if "LSST_FP" in doc["tag"]
    }
    dps = {"LSST_DP": dp, "LSST_FP": fp}
    for dp in message.content["photometry"]:
        assert dp["id"] in dps[dp["source"]]
