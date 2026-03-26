import json
from functools import partial
from pathlib import Path

import fastavro
import matplotlib.pyplot as plt
import mongomock
import numpy as np
import pymongo
import pytest

from ampel.alert.load.TarAlertLoader import TarAlertLoader
from ampel.contrib.hu.t0.ZiArchiveAdder import ZiArchiveAdder
from ampel.contrib.hu.t1.T1PositionalStreamCombine import T1PositionalStreamCombine
from ampel.log.AmpelLogger import DEBUG, AmpelLogger
from ampel.model.operator.AnyOf import AnyOf
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.ingest.ZiDataPointShaper import ZiDataPointShaperBase
from ampel.ztf.util.ZTFIdMapper import to_ztf_id

DATA_DIR = Path(__file__).parent / "data"
ZTF_TEST_NAME = "ZTF18abxhyqv"


# --------------------- copied from Ampel-ZTF --------------------- #


@pytest.fixture
def avro_packets():
    """
    4 alerts for a random AGN, widely spaced:

    ------------------ -------------------------- ------------------------
    candid             detection                  history
    ------------------ -------------------------- ------------------------
    673285273115015035 2018-11-05 06:50:48.001935 29 days, 22:11:31.004165
    879461413115015009 2019-05-30 11:04:25.996800 0:00:00
    882463993115015007 2019-06-02 11:08:09.003839 3 days, 0:03:43.007039
    885458643115015010 2019-06-05 11:00:26.997131 5 days, 23:56:01.000331
    ------------------ -------------------------- ------------------------
    """
    return partial(
        TarAlertLoader,
        file_path=str(Path(__file__).parent / "data" / f"{ZTF_TEST_NAME}.tar.gz"),
    )


@pytest.fixture
def raw_alert_dicts(avro_packets):
    def gen():
        for f in avro_packets():
            yield next(fastavro.reader(f))

    return gen


@pytest.fixture
def alerts(raw_alert_dicts):
    def gen():
        for d in raw_alert_dicts():
            yield ZiAlertSupplier.shape_alert_dict(d)

    return gen


# ----------------------------------------------------------------- #


def get_archive_result(*args, **kwargs):
    archive_result_file = (
        Path(__file__).parent / "data" / f"{ZTF_TEST_NAME}_archive_1arcmin.json"
    )
    with archive_result_file.open() as f:
        return json.load(f)


@pytest.fixture
def ztf_archive_adder(mock_context):
    logger = AmpelLogger.get_logger(console=dict(level=DEBUG))
    shaper = {"unit": "ZiDataPointShaper"}
    token = {"label": "ztf/archive/token", "value": ""}
    muxer = "ZiMongoMuxer"
    adder = ZiArchiveAdder(
        context=mock_context,
        logger=logger,
        augmenting_shaper=shaper,
        archive_token=token,
        mux=muxer,
        radius_arcsec=1,
    )
    adder.session.get = get_archive_result
    return adder


def test_ztf_archive_adder(alerts, ztf_archive_adder):
    """Just test functionality of the ZiArchiveAdder"""
    adder = ztf_archive_adder
    for alert in alerts():
        # convert to DataPoint
        dps = ZiDataPointShaperBase().process(alert.datapoints, alert.stock)

        i, c = adder.process(dps, stock_id=alert.stock)

        # check that all and only the correct archive points were added
        archive_result = get_archive_result()
        correct_archive_candids = [a["candid"] for a in archive_result]
        selected_candids = [ic["id"] for ic in c]
        assert all(
            np.isin(selected_candids, correct_archive_candids)
            | np.isin(selected_candids, [dp["id"] for dp in dps])
        )
        assert all(np.isin(correct_archive_candids, selected_candids))


def test_positional_stream_combine(alerts, ztf_archive_adder, monkeypatch, tmp_path):
    logger = AmpelLogger.get_logger(console=dict(level=DEBUG))
    sig1 = 0.1
    sig2 = 0.2
    nu1 = 1e5
    min_p = 0.9
    primary_tag = "ZTF1"
    secondary_tag = "ZTF2"

    monkeypatch.setattr(pymongo, "MongoClient", mongomock.MongoClient)

    t1 = T1PositionalStreamCombine(
        logger=logger,
        sigma1=sig1,
        sigma2=sig2,
        nu1=nu1,
        min_posterior=min_p,
        primary_tag=AnyOf(any_of=[primary_tag]),
        secondary_tag=AnyOf(any_of=[secondary_tag]),
        mongo_uri="mongodb://localhost:27017",
        database_name="test_db",
    )

    # patch tag
    ztf_archive_adder.tag = secondary_tag

    for i, alert in enumerate(alerts()):
        dps = ZiDataPointShaperBase().process(alert.datapoints, alert.stock)
        for dp in dps:
            dp["tag"] = [primary_tag, *list(dp["tag"])]
        _, archive_dps = ztf_archive_adder.process(dps, stock_id=alert.stock)
        for dp in archive_dps:
            dp["tag"] = [secondary_tag, *list(dp["tag"])]
        t1res = t1.combine(dps + archive_dps)

        fig, ax = plt.subplots()
        ax.scatter(
            [dp["body"]["ra"] for dp in dps if "ra" in dp["body"]],
            [dp["body"]["dec"] for dp in dps if "dec" in dp["body"]],
            label="alert",
            marker="x",
            zorder=10,
        )
        ax.scatter(
            [
                dp["body"]["ra"]
                for dp in archive_dps
                if "ra" in dp["body"] and dp["id"] in t1res.dps
            ],
            [
                dp["body"]["dec"]
                for dp in archive_dps
                if "dec" in dp["body"] and dp["id"] in t1res.dps
            ],
            label="archive",
            marker="s",
        )
        ax.set_aspect("equal")
        ax.legend()
        fig.savefig(tmp_path / f"{i}.png")
        plt.close()

        # check that the right source was selected
        assert to_ztf_id(t1res.meta["stock"]) == ZTF_TEST_NAME
        assert t1res.meta["p_association"] > min_p

        # check that all and only the correct archive points were added
        archive_result = get_archive_result()
        de_muxed_dp_ids = [
            1388394243115015014,
            1388290453115015012,
            1324402953115015007,
        ]
        correct_archive_candids = [
            a["candid"]
            for a in archive_result
            if (a["objectId"] == ZTF_TEST_NAME) and (a["candid"] not in de_muxed_dp_ids)
        ]
        selected_candids = [ic["id"] for ic in archive_dps if ic["id"] in t1res.dps]
        assert all(
            np.isin(selected_candids, correct_archive_candids)
            | np.isin(selected_candids, [dp["id"] for dp in dps])
        )
        assert all(np.isin(correct_archive_candids, selected_candids))
