import json
from pathlib import Path
from functools import partial

import numpy as np
import pytest
import fastavro

from ampel.contrib.hu.ingest.ZiArchiveAugmenter import ZiArchiveAugmenter
from ampel.log.AmpelLogger import DEBUG, AmpelLogger
from ampel.test.conftest import mock_context, _patch_mongo
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.dev.ZTFAlert import ZTFAlert
from ampel.ztf.ingest.ZiDataPointShaper import ZiDataPointShaperBase
from ampel.alert.load.TarAlertLoader import TarAlertLoader


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


def test_ztf_archive_augmenter(alerts, mock_context):
    for alert in alerts():
        # convert to DataPoint
        dps = ZiDataPointShaperBase().process(alert.datapoints, alert.stock)
        logger = AmpelLogger.get_logger(console=dict(level=DEBUG))
        shaper = {"unit": "ZiDataPointShaper"}
        token = {"label": "ztf/archive/token", "value": ""}
        muxer = ZiArchiveAugmenter(
            context=mock_context,
            logger=logger,
            tag="ZTF",
            augmenting_shaper=shaper,
            archive_token=token,
        )

        # patch archive call
        muxer.session.get = get_archive_result

        c, i = muxer.process(dps, stock_id=alert.stock)

        # check that all and only the correct archive points were added
        archive_result = get_archive_result()
        correct_archive_candids = [
            a["candid"] for a in archive_result if a["objectId"] == ZTF_TEST_NAME
        ]
        selected_candids = [ic["id"] for ic in c]
        assert all(
            np.isin(selected_candids, correct_archive_candids)
            | np.isin(selected_candids, [dp["id"] for dp in dps])
        )
        assert all(np.isin(correct_archive_candids, selected_candids))
