import json
from functools import partial
from pathlib import Path

import fastavro
import numpy as np
import pytest

from ampel.alert.load.TarAlertLoader import TarAlertLoader
from ampel.contrib.hu.augment.ZiArchiveAugmenter import ZiArchiveAugmenter
from ampel.log.AmpelLogger import DEBUG, AmpelLogger
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.ingest.ZiDataPointShaper import ZiDataPointShaperBase

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
        muxer = "ZiMongoMuxer"
        sig1 = 0.1
        sig2 = 0.2
        nu1 = 1e5
        min_p = 0.9
        augmenter = ZiArchiveAugmenter(
            context=mock_context,
            logger=logger,
            tag="ZTF",
            augmenting_shaper=shaper,
            archive_token=token,
            mux=muxer,
            min_posterior=min_p,
            nu1=nu1,  # sources per sqdeg
            sigma1=sig1,  # resolution in arcsec of primary survey
            sigma2=sig2,  # resolution in arcsec of secondary survey (both ZTF in this case)
        )

        # double check math
        rho1 = 4 * np.pi * nu1 / (np.pi / 180) ** 2
        sig1sig2 = (sig1**2 + sig2**2) * (np.pi / 180 / 3600) ** 2
        psi_rad = np.sqrt(np.log((1 / min_p - 1) / rho1 * 2 / sig1sig2) * 2 * sig1sig2)
        psi_arcsec = psi_rad * 180 * 3600 / np.pi
        rel_diff = (psi_arcsec - augmenter._radius_arcsec) / psi_arcsec
        assert abs(rel_diff) < 1e-15

        # patch archive call
        augmenter.session.get = get_archive_result

        c, i = augmenter.process(dps, stock_id=alert.stock)

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
