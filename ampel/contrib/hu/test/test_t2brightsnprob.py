import json
from os.path import dirname, join

import pytest

from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t2.T2BrightSNProb import T2BrightSNProb
from ampel.log.AmpelLogger import AmpelLogger
from ampel.view.LightCurve import LightCurve
from ampel.ztf.util.ZTFIdMapper import to_ampel_id

_ampel_flags = {
    0: "ZTF",
}

_ampel_photo_flags = _ampel_flags | {
    20: "PHOTOPOINT",
    21: "UPPERLIMIT",
    22: "ZTF_G",
    23: "ZTF_R",
    24: "ZTF_I",
    25: "ZTF_PRIV",
    26: "ZTF_PUB",
}


def _to_datapoint(d):
    return DataPoint(
        id=d.pop("_id"),
        body=d,
        tag=[_ampel_photo_flags[f] for f in d.pop("alTags") if f in _ampel_photo_flags],
        channel=[],
        meta=[],
    )


@pytest.fixture()
def lightcurve() -> LightCurve:
    with open(join(dirname(__file__), "lightcurve.ZTF18abmjvpb.json")) as f:
        blob = json.load(f)
        return LightCurve(
            compound_id=0,
            stock_id=to_ampel_id("ZTF18abmjvpb"),
            photopoints=[
                _to_datapoint(item["content"]["__jsonclass__"][1])
                for item in blob["ppo_list"]["__jsonclass__"][1]
            ],
            upperlimits=[
                _to_datapoint(item["content"]["__jsonclass__"][1])
                for item in blob["ulo_list"]["__jsonclass__"][1]
            ],
        )


def assert_equivalent(left, right):
    # split into optional and valued:  https://github.com/pytest-dev/pytest/pull/7710
    assert {k: v for k, v in left.items() if v is None} == {
        k: v for k, v in right.items() if v is None
    }
    assert {k: v for k, v in left.items() if v is not None} == pytest.approx(
        {k: v for k, v in right.items() if v is not None}
    )


def test_t2brightsnprob(lightcurve):
    monitor = T2BrightSNProb(logger=AmpelLogger.get_logger())  # type: ignore[call-arg]
    monitor.post_init()
    assert_equivalent(
        monitor.process(lightcurve),
        {
            "cut_pp": 0,
            "jd_det": 2458343.6521875,
            "jd_last": 2458373.62625,
            "ndet": 10,
            "mag_det": 19.96540069580078,
            "mag_last": 18.89119529724121,
            "t_lc": 29.97406250005588,
            "rb_med": 0.6449998319149017,
            "drb_med": None,
            "distnr_med": 2.2441248893737793,
            "magnr_med": 22.28499984741211,
            "classtar_med": 0.9884999990463257,
            "sgscore1_med": 0.5,
            "distpsnr1_med": 15.384442329406738,
            "sgscore2_med": 0.5,
            "distpsnr2_med": 17.98440933227539,
            "neargaia_med": 60.776649475097656,
            "maggaia_med": 17.757875442504883,
            "bool_pure": True,
            "t_predetect": 3.998587999958545,
            "bool_peaked": True,
            "jd_max": 2458360.7482292,
            "mag_peak": 18.341100692749023,
            "bool_rising": False,
            "bool_norise": False,
            "bool_hasgaps": False,
            "slope_rise_g": None,
            "slope_fall_g": None,
            "slope_rise_r": -0.06971756172380136,
            "slope_fall_r": 0.04324959751054604,
            "col_det": None,
            "col_last": None,
            "col_peak": None,
            "SNGuess": 6.749870969880001,
            "SNGuessBool": 1,
            "success": True,
        },
    )
