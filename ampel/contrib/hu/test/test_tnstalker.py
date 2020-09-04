import json
from os.path import dirname, join

import pytest

import ampel.contrib.hu.t3.TNSTalker as TNSTalker
from ampel.log.AmpelLogger import AmpelLogger
from ampel.util.legacy.json_v06 import object_hook


# Load mixed set of transients for check
@pytest.fixture
def t3_TNStvs_mixed():
    with open(join(dirname(__file__), "t3TransientTalker_testTV.json")) as f:
        return [json.loads(line, object_hook=object_hook) for line in f]


@pytest.fixture
def t3_TNStvs_good():
    with open(join(dirname(__file__), "t3TransientTalker_testsubmitTV.json")) as f:
        return [json.loads(line, object_hook=object_hook) for line in f]


@pytest.fixture
def testrunconfig():
    return {
        "logger": AmpelLogger.get_logger(),
        "submit_tns": False,
        "sandbox": False,
        "needed_catalogs": [],
    }


def test_run_t3_TnsTalker_selection(t3_TNStvs_mixed, testrunconfig):
    """
    Check whether correct TVs are accepted based on transient information
    """
    unit_test = TNSTalker.TNSTalker(**testrunconfig)
    out = [unit_test.accept_tview(tv) for tv in t3_TNStvs_mixed]
    assert len(out) == 30
    assert sum(out) == 5


@pytest.mark.xfail(reason="tries to connect to TNS")
def test_run_t3_TnsTalker_tnsnamefound(t3_TNStvs_good, testrunconfig, mocker):
    """
    Check whether it find all of these to be submitted. Note - can connect to remove TNS DB!
    """
    mock = mocker.patch(
        "ampel.contrib.hu.t3.TNSTalker.get_tnsname",
        side_effect=AssertionError("get_tnsname should not be called"),
    )
    unit_test = TNSTalker.TNSTalker(**testrunconfig)
    for tv in t3_TNStvs_good:
        tns_name, tns_internals, jup = unit_test.find_tns_name(tv)
        assert not mock.called
        assert tns_name is not None
