from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
import json
from os import environ
from os.path import dirname, join

import pytest

from ampel.contrib.hu.t3.TNSTalker import TNSTalker, TNSClient, TNS_BASE_URL_SANDBOX
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
    unit_test = TNSTalker(**testrunconfig)
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
    unit_test = TNSTalker(**testrunconfig)
    for tv in t3_TNStvs_good:
        tns_name, tns_internals, jup = unit_test.find_tns_name(tv)
        assert not mock.called
        assert tns_name is not None


def test_tnsclient():
    if not (api_key := environ.get("TNS_API_KEY")):
        raise pytest.skip("Test requires env var TNS_API_KEY")
    client = TNSClient(
        TNS_BASE_URL_SANDBOX,
        AmpelLogger.get_logger(),
        TNSToken(
            **{
                "id": 59228,
                "name": "ZTF_AMPEL_COMPLETE",
                "api_key": api_key,
            }
        ),
    )
    assert client.getInternalName("2018cow") == ('ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa', 'Got internal name response')
    assert client.search(244.000917, 22.268031) == (['SN2018cow'], 'Found TNS name(s)')
    assert client.getNames(244.000917, 22.268031) == ('SN2018cow', 'ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa')
