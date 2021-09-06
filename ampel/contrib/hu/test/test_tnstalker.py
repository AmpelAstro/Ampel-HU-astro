from os import environ

import pytest

from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.contrib.hu.t3.TNSTalker import TNSClient, TNS_BASE_URL_SANDBOX
from ampel.log.AmpelLogger import AmpelLogger

@pytest.fixture
def test_client():
    if not (api_key := environ.get("TNS_API_KEY")):
        raise pytest.skip("Test requires env var TNS_API_KEY")
    return TNSClient(
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

def test_tnsclient(test_client):
    assert test_client.getInternalName("2018cow") == ('ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa', 'Got internal name response')
    assert test_client.search(244.000917, 22.268031) == (['SN2018cow'], 'Found TNS name(s)')
    assert test_client.getNames(244.000917, 22.268031) == ('SN2018cow', 'ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa')


def test_tnsclient_backoff(test_client):
    import logging
    logging.basicConfig()
    for _ in range(12):
        assert test_client.search(244.000917, 22.268031) == (['SN2018cow'], 'Found TNS name(s)')
