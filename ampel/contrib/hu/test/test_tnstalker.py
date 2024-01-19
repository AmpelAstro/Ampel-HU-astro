from os import environ

import pytest

from ampel.contrib.hu.t3.tns.TNSToken import TNSToken
from ampel.contrib.hu.t3.TNSTalker import TNS_BASE_URL_SANDBOX, TNSClient
from ampel.log.AmpelLogger import AmpelLogger


@pytest.fixture
def tns_token():
    if not (api_key := environ.get("TNS_API_KEY")):
        raise pytest.skip("Test requires env var TNS_API_KEY")
    return TNSToken(
        id=59228,
        name="ZTF_AMPEL_COMPLETE",
        api_key=api_key,
    )


@pytest.fixture
def test_client(tns_token):
    return TNSClient(TNS_BASE_URL_SANDBOX, AmpelLogger.get_logger(), tns_token)


def test_tnsclient(test_client):
    assert test_client.getInternalName("2018cow") == (
        "ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa",
        "Got internal name response",
    )
    assert test_client.search(244.000917, 22.268031) == (
        ["SN2018cow"],
        "Found TNS name(s)",
    )
    assert test_client.getNames(244.000917, 22.268031) == (
        "SN2018cow",
        "ATLAS18qqn, ZTF18abcfcoo, Gaia18bqa",
    )


def test_tnsclient_backoff(test_client: TNSClient):
    import logging

    logging.basicConfig()
    for _ in range(12):
        assert test_client.search(244.000917, 22.268031) == (
            ["SN2018cow"],
            "Found TNS name(s)",
        )


@pytest.mark.asyncio
async def test_tnsclient_backoff_async(tns_token):
    from ampel.contrib.hu.t3.tns import TNSClient as TNSMirrorClient

    client = TNSMirrorClient(
        tns_token, timeout=120, maxParallelRequests=1, logger=AmpelLogger.get_logger()
    )
    ra, dec, matchradius = 244.000917, 22.268031, 5.0
    for _ in range(12):
        hits = [
            hit
            async for hit in client.search(
                **{"ra": ra, "dec": dec, "radius": matchradius, "units": "arcsec"}
            )
        ]
        assert hits
