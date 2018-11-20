
import pytest

@pytest.mark.xfail(reason="Credentials need to be set as env vars DESYCLOUD_USERNAME and DESYCLOUD_PASSWORD")
def test_desycloud():
    import requests
    from ampel.pipeline.config.AmpelArgumentParser import AmpelArgumentParser
    from ampel.pipeline.config.AmpelConfig import AmpelConfig

    parser = AmpelArgumentParser()
    parser.require_resource('desycloud')
    parser.parse_args(args=[])

    uri = AmpelConfig.get_config('resources.desycloud')
    dest = uri+'/AMPEL/helloworld.txt'
    content = 'hi there'

    # Upload, read, and delete a file with PUT/GET/DELETE. WebDAV is pretty
    # simple as long as you don't actually use any of its features.
    try:
        requests.put(dest, data=content).raise_for_status()

        rep = requests.get(dest)
        rep.raise_for_status()
        assert rep.text == content

    finally:
        rep = requests.delete(dest)
    rep.raise_for_status()

    assert not requests.get(dest).ok
