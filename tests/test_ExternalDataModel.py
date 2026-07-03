import secrets
from base64 import b64encode
from hashlib import md5
from unittest.mock import MagicMock

from pytest_mock import MockerFixture
from requests import Response

from ampel.model.ExternalDataModel import ExternalDataModel


def test_download(tmpdir, mocker: MockerFixture):

    content = secrets.token_bytes(36)
    md5_hash = md5(content).hexdigest()

    # do not actually call out to httpbin
    response = Response()
    response.raise_for_status = lambda: None
    response.iter_content = lambda chunk_size: [content]
    response.raw = MagicMock()
    get = mocker.patch(
        "ampel.model.ExternalDataModel.requests.get", return_value=response
    )

    # NB: use urlsafe altchars, the same encoding as httpbin assumes
    url = f"https://httpbin.org/base64/{b64encode(content, altchars=b'-_').decode()}"

    m = ExternalDataModel(
        name=str(tmpdir / "test.txt"),
        url=url,
        md5=md5_hash,
    )
    path = m.local_path()
    assert get.call_count == 1
    assert path.exists()

    assert path.read_bytes() == content

    m._cache.clear()
    path = m.local_path()
    assert path.exists()
    assert get.call_count == 1  # no new download

    path.unlink()
    assert not path.exists()

    m.local_path()
    assert path.read_bytes() == content
    assert get.call_count == 2  # new download after file was deleted
