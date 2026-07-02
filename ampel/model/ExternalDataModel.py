import os
from functools import cache, lru_cache
from hashlib import md5 as _md5
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, ClassVar
from urllib.parse import urlparse

import requests
from pydantic import Field

from ampel.base.AmpelBaseModel import AmpelBaseModel

if TYPE_CHECKING:
    from glide_sync import GlideClient


class ExternalDataModel(AmpelBaseModel):
    name: Annotated[str, Field(description="Name of the file")]
    url: Annotated[str, Field(description="URL of the file")]
    md5: Annotated[
        str, Field(description="MD5 hash of the file", pattern=r"^[a-fA-F0-9]{32}$")
    ]

    _cache: ClassVar[dict[tuple[str, str, str], Path]] = {}

    @staticmethod
    @cache
    def _get_valkey() -> "GlideClient | None":
        url = os.getenv("VALKEY_URL")
        if not url:
            return None
        parts = urlparse(url)
        if (hostname := parts.hostname) is None:
            return None

        try:
            from glide_sync import (  # noqa: PLC0415
                GlideClient,
                GlideClientConfiguration,
                NodeAddress,
            )
        except ImportError:
            return None

        return GlideClient.create(
            GlideClientConfiguration(
                [NodeAddress(host=hostname, port=parts.port or 6379)]
            )
        )

    @staticmethod
    @lru_cache(maxsize=1024)
    def _file_md5(path: Path) -> str:
        h = _md5()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    def local_path(self) -> Path:
        """
        Get the local path to the file, downloading it if necessary."""
        key = (self.name, self.url, self.md5)
        cached = self._cache.get(key)
        if cached and cached.exists() and self._file_md5(cached) == self.md5:
            return cached

        path = Path(self.name)
        if path.exists():
            if self._file_md5(path) == self.md5:
                self._cache[key] = path
                return path
            raise ValueError(
                f"{path} exists, but MD5 does not match expected value for {self.url}. "
                "Remove the file or change the name to download a new version."
            )

        path.parent.mkdir(parents=True, exist_ok=True)
        valkey = self._get_valkey()
        if valkey and (content := valkey.get(self.url)):
            with path.open("wb") as f:
                f.write(content)
        else:
            with (
                requests.get(self.url, stream=True) as response,
                path.open("wb") as out,
            ):
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        out.write(chunk)
            if valkey:
                # NB: values must be smaller than 512MB
                valkey.set(self.url, path.read_bytes())

        if self._file_md5(path) != self.md5:
            path.unlink(missing_ok=True)
            if valkey:
                valkey.delete([self.url])
            raise ValueError(
                f"MD5 mismatch for {self.url}: expected {self.md5}, got {self._file_md5(path)}"
            )

        self._cache[key] = path
        return path
