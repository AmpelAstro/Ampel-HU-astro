#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TransientInfoDumper.py
# License           : BSD-3-Clause
# Author            : Jakob van Santen <jakob.van.santen@desy.de>
# Date              : 15.08.2018
# Last Modified Date: 02.05.2023
# Last Modified By  : Simeon Reusch <simeon.reusch@desy.de>
import uuid
from collections.abc import Generator
from gzip import GzipFile
from io import BytesIO
from urllib.parse import ParseResult, urlparse, urlunparse
from xml.etree import ElementTree

import requests
from requests.auth import HTTPBasicAuth

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import T3Send, UBson
from ampel.util.json import AmpelEncoder
from ampel.view.TransientView import TransientView


def strip_auth_from_url(url):
    try:
        auth = requests.utils.get_auth_from_url(url)
        scheme, netloc, path, params, query, fragment = urlparse(url)
        netloc = netloc[(netloc.index("@") + 1) :]
        url = urlunparse(ParseResult(scheme, netloc, path, params, query, fragment))
        return url, auth
    except KeyError:
        return url, None


def strip_path_from_url(url):
    scheme, netloc, *_ = urlparse(url)
    return urlunparse(ParseResult(scheme, netloc, "/", "", "", ""))


class TransientViewDumper(AbsPhotoT3Unit):
    """"""

    version = 0.1
    require = ("desycloud",)

    # If this is passed, files are always saved locally
    outputfile: None | str = None

    desycloud_auth: NamedSecret[dict] = NamedSecret[dict](label="desycloud")
    desycloud_folder: str = "dumps"
    desycloud_filename: str = str(uuid.uuid1())

    def post_init(self) -> None:
        if not self.outputfile:
            self.buffer = BytesIO()
            self.outfile = GzipFile(
                filename=self.desycloud_filename + ".json",
                fileobj=self.buffer,
                mode="w",
            )
            self.path = (
                f"/AMPEL/{self.desycloud_folder}/"
                + self.desycloud_filename
                + ".json.gz"
            )
            self.session = requests.Session()
            assert self.resource
            self.webdav_base = self.resource["desycloud"]
            self.ocs_base = (
                strip_path_from_url(self.resource["desycloud"])
                + "/ocs/v1.php/apps/files_sharing/api/v1"
            )
        else:
            self.outfile = GzipFile(self.outputfile + ".json.gz", mode="w")

        # don't bother preserving immutable types
        self.encoder = AmpelEncoder(lossy=True)

    def process(
        self, transients: Generator[TransientView, T3Send, None], t3s: T3Store
    ) -> UBson | UnitResult:
        count = 0
        for count, tran_view in enumerate(transients, 1):  # noqa: B007
            self.outfile.write(self.encoder.encode(tran_view).encode("utf-8"))
            self.outfile.write(b"\n")
        self.outfile.close()
        self.logger.info(f"Total number of transients written: {count}")
        if self.outputfile:
            self.logger.info(self.outputfile + ".json.gz")
        else:
            assert isinstance(self.buffer, BytesIO)
            mb = len(self.buffer.getvalue()) / 2.0**20
            self.logger.info(f"{mb:.1f} MB of gzipped JSONy goodness")
            auth = HTTPBasicAuth(**self.desycloud_auth.get())

            self.session.put(
                self.webdav_base + self.path,
                data=self.buffer.getvalue(),
                auth=auth,
            ).raise_for_status()
            response = self.session.post(
                self.ocs_base + "/shares",
                data=dict(path=self.path, shareType=3),
                auth=auth,
                headers={"OCS-APIRequest": "true"},  # i'm not a CSRF attack, i swear
            )
            if response.ok and (
                element := ElementTree.fromstring(response.text).find("data/url")
            ):
                if element.text:
                    self.logger.info(element.text)
            else:
                response.raise_for_status()
        return None
