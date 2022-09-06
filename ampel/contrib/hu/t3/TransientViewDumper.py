#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/TransientInfoDumper.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                15.08.2018
# Last Modified Date:  15.08.2018
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import uuid, requests
from requests.auth import HTTPBasicAuth
from gzip import GzipFile
from io import BytesIO
from collections.abc import Generator
from urllib.parse import ParseResult, urlparse, urlunparse
from xml.etree import ElementTree
from ampel.types import UBson, T3Send
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit
from ampel.secret.NamedSecret import NamedSecret
from ampel.util.json import AmpelEncoder
from ampel.view.SnapView import SnapView
from ampel.view.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult


def strip_auth_from_url(url):
    try:
        auth = requests.utils.get_auth_from_url(url)
        scheme, netloc, path, params, query, fragment = urlparse(url)
        netloc = netloc[(netloc.index("@") + 1):]
        url = urlunparse(ParseResult(scheme, netloc, path, params, query, fragment))
        return url, auth
    except KeyError:
        return url, None


def strip_path_from_url(url):
    scheme, netloc, path, params, query, fragment = urlparse(url)
    return urlunparse(ParseResult(scheme, netloc, "/", None, None, None))


class TransientViewDumper(AbsT3ReviewUnit):
    """"""

    version = 0.1
    resources = ("desycloud",)

    outputfile: None | str = None
    desycloud_auth: NamedSecret[dict] = NamedSecret(label="desycloud")

    def post_init(self) -> None:
        if not self.outputfile:
            self.outfile = GzipFile(fileobj=BytesIO(), mode="w")
            self.path = "/AMPEL/dumps/" + str(uuid.uuid1()) + ".json.gz"
            self.session = requests.Session()
            assert self.resource
            self.webdav_base = self.resource["desycloud"]
            self.ocs_base = (
                strip_path_from_url(self.resource["desycloud"]) +
                "/ocs/v1.php/apps/files_sharing/api/v1"
            )
        else:
            self.outfile = GzipFile(self.outputfile + ".json.gz", mode="w")
        # don't bother preserving immutable types
        self.encoder = AmpelEncoder(lossy=True)


    def process(self, transients: Generator[SnapView, T3Send, None], t3s: T3Store) -> UBson | UnitResult:

        count = 0
        for count, tran_view in enumerate(transients, 1):
            self.outfile.write(self.encoder.encode(tran_view).encode("utf-8"))
            self.outfile.write(b"\n")

        self.outfile.flush()
        self.logger.info("Total number of transient printed: %i" % count)
        if self.outputfile:
            self.outfile.close()
            self.logger.info(self.outputfile + ".json.gz")
        else:
            assert isinstance(self.outfile.fileobj, BytesIO)
            mb = len(self.outfile.fileobj.getvalue()) / 2.0 ** 20
            self.logger.info("{:.1f} MB of gzipped JSONy goodness".format(mb))
            auth = HTTPBasicAuth(**self.desycloud_auth.get())
            self.session.put(
                self.webdav_base + self.path,
                data=self.outfile.fileobj.getvalue(),
                auth=auth,
            ).raise_for_status()
            response = self.session.post(
                self.ocs_base + "/shares",
                data=dict(path=self.path, shareType=3),
                auth=auth,
                headers={"OCS-APIRequest": "true"},  # i'm not a CSRF attack, i swear
            )
            self.outfile.close()
            if response.ok and (
                element := ElementTree.fromstring(response.text).find("data/url")
            ):
                if element.text:
                    self.logger.info(element.text)
            else:
                response.raise_for_status()
        return None
