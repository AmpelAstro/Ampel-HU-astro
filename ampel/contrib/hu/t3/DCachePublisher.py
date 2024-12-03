#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/DCachePublisher.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                27.02.2020
# Last Modified Date:  27.02.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>


import datetime
import gzip
import io
import os
import sys
import time
import traceback
from collections.abc import Callable, Generator, Iterable
from functools import partial
from typing import Any, TypeVar
from xml.etree import ElementTree

import backoff
import pytz
import requests

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.struct.T3Store import T3Store
from ampel.util.json import AmpelEncoder
from ampel.view.SnapView import SnapView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id

_T = TypeVar("_T")


class DCachePublisher(AbsT3Unit):
    """
    Publish TransientViews to DCache in gzipped JSON format

    Each run writes a manifest file to {baseDir}/{channel}/manifest/latest.json.gz
    The manifest has the form:
    {
        "time": timestamp,
        "previous": url of previous manifest (or null),
        "updated": list of TransientView URLs
    }

    To build a list of TransientViews updated since some date, read latest.json.gz,
    consume all URLs in "updated", and then repeat with the manifest pointed to
    by "previous," repeating until the "time" of the current manifest is larger
    than the target date.

    Obtain an authorization macaroon with e.g.::

        python -m ampel.contrib.hu.t3.DCachePublisher macaroon -u 'USER:PASS' --validity PT1h
    """

    require = ("dcache",)

    channel: str
    dry_run: bool = False
    base_dir: str = "/ampel/ztf/transient-views"
    max_parallel_requests: int = 8
    authz: NamedSecret[str] = NamedSecret(label="dcache/macaroon")

    def post_init(self) -> None:
        self.updated_urls: list[str] = []
        self.existing_paths: set[tuple[str, ...]] = set()
        self.dt = 0.0

        # don't bother preserving immutable types
        self.encoder = AmpelEncoder(lossy=True)
        assert self.resource
        self.base_dest = self.resource["dcache"] + self.base_dir

    def serialize(self, tran_view: SnapView | dict[str, Any]) -> bytes:
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="w") as f:
            f.write(self.encoder.encode(tran_view).encode("utf-8"))
        return buf.getvalue()

    def create_directory(
        self, session: requests.Session, path_parts: list[str]
    ) -> None:
        path = []
        while len(path_parts) > 0:
            # NB: os.path.join drops any parts before those that start with /
            path.append(path_parts.pop(0).lstrip("/"))
            if tuple(path) in self.existing_paths:
                continue

            if not self.dry_run:
                url = os.path.join(self.base_dest, *path)
                if not self.exists(session, url):
                    resp = session.request("MKCOL", url)
                    # MKCOL on an existing directory returns 'Method not allowed'
                    if not (resp.status_code < 400 or resp.status_code == 405):
                        resp.raise_for_status()
            self.existing_paths.add(tuple(path))

    def put(self, session: requests.Session, url, data, timeout=1.0):
        if self.dry_run:
            return
        session.put(url, data=data)

    def exists(self, session: requests.Session, url: str) -> bool:
        if self.dry_run:
            return False
        resp = session.head(url)
        return resp.status_code in {200, 201, 204}

    @staticmethod
    def is_permanent_error(exc: Exception):
        if isinstance(exc, requests.HTTPError):
            return exc.response.status_code not in {403, 423, 500}
        return isinstance(exc, requests.ConnectionError)

    def _on_backoff(self, details):
        exc_typ, exc, _ = sys.exc_info()
        err = (
            traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
            if exc
            else None
        )
        self.logger.warn(
            "backoff",
            extra={
                "exc": err,
                "target_args": details["args"],
                "wait": details["wait"],
                "tries": details["tries"],
            },
        )

    def _on_giveup(self, details):
        exc_typ, exc, _ = sys.exc_info()
        err = (
            traceback.format_exception_only(exc_typ, exc)[-1].rstrip("\n")
            if exc
            else None
        )
        self.logger.warn(
            "gave up",
            extra={
                "exc": err,
                "target_args": details["args"],
                "tries": details["tries"],
            },
        )

    def publish_transient(self, session: requests.Session, tran_view: SnapView) -> str:
        assert isinstance(tran_view.id, int)
        ztf_name = to_ztf_id(tran_view.id)
        # group such that there are a maximum of 17576 in the same directory
        prefix = [ztf_name[:7], ztf_name[7:9]]

        self.create_directory(session, [self.channel, *prefix])
        base_dir = os.path.join(self.base_dest, self.channel, *prefix)
        fname = base_dir + f"/{ztf_name}.json.gz"
        self.put(session, fname, data=self.serialize(tran_view))

        return fname

    def send_requests(
        self,
        unbound_tasks: Iterable[Callable[[requests.Session], _T]],
    ) -> Iterable[_T]:
        """
        Send a batch of requests

        :param unbound_tasks: sequence of callables that take a single argument
            that is a requests.Session
        """
        session = requests.Session()
        session.headers["Authorization"] = f"BEARER {self.authz.get()}"
        # verify should be a bool or a path to a CA cert dir
        session.verify = os.environ.get("X509_CERT_DIR", True)
        # robustify AsyncClient.request() against transitory failures
        session.request = backoff.on_exception(  # type: ignore[method-assign]
            backoff.expo,
            requests.exceptions.RequestException,
            logger=None,
            giveup=self.is_permanent_error,
            on_giveup=self._on_giveup,
            on_backoff=self._on_backoff,
        )(session.request)

        return (task(session) for task in unbound_tasks)

    def process(
        self, gen: Generator[SnapView, JournalAttributes, None], t3s: T3Store
    ) -> None:
        """
        Publish a transient batch
        """
        t0 = time.time()
        updated_urls = list(
            self.send_requests(
                [
                    partial(self.publish_transient, tran_view=tran_view)
                    for tran_view in gen
                ]
            )
        )
        self.dt += time.time() - t0

        self.logger.info(f"Published {len(updated_urls)} transients in {self.dt:.1f} s")

        if not self.dry_run:
            self.send_requests([self.publish_manifest])

    def publish_manifest(self, session: requests.Session) -> None:
        prefix = os.path.join(self.base_dest, self.channel)
        manifest_dir = os.path.join(prefix, "manifest")
        # create any directories that descend from the base path
        self.create_directory(
            session,
            list(
                os.path.split(
                    prefix[len(os.path.commonprefix((self.base_dest, manifest_dir))) :]
                )
            ),
        )

        latest = os.path.join(manifest_dir, "latest.json.gz")
        previous = None

        # Move the previous manifest aside.
        response = session.request("PROPFIND", latest)
        if (
            response.status_code < 400
            and (
                element := ElementTree.fromstring(response.text).find(
                    "**/{DAV:}prop/{DAV:}creationdate"
                )
            )
            and (date := element.text)
        ):
            previous = os.path.join(manifest_dir, date + ".json.gz")
            self.logger.info(f"Moving {latest} to {previous}")
            (
                session.request(
                    "MOVE",
                    latest,
                    headers={"Destination": previous},
                )
            ).raise_for_status()

        # Write a new manifest with a link to previous manifest.
        (
            session.put(
                latest,
                data=self.serialize(
                    {
                        "time": datetime.datetime.now(tz=pytz.UTC).isoformat(),
                        "previous": previous,
                        "updated": self.updated_urls,
                    }
                ),
            )
        ).raise_for_status()

        # Create and post a read-only macaroon for the prefix path
        response = session.post(
            prefix + "/",
            headers={"Content-Type": "application/macaroon-request"},
            json={
                "caveats": ["activity:DOWNLOAD,LIST,READ_METADATA"],
                "validity": "PT48h",
            },
        )
        authz = response.json()
        session.put(
            os.path.join(prefix, "macaroon"),
            data=authz["macaroon"],
        )
        self.logger.info(f"Access token: {authz}")


def get_macaroon(
    user: requests.auth.HTTPBasicAuth,
    path: str,
    validity: str,
    activity: None | list[str] = None,
    ip: None | list[str] = None,
    host: str = "globe-door.ifh.de",
    port: int = 2880,
):
    """
    Request a DCache macaroon
    """
    caveats = []
    if activity:
        caveats.append(f'activity:{",".join(activity)}')
    if ip:
        caveats.append(f'ip:{",".join(ip)}')
    body: dict[str, Any] = {"validity": validity}
    if caveats:
        body["caveats"] = caveats

    session = requests.Session()
    # verify should be a bool or a path to a CA cert dir
    session.verify = os.environ.get("X509_CERT_DIR", True)
    resp = session.post(
        f"https://{host}:{port}{path}",
        headers={"Content-Type": "application/macaroon-request"},
        json=body,
        auth=user,
    )
    return (resp.json())["macaroon"]


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True

    p = subparsers.add_parser("macaroon")
    p.set_defaults(command=get_macaroon)
    p.add_argument(
        "--path",
        default="/pnfs/ifh.de/acs/ampel/ztf/transient-views/",
        help="path to request macaroon for",
    )
    p.add_argument("--host", default="globe-door.ifh.de")
    p.add_argument("--port", type=int, default=2880)
    p.add_argument(
        "-u",
        "--user",
        type=lambda x: requests.auth.HTTPBasicAuth(*(x.split(":"))),
        help="username:password for basic auth",
    )
    p.add_argument(
        "--validity", default="PT1m", help="token lifetime as ISO 8601 duration"
    )
    p.add_argument(
        "--activity",
        nargs="+",
        choices=(
            "DOWNLOAD",
            "UPLOAD",
            "DELETE",
            "MANAGE",
            "LIST",
            "READ_METADATA",
            "UPDATE_METADATA",
        ),
        help="allowed actions on PATH",
    )
    p.add_argument(
        "--ip", nargs="+", help="netmasks of hosts that may use the macaroon"
    )

    opts = parser.parse_args()
    kwargs = dict(opts.__dict__)
    kwargs.pop("command")
    print(opts.command(**kwargs))  # noqa: T201
