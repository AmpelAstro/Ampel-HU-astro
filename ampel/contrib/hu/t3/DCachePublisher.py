#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/DCachePublisher.py
# License:             BSD-3-Clause
# Author:              Jakob van Santen <jakob.van.santen@desy.de>
# Date:                27.02.2020
# Last Modified Date:  27.02.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>


import asyncio, datetime, gzip, io, os, ssl, sys, time, traceback, backoff, pytz
from contextlib import nullcontext
from functools import partial
from typing import Any
from collections.abc import Generator
from xml.etree import ElementTree
from aiohttp import ClientSession, helpers, TCPConnector
from aiohttp.client_exceptions import (
    ClientConnectionError,
    ClientConnectorError,
    ClientConnectorSSLError,
    ClientResponseError,
    ServerDisconnectedError,
)

from ampel.types import UBson
from ampel.struct.T3Store import T3Store
from ampel.abstract.AbsT3ReviewUnit import AbsT3ReviewUnit
from ampel.secret.NamedSecret import NamedSecret
from ampel.util.json import AmpelEncoder
from ampel.view.SnapView import SnapView
from ampel.ztf.util.ZTFIdMapper import to_ztf_id
from ampel.struct.UnitResult import UnitResult
from ampel.struct.JournalAttributes import JournalAttributes


class DCachePublisher(AbsT3ReviewUnit):
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

    async def create_directory(self, request, path_parts) -> None:
        path = []
        while len(path_parts) > 0:
            # NB: os.path.join drops any parts before those that start with /
            path.append(path_parts.pop(0).lstrip("/"))
            if tuple(path) in self.existing_paths:
                continue

            if not self.dry_run:
                url = os.path.join(self.base_dest, *path)
                if not await self.exists(request, url):
                    resp = await request("MKCOL", url, raise_for_status=False)
                    # MKCOL on an existing directory returns 'Method not allowed'
                    if not (resp.status < 400 or resp.status == 405):
                        resp.raise_for_status()
            self.existing_paths.add(tuple(path))

    async def put(self, request, url, data, timeout=1.0):
        if self.dry_run:
            return
        await request("PUT", url, data=data)

    async def exists(self, request, url):
        if self.dry_run:
            return
        resp = await request("HEAD", url, raise_for_status=False)
        return resp.status in {200, 201, 204}

    @staticmethod
    def is_permanent_error(exc):
        if isinstance(exc, ClientResponseError):
            return exc.code not in {403, 423, 500}
        elif isinstance(exc, ClientConnectorSSLError):
            return True
        else:
            return False

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

    async def publish_transient(self, request, tran_view: SnapView) -> str:
        assert isinstance(tran_view.id, int)
        ztf_name = to_ztf_id(tran_view.id)
        # group such that there are a maximum of 17576 in the same directory
        prefix = [ztf_name[:7], ztf_name[7:9]]

        await self.create_directory(request, [self.channel, *prefix])
        base_dir = os.path.join(self.base_dest, self.channel, *prefix)
        fname = base_dir + f"/{ztf_name}.json.gz"
        await self.put(request, fname, data=self.serialize(tran_view))

        return fname

    async def send_requests(self, unbound_tasks):
        """
        Send a batch of requests

        :param unbound_tasks: sequence of callables that take a single argument
            that is an aiohttp request object
        """
        ssl_context = ssl.create_default_context(
            capath=os.environ.get("X509_CERT_DIR", None)
        )
        async with TCPConnector(
            limit=self.max_parallel_requests, ssl_context=ssl_context
        ) as connector:
            async with ClientSession(
                connector=connector,
                headers={"Authorization": f"BEARER {self.authz.get()}"},
                raise_for_status=True,
            ) as session:
                # ClientSession.request is a normal function returning an
                # async context manager. While this is kind of equivalent
                # to a coroutine, it can't be decorated with backoff. Wrap
                # it so it can be.
                # Note that it is important to not actually enter the
                # context manager here, as this would close the response
                # on exit, causing read() to raise ClientConnectionError.
                async def async_request(*args, **kwargs):
                    return await session.request(*args, **kwargs)

                # robustify Session.request() against transitory failures
                request = backoff.on_exception(
                    backoff.expo,
                    (
                        asyncio.TimeoutError,
                        ClientResponseError,
                        ClientConnectorError,
                        ClientConnectionError,
                        ServerDisconnectedError,
                    ),
                    logger=None,
                    giveup=self.is_permanent_error,
                    on_giveup=self._on_giveup,
                    on_backoff=self._on_backoff,
                )(async_request)

                tasks = [task(request) for task in unbound_tasks]
                return await asyncio.gather(*tasks)

    def process(
        self, gen: Generator[SnapView, JournalAttributes, None], t3s: T3Store
    ) -> UBson | UnitResult:
        """
        Publish a transient batch
        """
        loop = asyncio.get_event_loop()
        t0 = time.time()
        self.updated_urls += loop.run_until_complete(
            self.send_requests(
                [
                    partial(self.publish_transient, tran_view=tran_view)
                    for tran_view in gen
                ]
            )
        )
        self.dt += time.time() - t0

        self.logger.info(
            f"Published {len(self.updated_urls)} transients in {self.dt:.1f} s"
        )

        if not self.dry_run:
            asyncio.get_event_loop().run_until_complete(
                self.send_requests([self.publish_manifest])
            )

        return None

    async def publish_manifest(self, request) -> None:
        prefix = os.path.join(self.base_dest, self.channel)
        manifest_dir = os.path.join(prefix, "manifest")
        # create any directories that descend from the base path
        await self.create_directory(
            request,
            list(
                os.path.split(
                    prefix[len(os.path.commonprefix((self.base_dest, manifest_dir))) :]
                )
            ),
        )

        latest = os.path.join(manifest_dir, "latest.json.gz")
        previous = None

        # Move the previous manifest aside.
        async with await request(
            "PROPFIND", latest, raise_for_status=False
        ) as response:
            if response.status < 400 and (
                element := ElementTree.fromstring(await response.text()).find(
                    "**/{DAV:}prop/{DAV:}creationdate"
                )
            ):
                if date := element.text:
                    previous = os.path.join(manifest_dir, date + ".json.gz")
                    self.logger.info(f"Moving {latest} to {previous}")
                    await request(
                        "MOVE",
                        latest,
                        headers={"Destination": previous},
                    )

        # Write a new manifest with a link to previous manifest.
        await request(
            "PUT",
            latest,
            data=self.serialize(
                {
                    "time": datetime.datetime.now(tz=pytz.UTC).isoformat(),
                    "previous": previous,
                    "updated": self.updated_urls,
                }
            ),
        )

        # Create and post a read-only macaroon for the prefix path
        async with await request(
            "POST",
            prefix + "/",
            headers={"Content-Type": "application/macaroon-request"},
            json={
                "caveats": ["activity:DOWNLOAD,LIST,READ_METADATA"],
                "validity": "PT48h",
            },
        ) as response:
            authz = await response.json()
        await request(
            "PUT",
            os.path.join(prefix, "macaroon"),
            data=authz["macaroon"],
        )
        self.logger.info(f"Access token: {authz}")


def get_macaroon(
    user: helpers.BasicAuth,
    path: str,
    validity: str,
    activity: list[str] = [],
    ip: list[str] = [],
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

    async def fetch():
        ssl_context = ssl.create_default_context(
            capath=os.environ.get("X509_CERT_DIR", None)
        )
        async with TCPConnector(ssl_context=ssl_context) as connector:
            async with ClientSession(auth=user, connector=connector) as session:
                resp = await session.post(
                    f"https://{host}:{port}{path}",
                    headers={"Content-Type": "application/macaroon-request"},
                    json=body,
                )
                return (await resp.json())["macaroon"]

    return asyncio.get_event_loop().run_until_complete(fetch())


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
        type=lambda x: helpers.BasicAuth(*(x.split(":"))),
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
    print(opts.command(**kwargs))
