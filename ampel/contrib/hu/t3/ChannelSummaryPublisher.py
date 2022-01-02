#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t3/ChannelSummaryPublisher.py
# License:             BSD-3-Clause
# Author:              m. giomi <matteo.giomi@desy.de>
# Date:                13.11.2018
# Last Modified Date:  16.12.2020
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import datetime, backoff, requests
from io import BytesIO, StringIO
from astropy.time import Time
from pytz import timezone
from requests.auth import HTTPBasicAuth
from typing import Any, Dict, Optional, Set, Generator, Union

from ampel.types import ChannelId, UBson, T3Send
from ampel.view.T3Store import T3Store
from ampel.secret.NamedSecret import NamedSecret
from ampel.util.json import AmpelEncoder, load
from ampel.view.TransientView import TransientView
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit


class ChannelSummaryPublisher(AbsPhotoT3Unit):
    """
    Create a json file with summary statistics for the channel. For the transients
    detected in the last N days, this json file contains, i.e. coords, RB score,
    first detection, latest detection, and the total number of transients detected
    by the channel. The summary info for each transient is taken from
    T2LightCurveSummary.
    """

    dry_run: bool = False
    base_url: str = "https://desycloud.desy.de/remote.php/webdav/AMPEL/ZTF"
    auth: NamedSecret[list] = NamedSecret(label="desycloud/valery")


    def post_init(self) -> None:

        self.summary: Dict[str, Any] = {}
        self._jd_range = [float("inf"), -float("inf")]
        self._channels: Set[ChannelId] = set()
        self.session = requests.Session()
        self.session.auth = HTTPBasicAuth(*self.auth.get())


    def extract_from_transient_view(
        self, tran_view: TransientView
    ) -> Optional[Dict[str, Any]]:
        """
        given transient view object return a dictionary
        with the desired metrics
        """
        out: Dict[str, Any] = {}
        assert tran_view.stock
        if names := tran_view.stock.get("name"):
            out["ztf_name"] = next(
                name
                for name in names
                if isinstance(name, str) and name.startswith("ZTF")
            )
            out["tns_names"] = tuple(
                name
                for name in names
                if isinstance(name, str) and name.startswith("TNS")
            )

        # incorporate T2LightCurveSummary
        if summary := tran_view.get_t2_body(unit="T2LightCurveSummary"):
            assert isinstance(summary, dict)
            out.update(summary)
            last_detection = summary["last_detection"]
            if last_detection < self._jd_range[0]:
                self._jd_range[0] = last_detection
            if last_detection > self._jd_range[1]:
                self._jd_range[1] = last_detection

        return out


    def process(self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None) -> Union[UBson, UnitResult]:
        """
        load the stats from the alerts
        """
        for tran_view in gen:
            assert tran_view.stock
            if len(channels := tran_view.stock.get("channel") or []) != 1:
                raise ValueError("Only single-channel views are supported")
            info_dict = self.extract_from_transient_view(tran_view)
            if not info_dict:
                continue
            key = info_dict.pop("ztf_name")
            self.summary[key] = info_dict
            self._channels.add(channels[0])

        self.done()
        return None


    @backoff.on_exception(
        backoff.expo,
        (TimeoutError, requests.exceptions.HTTPError),
        giveup=lambda exc: isinstance(exc, requests.exceptions.HTTPError) and
        exc.response.status_code not in {400, 403, 405, 423, 500}
    )
    def done(self) -> None:
        """"""
        if len(self._channels) == 0:
            return
        elif len(self._channels) > 1:
            raise ValueError(
                f"Got multiple channels ({list(self._channels)}) in summary"
            )

        # Find the date of the most recent observation, in Pacific time
        timestamp = Time(self._jd_range[-1], format="jd").to_datetime(
            timezone("US/Pacific")
        )
        # If before noon, it belongs to the night that started yesterday
        if timestamp.hour < 12:
            timestamp -= datetime.timedelta(days=1)
        filename = timestamp.strftime("channel-summary-%Y%m%d.json")

        channel = list(self._channels)[0]
        basedir = f"{self.base_url}/{channel}"
        rep = self.session.head(basedir)
        if not (rep.ok or self.dry_run):
            self.session.request("MKCOL", basedir).raise_for_status()
        try:
            rep = self.session.get(f"{basedir}/{filename}")
            rep.raise_for_status()
            partial_summary = next(load(BytesIO(rep.content)))
            partial_summary.update(self.summary)
            self.summary = partial_summary
        except (requests.exceptions.HTTPError, StopIteration):
            pass

        outfile = StringIO()
        outfile.write(AmpelEncoder(lossy=True).encode(self.summary))
        outfile.write("\n")
        mb = len(outfile.getvalue().encode()) / 2.0 ** 20
        self.logger.info(
            "{}: {} transients {:.1f} MB".format(filename, len(self.summary), mb)
        )
        if not self.dry_run:
            self.session.put(
                "{}/{}".format(basedir, filename), data=outfile.getvalue()
            ).raise_for_status()
