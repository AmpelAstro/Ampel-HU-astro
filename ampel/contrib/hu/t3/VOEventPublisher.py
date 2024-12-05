#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t3/VOEventPublisher.py
# License:             BSD-3-Clause
# Author:              jno <jnordin@physik.hu-berlin.de>
# Date:                26.09.2022
# Last Modified Date:  26.09.2022
# Last Modified By:    jno <jnordin@physik.hu-berlin.de>

import datetime
from collections.abc import Generator
from typing import TYPE_CHECKING, Any, Literal

import voeventparse as vp
from astropy.time import Time

from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.struct.T3Store import T3Store
from ampel.types import T3Send
from ampel.util.mappings import get_by_path
from ampel.view.TransientView import TransientView

if TYPE_CHECKING:
    pass


class VOEventPublisher(AbsPhotoT3Unit):
    """

    Unit for creating a VOEvent pased on T2output and transient LightCurve.

    Selection of fields to save. Matches structure of t2document result dict, e.g.:
    why_schema = { { 't2_unit_name'  : {
                'voevent_label' : ['path','to','val'],
                'voevent_label_2' : ['other','path']
            },
        } }

    """

    # Characteristics of the VOEvent
    voevent_ivorn: str = "AMPEL/dev"
    voevent_role: str = "test"
    voevent_stream: str = "myFastDecliners"
    vovent_streamid: int = 1

    # Selection of fields to output ("WHY?")
    name_filter: dict[str, str] = {"ZTF name": "ZTF", "TNS ID": "TNS"}
    # Which datapoints to include in submission?
    which_photometry: Literal["first", "last"] = "first"
    # Schema for state dependent T2s (one row for each)
    why_schema: dict[str, Any]
    # Temporary file name
    fname = "voevent.xml"

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: None | T3Store = None
    ) -> None:
        """
        Loop through provided TransientViews and extract data according to the
        configured schema.
        """

        for k, tran_view in enumerate(gen, 1):
            # Chategorized by stock (internal ampel ID) and channel
            assert tran_view.stock is not None
            channels = tran_view.stock.get("channel")
            assert channels is not None
            channel = (
                channels[0]
                if isinstance(channels, list | tuple) and len(channels) == 1
                else "/".join([str(c) for c in channels])
            )

            # Get photometry from first|last visit
            dps = tran_view.get_photopoints()
            assert dps is not None
            dps = sorted(dps, key=lambda d: d["body"]["jd"])
            if self.which_photometry == "first":
                dp = dps[0]["body"]
            elif self.which_photometry == "last":
                dp = dps[-1]["body"]

            # Collect t2 information.
            # Only uses latest state. See TransientTablePublisher for
            # possible extension to state.
            t2dict: dict[str, Any] = {}
            for t2unit, table_entries in self.why_schema.items():
                for t2res in tran_view.get_t2_views(unit=t2unit):
                    for label, path in table_entries.items():
                        assert t2res.body
                        body = t2res.body[-1]
                        assert isinstance(body, dict)
                        if result := get_by_path(body, path):
                            t2dict[label] = result
            if len(t2dict.keys()) == 0:
                continue

            # Generate VOEvent
            v = vp.Voevent(
                stream=self.voevent_stream,
                stream_id=self.vovent_streamid,
                role=self.voevent_role,
            )
            vp.set_who(v, datetime.datetime.now(tz=datetime.timezone.utc))
            vp.set_author(
                v,
                title=f"Results from {channel} channel",
                contactName="ampel-info@desy.de",
            )
            v.Description = "Generated through https://github.com/AmpelProject/Ampel-HU-astro/ampel/contrib/hu/t3/VOEventPublisher.py"

            v.What.append(vp.Param(name="mag", value=dp["magpsf"], ucd="phot.mag"))

            dt = Time(dp["jd"], format="jd").to_datetime()
            dt = dt.replace(tzinfo=datetime.timezone.utc)

            vp.add_where_when(
                v,
                coords=vp.Position2D(
                    ra=dp["ra"],
                    dec=dp["dec"],
                    err=0,
                    units="deg",
                    system=vp.definitions.sky_coord_system.utc_fk5_geo,
                ),
                obs_time=dt,
                observatory_location="Palomar P48 / ZTF camera",
            )

            # Add collected results
            vp.add_how(v, descriptions=[f"{k}: {v}" for k, v in t2dict.items()])
            vp.add_why(v)
            v.Why.Description = "Selected based on AMPEL T2s: {}".format(
                " ".join(self.why_schema.keys())
            )

            with open(self.fname, "wb") as fw:
                vp.dump(v, fw)
