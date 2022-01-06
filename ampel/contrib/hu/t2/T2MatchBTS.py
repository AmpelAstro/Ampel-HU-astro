#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2MatchBts.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                01.06.2021
# Last Modified Date:  13.12.2021
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Union
import requests
import backoff
import pandas as pd
import io
from datetime import datetime

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.view.LightCurve import LightCurve

#from ampel.view.T2DocView import T2DocView
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper
from ampel.enum.DocumentCode import DocumentCode
from bson import tz_util


class T2MatchBTS(AbsLightCurveT2Unit):
    """

    Add information from the BTS explorer page.

    Warning: This will be synced once at init, and the unit output thus depends on
    external sources!
    """

    # Query used to get data from https://sites.astro.caltech.edu/ztf/bts/bts.php
    bts_query: str = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?f=s&format=csv"


    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: e.response.status_code not in {503, 429},
        max_time=60,
    )
    def post_init(self)->None:
        """
        Obtain a synced copy of the BTS explorer output.
        """

        self.bts_df = None

        req = requests.get(self.bts_query)
        req.raise_for_status()

        if req.ok:
            df = pd.read_csv(io.StringIO(req.content.decode()))
            df['synced_at'] = datetime.now(tz_util.utc).timestamp()
            cols = df.columns
            newcols = {col:'bts_'+col for col in cols}
            df.rename(columns=newcols, inplace=True)
            self.bts_df = df



    def process(self, light_curve: LightCurve) -> Union[UBson, UnitResult]:
        """
        Check whether transient exists in the bts df.
        If so, return content

        :returns: dict with content of BTS explorer for transient, together with timestamp.
        E.g.:
            {'bts_ZTFID': 'ZTF21aarqkes',
             'bts_IAUID': 'SN2021hpr',
             'bts_RA': '10:16:38.62',
             'bts_Dec': '+73:24:01.8',
             'bts_peakt': 1324.68,
             'bts_peakfilt': 'g',
             'bts_peakmag': 14.1805,
             'bts_peakabs': '-18.92',
             'bts_duration': '24.527',
             'bts_rise': '11.198',
             'bts_fade': '13.329',
             'bts_type': 'SN Ia',
             'bts_redshift': '0.00934',
             'bts_b': 39.45051306,
             'bts_A_V': 0.067,
             'bts_synced_at': 1622621894.849266}
        """

        # Get name used by BTS
        ztf_name = ZTFIdMapper.to_ext_id(light_curve.stock_id)

        if self.bts_df is None:
            return UnitResult(doc_code=DocumentCode.MISSING_INFO)


        match = self.bts_df[self.bts_df['bts_ZTFID'] == ztf_name].to_dict(orient='index')

        if len(match) == 0:
            # In case of now match, only returned timestamp when check was made
            return {'bts_synced_at' : self.bts_df['bts_synced_at'][0]}

        # Otherwise, return full match dictionary. Assuming unique BTS match, otherwise first entry is retrieved
        return list(match.values())[0]
