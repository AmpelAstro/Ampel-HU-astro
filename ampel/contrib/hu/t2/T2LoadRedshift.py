#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2LoadRedshift.py
# License:             BSD-3-Clause
# Author:              alice.townsend@physik.hu-berlin.de
# Date:                06.04.2023
# Last Modified Date:  06.04.2023
# Last Modified By:    alice.townsend@physik.hu-berlin.de

from datetime import datetime

import pandas as pd
from bson import tz_util

# from ampel.view.T2DocView import T2DocView
from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.LightCurve import LightCurve
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper


class T2LoadRedshift(AbsLightCurveT2Unit):
    """

    Add redshifts from external .csv
    """

    # Path to file
    df_path: str = "/Users/alicetownsend/FPlist_all_analysis/desi/pandas_db2.csv"

    def post_init(self) -> None:
        """
        Obtain a recent copy of df.
        """

        self.z_df = None

        df = pd.read_csv(self.df_path)
        df["synced_at"] = datetime.now(tz_util.utc).timestamp()
        cols = df.columns
        newcols = {col: "T2LoadRedshift_" + col for col in cols}
        df.rename(columns=newcols, inplace=True)
        self.z_df = df

    def process(self, light_curve: LightCurve) -> UBson | UnitResult:
        """
        Check whether transient exists in the df.
        If so, return content
        """

        # Get name used by BTS
        assert isinstance(light_curve.stock_id, int)
        ztf_name = ZTFIdMapper.to_ext_id(light_curve.stock_id)

        if self.z_df is None:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        match = next(
            iter(
                self.z_df[self.z_df["T2LoadRedshift_ztfid"] == ztf_name]
                .to_dict(orient="index")
                .values()
            ),
            None,
        )
        if match is None:
            # In case of no match, only returned timestamp when check was made
            return {
                "T2LoadRedshift_synced_at": self.z_df["T2LoadRedshift_synced_at"][0]
            }
        # Otherwise, return full match dictionary. Assuming unique BTS match, otherwise first entry is retrieved
        return {str(k): v for k, v in match.items()}
