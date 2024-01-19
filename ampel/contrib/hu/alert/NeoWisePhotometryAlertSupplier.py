#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-ZTF/ampel/ztf/alert/NeoWisePhotometryAlertSupplier.py
# License           : BSD-3-Clause
# Author            :
# Date              :
# Last Modified Date:
# Last Modified By  :

from typing import Any, Literal, Dict, List, Optional, Sequence, Any, Tuple
import sys, os
from bson import encode
from hashlib import blake2b
from ampel.ztf.util.ZTFIdMapper import to_ampel_id
from ampel.view.ReadOnlyDict import ReadOnlyDict
from ampel.alert.BaseAlertSupplier import BaseAlertSupplier
from ampel.alert.AmpelAlert import AmpelAlert
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.types import Tag
import numpy as np
import pandas as pd
import astropy
import json
from astropy.io import ascii
from io import BytesIO


class NeoWisePhotometryAlertSupplier(BaseAlertSupplier):
    """
    Iterable class that, for each transient name  provided by the underlying alert_loader
    returns a PhotoAlert instance.

    Assumes file format according to that provided by J Necker.

    Some steps known to have been applied are to have a ref flux subtracted and observations within each 6-month visit
    averaged such that the std can be used as error.
    """

    stat_pps: int = 0
    stat_uls: int = 0
    #    force_pos_flux : bool = False

    dpid: Literal["hash", "inc"] = "hash"
    #    external_directory: Optional[ str ]
    #    deserialize: None | Literal["avro", "json"]

    def __init__(self, **kwargs) -> None:
        kwargs["deserialize"] = "json"
        super().__init__(**kwargs)
        self.counter = 0 if self.dpid == "hash" else 1

    def __next__(self) -> AmpelAlert:
        """
        :returns: a dict with a structure that AlertProcessor understands
        :raises StopIteration: when alert_loader dries out.
        :raises AttributeError: if alert_loader was not set properly before this method is called
        """
        d = self._deserialize(next(self.alert_loader))  # type: ignore

        # assure that a timewise lightcurve is in the data
        while "timewise_lightcurve" not in d[1]:
            d = self._deserialize(next(self.alert_loader))

        transient_name = d[0]
        df = pd.DataFrame.from_dict(d[1]["timewise_lightcurve"])

        # Some units assume jd, and conversion not yet completed
        df["jd"] = df["mean_mjd"] + 2400000.5
        df["W1_mag_ul"].replace(0.0, "False", inplace=True)
        df["W2_mag_ul"].replace(0.0, "False", inplace=True)
        #        df["W1_flux_ul"].replace(0.0, "False", inplace=True)
        #        df["W2_flux_ul"].replace(0.0, "False", inplace=True)
        df["W1_flux_density_ul"].replace(0.0, "False", inplace=True)
        df["W2_flux_density_ul"].replace(0.0, "False", inplace=True)
        df["W1_mag_ul"].replace(1.0, "True", inplace=True)
        df["W2_mag_ul"].replace(1.0, "True", inplace=True)
        df["W1_flux_density_ul"].replace(1.0, "True", inplace=True)
        df["W2_flux_density_ul"].replace(1.0, "True", inplace=True)

        if "timewise_metadata" in d[1].keys():
            # calculate reduced chi2
            timewise_metadata = d[1]["timewise_metadata"]
            for b in ["W1", "W2"]:
                timewise_metadata[f"{b}_red_chi2"] = (
                    timewise_metadata[f"{b}_chi2_to_med_flux_density"]
                    / (timewise_metadata[f"{b}_N_datapoints_flux_density"] - 1)
                    if timewise_metadata[f"{b}_N_datapoints_flux_density"] > 1
                    else np.nan
                )

            df[list(timewise_metadata.keys())] = pd.DataFrame(
                list([timewise_metadata.values()]), index=df.index
            )
            selected_columns_W1 = [
                "mean_mjd",
                "W1_mean_mag",
                "W1_mag_rms",
                "W1_mag_ul",
                "W1_mean_flux_density",
                "W1_flux_density_rms",
                "W1_flux_density_ul",
                "jd",
            ] + [col for col in list(timewise_metadata) if "W1" in col]
            selected_columns_W2 = [
                "mean_mjd",
                "W2_mean_mag",
                "W2_mag_rms",
                "W2_mag_ul",
                "W2_mean_flux_density",
                "W2_flux_density_rms",
                "W2_flux_density_ul",
                "jd",
            ] + [col for col in list(timewise_metadata) if "W2" in col]
        else:
            selected_columns_W1 = [
                "mean_mjd",
                "W1_mean_mag",
                "W1_mag_rms",
                "W1_mag_ul",
                "W1_mean_flux_density",
                "W1_flux_density_rms",
                "W1_flux_density_ul",
                "jd",
            ]
            selected_columns_W2 = [
                "mean_mjd",
                "W2_mean_mag",
                "W2_mag_rms",
                "W2_mag_ul",
                "W2_mean_flux_density",
                "W2_flux_density_rms",
                "W2_flux_density_ul",
                "jd",
            ]
        df_W1 = df[selected_columns_W1].copy()
        df_W2 = df[selected_columns_W2].copy()

        df_W1.rename(
            columns={
                "W1_mean_flux_density": "mean_flux",
                "W1_flux_density_rms": "flux_rms",
                "W1_mean_mag": "mean_mag",
                "W1_mag_rms": "mag_rms",
                "W1_mag_ul": "mag_ul",
                "W1_flux_density_ul": "flux_ul",
            },
            inplace=True,
        )
        df_W2.rename(
            columns={
                "W2_mean_flux_density": "mean_flux",
                "W2_flux_density_rms": "flux_rms",
                "W2_mean_mag": "mean_mag",
                "W2_mag_rms": "mag_rms",
                "W2_mag_ul": "mag_ul",
                "W2_flux_density_ul": "flux_ul",
            },
            inplace=True,
        )

        df_W1.columns = df_W1.columns.str.replace("W1_", "")
        df_W2.columns = df_W2.columns.str.replace("W2_", "")

        if "W1_mag_Npoints" in df.columns:
            df_W1["mag_Npoints"] = df["W1_mag_Npoints"]
            df_W2["mag_Npoints"] = df["W2_mag_Npoints"]
            df_W1["flux_density_Npoints"] = df["W1_flux_density_Npoints"]
            df_W2["flux_density_Npoints"] = df["W2_flux_density_Npoints"]
        if "ra" in d[1].keys():
            df_W1["ra"] = d[1]["ra"]
            df_W1["dec"] = d[1]["dec"]
            df_W2["ra"] = d[1]["ra"]
            df_W2["dec"] = d[1]["dec"]

        df_W1 = df_W1.astype({"flux_ul": str, "mag_ul": str})
        df_W2 = df_W2.astype({"flux_ul": str, "mag_ul": str})

        # Define the standard fields
        df_W1["filter"] = "Wise_W1"
        #        df_W1['zp']  = self.zp_W1
        df_W1["magpsf"] = df_W1["mean_mag"]
        df_W1["sigmapsf"] = df_W1["mag_rms"]
        df_W1["programid"] = 1  #  Faking this to use the standard ingester
        df_W1["fid"] = 1  #  Faking this to use the standard ingester
        ipos = df_W1["magpsf"] < 999
        df_W1["isdiffpos"] = "f"
        df_W1.loc[ipos, "isdiffpos"] = "t"
        #        df_W1['flux_ul'] = str(df_W1['flux_ul'])

        df_W2["filter"] = "Wise_W2"
        #        df_W2['zp']  = self.zp_W2
        df_W2["magpsf"] = df_W2["mean_mag"]
        df_W2["sigmapsf"] = df_W2["mag_rms"]
        df_W2["programid"] = 1  #  Faking this to use the standard ingester
        df_W2["fid"] = 2  #  Faking this to use the standard ingester
        ipos = df_W2["magpsf"] < 999
        df_W2["isdiffpos"] = "f"
        df_W2.loc[ipos, "isdiffpos"] = "t"
        #        df_W2['flux_ul'] = str(df_W2['flux_ul'])

        # Ingester makes use of the readout quadrant ID when looking for superceded data.
        # Will arbitrary set this according to the different filters, as these are often set to the same date
        df_W1["rcid"] = 1
        df_W2["rcid"] = 2

        all_ids = b""
        pps = []
        for index, row in df_W1.iterrows():
            pp = dict(row)
            pp_hash = blake2b(encode(pp), digest_size=7).digest()
            if self.counter:
                pp["candid"] = self.counter
                self.counter += 1
            else:
                pp["candid"] = int.from_bytes(pp_hash, byteorder=sys.byteorder)

            all_ids += pp_hash
            pps.append(ReadOnlyDict(pp))
        for index, row in df_W2.iterrows():
            pp = dict(row)
            pp_hash = blake2b(encode(pp), digest_size=7).digest()
            if self.counter:
                pp["candid"] = self.counter
                self.counter += 1
            else:
                pp["candid"] = int.from_bytes(pp_hash, byteorder=sys.byteorder)

            all_ids += pp_hash
            pps.append(ReadOnlyDict(pp))

        if not pps:
            return self.__next__()

        # Update stats
        #        self.stat_pps += len(pps)
        ##        t = tuple(pps)

        return AmpelAlert(
            id=int.from_bytes(  # alert id
                blake2b(all_ids, digest_size=7).digest(), byteorder=sys.byteorder
            ),
            stock=int(transient_name),  # internal ampel id
            datapoints=tuple(pps),
        )
