#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2LSPhotoZTap.py
# License           : BSD-3-Clause
# Author            : jnordin
# Date              : 20.04.2021
# Last Modified Date: 21.04.2021
# Last Modified By  : jnordin

from collections import OrderedDict
from collections.abc import Sequence
from functools import cached_property, partial
from io import BytesIO
from math import acos, cos, pi, sin
from typing import Any
from urllib.parse import urlparse, urlunparse

import backoff
import numpy as np
import requests
from astropy.io.votable import parse_single_table
from astropy.table import Table
from pandas import read_csv

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.enum.DocumentCode import DocumentCode
from ampel.secret.NamedSecret import NamedSecret
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


def convert(inp, outfmt="pandas", verbose=False, **kwargs):
    """
    *** Taken from datalab dl/helpers/util/convert ***
    (to avoid pulling in _all_ the dependencies)

    Convert input `inp` to a data structure defined by `outfmt`.

    Parameters
    ----------
    inp : str
        String representation of the result of a query. Usually this
        is a CSV-formatted string, but can also be, e.g. an
        XML-formatted votable (as string)

    outfmt : str
        The desired data structure for converting `inp` to. Default:
        'pandas', which returns a Pandas dataframe. Other available
        conversions are:

          string - no conversion
          array - Numpy array
          structarray - Numpy structured array (also called record array)
          table - Astropy Table
          votable - Astropy VOtable

        For outfmt='votable', the input string must be an
        XML-formatted string. For all other values, as CSV-formatted
        string.

    verbose : bool
        If True, print status message after conversion. Default: False

    kwargs : optional params
        Will be passed as **kwargs to the converter method.


    Example
    -------
    Convert a CSV-formatted string to a Pandas dataframe

    .. code-block:: python

       arr = convert(inp,'array')
       arr.shape  # arr is a Numpy array

       df = convert(inp,outfmt='pandas')
       df.head()  # df is as Pandas dataframe, with all its methods

       df = convert(inp,'pandas',na_values='Infinity') # na_values is a kwarg; adds 'Infinity' to list of values converter to np.inf

    """
    # When there are duplicate column names in the table, it would not work when converting to Astropy Table and Votable, so we
    # have to add '_n' as an identifier to the duplicate column names.
    index = inp.find("\n")
    header = inp[0:index]
    inp = inp[index + 1 :]
    list = header.split(",")
    col_dict: dict[str, int] = {}
    new_s = ""
    for l in list:
        if l in col_dict:
            n = col_dict[l]
            col_dict[l] = n + 1
            new_s += l + "_" + str(n) + ","
        else:
            new_s += l + ","
            col_dict[l] = 1
    inp = new_s[:-1] + "\n" + inp

    # map outfmt container types to a tuple:
    # (:func:`queryClient.query()` fmt-value, descriptive title,
    # processing function for the result string)
    mapping = OrderedDict(
        [
            (
                "string",
                ("csv", "CSV formatted table as a string", lambda x: x.getvalue()),
            ),
            (
                "array",
                (
                    "csv",
                    "Numpy array",
                    partial(np.loadtxt, unpack=False, skiprows=1, delimiter=","),
                ),
            ),
            (
                "structarray",
                (
                    "csv",
                    "Numpy structured / record array",
                    partial(np.genfromtxt, dtype=float, delimiter=",", names=True),
                ),
            ),
            ("pandas", ("csv", "Pandas dataframe", read_csv)),
            ("table", ("csv", "Astropy Table", partial(Table.read, format="csv"))),
            ("votable", ("votable", "Astropy VOtable", parse_single_table)),
        ]
    )

    if isinstance(inp, bytes):
        b = BytesIO(inp)
    elif isinstance(inp, str):
        b = BytesIO(inp.encode())
    else:
        raise TypeError("Input must be of bytes or str type.")

    output = mapping[outfmt][2](b, **kwargs)

    if isinstance(output, bytes):
        output = output.decode()

    if verbose:
        print(f"Returning {mapping[outfmt][1]}")  # noqa: T201

    return output


class T2LSPhotoZTap(AbsPointT2Unit):
    """
    Query the NOIR DataLab service for photometric redshifts from the
    Legacy Survey.

    Other queries can in principle be made as long as the string format parameters are the same (ra, dec, match_radius).

    """

    # Astro DataLab user id
    datalab_user: NamedSecret[str]
    datalab_pwd: NamedSecret[str]

    # Match parameters
    match_radius: float = 10  # in arcsec

    # Query. Candidate position and radius will be added
    query: str = "SELECT ra, dec, photo_z.z_phot_median, photo_z.z_phot_mean, photo_z.z_phot_std, photo_z.z_phot_l68, z_phot_u68, photo_z.z_spec, tractor.dered_mag_g, tractor.dered_mag_r, tractor.dered_mag_z, tractor.dered_mag_w1, tractor.dered_mag_w2 , tractor.dered_mag_w3, tractor.dered_mag_w4, tractor.snr_g, tractor.snr_r, tractor.snr_z, tractor.snr_w1, tractor.snr_w2, tractor.snr_w3, tractor.snr_w4 FROM ls_dr9.tractor as tractor JOIN ls_dr9.photo_z as photo_z on photo_z.ls_id = tractor.ls_id WHERE 't' = Q3C_RADIAL_QUERY(ra, dec,%.6f,%.6f,%.6f)"

    # run only on first datapoint by default
    # NB: this assumes that docs are created by DualPointT2Ingester
    ingest: dict = {"eligible": {"pps": "first"}}

    # Path to noir queries
    datalab_query_url: str = "https://datalab.noirlab.edu/query"

    @cached_property
    def session(self) -> requests.Session:
        """Obtain a session with an auth token"""
        session = requests.Session()
        parts = urlparse(self.datalab_query_url)
        response = session.get(
            urlunparse((parts.scheme, parts.netloc, "/auth/login", "", "", "")),
            params={  # type: ignore
                "username": self.datalab_user,
                "profile": "default",
                "debug": "False",
            },
            headers={"X-DL-Password": self.datalab_pwd},  # type: ignore
        )
        response.raise_for_status()
        session.headers.update({"X-DL-AuthToken": response.text})
        return session

    @backoff.on_exception(
        backoff.expo,
        requests.ConnectionError,
        max_tries=5,
        factor=10,
    )
    @backoff.on_exception(
        backoff.expo,
        requests.HTTPError,
        giveup=lambda e: isinstance(e, requests.HTTPError)
        and e.response.status_code not in {503, 429},
        max_time=60,
    )
    def _astrolab_query(self, ra: float, dec: float) -> Sequence[dict[str, Any]]:
        self.logger.debug(f"Querying {ra} {dec}")

        r = self.session.get(
            f"{self.datalab_query_url}/query",
            params={
                "sql": self.query % (ra, dec, self.match_radius / 3600),
                "ofmt": "csv",
                "async": "False",
            },
            timeout=300,
        )
        if not r.ok:
            self.logger.debug(f"DL query failed at {ra} {dec}" % (ra, dec))
            return []

        # First convert to string and then to dict
        ret_dict = convert(str(r.content.decode()), "pandas").to_dict(orient="records")

        self.logger.debug(f"Got {len(ret_dict)} matches")

        return ret_dict

    def add_separation(
        self, match_dict: Sequence[dict[str, Any]], target_ra: float, target_dec: float
    ) -> Sequence[dict[str, Any]]:
        """
        Iterate through catalog entries (dict) and add separation to target.
        """
        c = pi / 180

        for el in match_dict:
            if "dec" in el and "ra" in el:
                el["dist2transient"] = (
                    acos(
                        sin(target_dec * c) * sin(el["dec"] * c)
                        + cos(target_dec * c)
                        * cos(el["dec"] * c)
                        * cos((target_ra - el["ra"]) * c)
                    )
                    * 206264.8062
                )  # to arcsecs
            else:
                el["dist2transient"] = None

        return match_dict

    def process(self, datapoint: DataPoint) -> UBson | UnitResult:
        return_all: bool = (
            False  # whether to return all matches or closest match (default)
        )
        """
        Query DataLab through the unit query string combined with 
        transient position. 

        :returns: 
        {
        0:  {'ra': 247.033806830109,
             'dec': 63.8236957439316,
             'z_phot_median': 0.364699,
             'z_phot_mean': 0.640077,
             'z_phot_l68': 0.14177,
             'z_phot_u68': 1.30389,
             'z_spec': -99,
             ...,
             'dist2transient': 0.30846301868050724},
        1:  { ... },
        }


        Note that, when a match is found, the distance of the lightcurve object
        to the match counterpart is also returned as the 'dist2transient' key.
        """

        try:
            transient_ra = datapoint["body"]["ra"]
            transient_dec = datapoint["body"]["dec"]
        except KeyError:
            ##return T2RunState.MISSING_INFO
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        # Query Datalab
        match_list = self._astrolab_query(transient_ra, transient_dec)

        # Add separation between target and query detection
        if len(match_list) > 0:
            match_list = self.add_separation(match_list, transient_ra, transient_dec)

        # Return a T2 result (dict-like)
        if len(match_list) > 0:
            if return_all:
                return {f"T2LSPhotoZTap{k}": item for k, item in enumerate(match_list)}
            min_dist = min(match_list, key=lambda x: x["dist2transient"])
            return {"T2LSPhotoZTap": min_dist}
        return {"T2LSPhotoZTap": None}
