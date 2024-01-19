#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t2/T2LSPhotoZTap.py
# License           : BSD-3-Clause
# Author            : jnordin
# Date              : 20.04.2021
# Last Modified Date: 21.04.2021
# Last Modified By  : jnordin

from typing import Any, Sequence, TYPE_CHECKING, Union

from astropy.table import Table

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint

from ampel.struct.UnitResult import UnitResult
from ampel.enum.DocumentCode import DocumentCode
from ampel.types import UBson

from math import cos, sin, acos, pi

# Datalab
from dl import authClient

from collections import OrderedDict
import numpy as np
from pandas import read_csv
from astropy.table import Table
from astropy.io.votable import parse_single_table
from functools import partial
from io import BytesIO

import backoff
import requests


def convert(inp, outfmt="pandas", verbose=False, **kwargs):
    """
    *** Taken from datalab dl/helpers/util/convert ***
    (to avoid loading alld ependencies)

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
        if l in col_dict.keys():
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
        print("Returning %s" % mapping[outfmt][1])

    return output


class T2LSPhotoZTap(AbsPointT2Unit):
    """
    Query the NOIR DataLab service for photometric redshifts from the
    Legacy Survey.

    Other queries can in principle be made as long as the string format parameters are the same (ra, dec, match_radius).

    """

    # Astro DataLab user id
    datalab_user: str
    datalab_pwd: str
    ##datalab_str : Secret

    # Match parameters
    match_radius: float = 10  # in arcsec

    # Query. Candidate position and radius will be added
    query: str = "SELECT ra, dec, photo_z.z_phot_median, photo_z.z_phot_mean, photo_z.z_phot_l68, z_phot_u68, photo_z.z_spec, tractor.dered_mag_g, tractor.dered_mag_r, tractor.dered_mag_z, tractor.dered_mag_w1, tractor.dered_mag_w2 , tractor.dered_mag_w3, tractor.dered_mag_w4, tractor.snr_g, tractor.snr_r, tractor.snr_z, tractor.snr_w1, tractor.snr_w2, tractor.snr_w3, tractor.snr_w4 FROM ls_dr8.tractor as tractor JOIN ls_dr8.photo_z as photo_z on photo_z.ls_id = tractor.ls_id WHERE 't' = Q3C_RADIAL_QUERY(ra, dec,%.6f,%.6f,%.6f)"

    # run only on first datapoint by default
    # NB: this assumes that docs are created by DualPointT2Ingester
    ingest: dict = {"eligible": {"pps": "first"}}

    # Path to noir queries
    datalab_query_url: str = "https://datalab.noirlab.edu/query"

    def post_init(self) -> None:
        # obtain security token
        self.token = authClient.login(self.datalab_user, self.datalab_pwd)

    #        self.token = authClient.login(self.datalab_user, self.datalab_str)

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
    def _astrolab_query(
        self, ra: float, dec: float
    ) -> Sequence[
        dict[str, Any]
    ]:  # Does one need to add List[None] here for empty returns?
        self.logger.debug("Querying %s %s" % (ra, dec))

        # Original qery, using the wrong path (and also a lot of dependencies)
        ## should be possible to adjust this:
        ## qc = queryClient.queryClient()
        ## qc.set_svc_url( 'https://datalab.noirlab.edu/query' )
        ## But this yiels a split error
        # Old query
        # ret = queryClient.query(self.token, self.query % (ra, dec, float(self.match_radius) / 3600) )
        # ret_dict = convert(ret,'pandas').to_dict(orient='records')
        # New manual:
        headers = {"X-DL-AuthToken": (self.token)}
        sql = self.query % (ra, dec, float(self.match_radius) / 3600)
        qfmt = "csv"
        async_ = False
        timeout = 300
        dburl = "%s/query?sql=%s&ofmt=%s&async=%s" % (
            self.datalab_query_url,
            sql,
            qfmt,
            async_,
        )
        r = requests.get(dburl, headers=headers, timeout=timeout)
        if not r.ok:
            self.logger.debug("DL query failed at %s %s" % (ra, dec))
            return []

        # First convert to string and then to dict
        ret_dict = convert(str(r.content.decode()), "pandas").to_dict(orient="records")

        self.logger.debug("Got %s matches" % (len(ret_dict)))

        return ret_dict

    def add_separation(
        self, match_dict: Sequence[dict[str, Any]], target_ra: float, target_dec: float
    ) -> Sequence[dict[str, Any]]:
        """
        Iterate through catalog entries (dict) and add separation to target.
        """
        c = pi / 180

        for el in match_dict:
            if "dec" in el.keys() and "ra" in el.keys():
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

    def process(self, datapoint: DataPoint) -> Union[UBson, UnitResult]:
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
            if return_all == True:
                return {
                    "T2LSPhotoZTap{}".format(k): item
                    for k, item in enumerate(match_list)
                }
            else:
                min_dist = min(match_list, key=lambda x: x["dist2transient"])
                return {"T2LSPhotoZTap": min_dist}
        else:
            return {"T2LSPhotoZTap": None}
