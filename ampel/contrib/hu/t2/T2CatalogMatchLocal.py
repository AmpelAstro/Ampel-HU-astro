#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel-hu-astro/ampel/contrib/hu/t2/T2CatalogMatchLocal.py
# License:             BSD-3-Clause
# Author:              matteo.giomi@desy.de
# Date:                24.08.2018
# Last Modified Date:  2.12.2021
# Last Modified By:    jn <jnordin@physik.hu-berlin.de>

from typing import Any, ClassVar, Literal
from functools import cached_property
from ampel.types import UBson

from extcats import CatalogQuery
from pymongo import MongoClient
from pymongo.errors import AutoReconnect
import backoff


from astropy.coordinates import SkyCoord
from astropy.table import Table
from extcats.catquery_utils import get_closest, get_distances
from numpy import asarray, degrees

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult

from ampel.model.DPSelection import DPSelection

from ampel.ztf.t2.T2CatalogMatch import CatalogModel



class ExtcatsUnit:
    """
    Interface to standard extcats
    """

#    require = ("extcats",)
    extcat_path: None | str = None

    max_query_time: float = 300

    @cached_property
    def catq_client(self):
        # As this is for local use we assume authorization is not needed
        #return MongoClient(self.resource["extcats"], **self.extcats.auth.get())
        return MongoClient(self.extcat_path)

    def get_extcats_query(self, catalog, logger, ra_key="ra", dec_key="dec", *args, **kwargs):
        q = CatalogQuery.CatalogQuery(
            catalog,
            *args,
            **kwargs,
            ra_key=ra_key,
            dec_key=dec_key,
            logger=logger,
            dbclient=self.catq_client,
        )
        decorate = backoff.on_exception(
            backoff.expo,
            AutoReconnect,
            logger=logger,
            max_time=self.max_query_time,
        )
        for method in 'binaryserach', 'findclosest', 'findwithin':
            setattr(
                q,
                method,
                decorate(getattr(q, method))
            )
        return q


class ExtcatsCatalogModel(CatalogModel):
    catq_kwargs: dict[str,Any] = {}
    use: Literal["extcats"]

class T2CatalogMatchLocal(ExtcatsUnit, AbsPointT2Unit):
    """
    Cross matches the position of a transient to those of sources in a set of catalogs

    :param catalogs: each value specifies a catalog in extcats or catsHTM format
     and the query parameters.

    .. note:: This requires access to a mongod hosting extcats databases and a
      filesystem containing catsHTM catalogs. For a similar unit that uses hosted
      services, see :class:`~ampel.ztf.t2.T2CatalogMatch`.
    """

    # run only on first datapoint by default
    eligible: ClassVar[DPSelection] = DPSelection(filter='PPSFilter', sort='jd', select='first')

    # Each value specifies a catalog in extcats or catsHTM format and the query parameters
    catalogs: dict[str, ExtcatsCatalogModel]

    # Default behavior is to return the closest match, but can be switched to returning all
    closest_match: bool = True

    # Extcat config from context. Could add auhtorization
    require = ("extcats",)



    def post_init(self):

        # Select extcat
        assert self.resource
        self.extcat_path = self.resource["extcats"]

        # dict for catalog query objects
        self.catq_objects = {}

    def init_extcats_query(self, catalog, **kwargs) -> "CatalogQuery":
        """
        Return the extcats.CatalogQuery object corresponding to the desired
        catalog. Repeated requests to the same catalog will not cause new duplicated
        CatalogQuery instances to be created.
        """

        # check if you have already init this peculiar query
        if not (catq := self.catq_objects.get(catalog)):
            catq = self.get_extcats_query(catalog, self.logger, **kwargs)
            self.catq_objects[catalog] = catq
        return catq

    def process(self, datapoint: DataPoint) -> UBson | UnitResult:
        """
        :returns: example of a match in SDSS but not in NED:

        {
            'SDSS_spec': {
                'z': 0.08820018172264099,
                'bptclass': 2.0,
                'subclass': '',
                'dist2transient': 1.841666956181802e-09}
            },
            'NED': False
        }

        Note that, when a match is found, the distance of the lightcurve object
        to the catalog counterpart is also returned as the 'dist2transient' key.
        """

        try:
            transient_ra = datapoint["body"]["ra"]
            transient_dec = datapoint["body"]["dec"]
        except KeyError:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        self.logger.debug(
            "Transient position (ra, dec): {transient_ra:.4f}, {transient_dec:.4f} deg"
        )

        # initialize the catalog quer(ies). Use instance variable to aviod duplicates
        out_dict: dict[str, Any] = {}
        for catalog, cat_opts in self.catalogs.items():

            src = None
            self.logger.debug(
                f"Loading catalog {catalog} using options: {str(cat_opts)}"
            )
            # find out how ra/dec are called in the catalog
            catq_kwargs = cat_opts.catq_kwargs
            if catq_kwargs is None:
                ra_key, dec_key = "ra", "dec"
            else:
                ra_key, dec_key = (
                    catq_kwargs.get("ra_key", "ra"),
                    catq_kwargs.get("dec_key", "dec"),
                )


            # how do you want to support the catalog?
            if cat_opts.use == "extcats":

                # get the catalog query object and do the query
                catq = self.init_extcats_query(catalog, **cat_opts.catq_kwargs)

                if self.closest_match:
                    src, dist = catq.findclosest(
                        transient_ra,
                        transient_dec,
                        cat_opts.rs_arcsec,
                        pre_filter=cat_opts.pre_filter,
                        post_filter=cat_opts.post_filter,
                    )
                    # Sinlge row, cannot just add column
                    if src is not None:
                        src = [dict(src)]
                        src[0]["dist2transient"] = dist
                else:
                    src = catq.findwithin(
                        transient_ra,
                        transient_dec,
                        cat_opts.rs_arcsec,
                        pre_filter=cat_opts.pre_filter,
                        post_filter=cat_opts.post_filter,
                    )
                    if src is not None:
                        # get distances to all sources
                        dist = get_distances(
                            transient_ra,
                            transient_dec,
                            src,
                            ra_key,
                            dec_key,
                        )
                        src["dist2transient"] = dist

            elif cat_opts.use == "catsHTM":

                # catshtm needs coordinates in radians
                transient_coords = SkyCoord(transient_ra, transient_dec, unit="deg")
                srcs, colnames, colunits = self.catshtm.cone_search(
                    catalog,
                    transient_coords.ra.rad,
                    transient_coords.dec.rad,
                    cat_opts.rs_arcsec,
                )

                if len(srcs) > 0:   # JN: Always true if at least one match is done?

                    # format to astropy Table
                    srcs_tab = Table(asarray(srcs), names=colnames)

                    # convert distances to degrees (catsHTM stuff is in radians)
                    srcs_tab[ra_key] = degrees(srcs_tab[ra_key])
                    srcs_tab[dec_key] = degrees(srcs_tab[dec_key])

                    if self.closest_match:
                        # get the closest source and its distance (catsHTM stuff is in radians)
                        src, dist = get_closest(
                            transient_coords.ra.degree,
                            transient_coords.dec.degree,
                            srcs_tab,
                            ra_key,
                            dec_key,
                        )
                        # Sinlge row, cannot just add column
                        if src is not None:
                            src = [dict(src)]
                            src[0]["dist2transient"] = dist

                    else:
                        # get distances to all sources
                        src = srcs_tab
                        dist = get_distances(
                            transient_ra,
                            transient_dec,
                            src,
                            ra_key,
                            dec_key,
                        )
                        src["dist2transient"] = dist
            else:
                raise ValueError(
                    f"use option can not be {cat_opts.use} for catalog {catalog}. "
                    f"valid are 'extcats' or 'catsHTM'"
                )

            if src is not None:
                # if you found a cp add the required field from the catalog:
                # if keys_to_append argument is given or if it is equal to 'all'
                # then take all the columns in the catalog. Otherwise only add the
                # requested ones.

                keys_to_append = set(
                    cat_opts.keys_to_append if cat_opts.keys_to_append else src.colnames
                )
                if cat_opts.use == "extcats":
                    keys_to_append.difference_update({"_id", "pos"})
                # Always include dist2transient
                keys_to_append.add("dist2transient")



                out_dict[catalog] = []


                # Can skip full loop if we know there are no errors or weird formats?
                for one_src in src:

                    to_add = {}
                    for field in keys_to_append:
                        try:
                            val = one_src[field].tolist()
                        except AttributeError:
                            val = one_src[field]
                        except KeyError:
                            continue
                        to_add[field] = val
                    out_dict[catalog].append(to_add)

                # To be backwards compatible we directly return dict if only one match was requested, a list otherwise
                if self.closest_match:
                    # Above checks should have returned a dict with exactly one entry
                    assert len(out_dict[catalog]) == 1
                    out_dict[catalog] = out_dict[catalog][0]

            else:
                self.logger.debug(
                    f"no match found in catalog {catalog} within "
                    f"{cat_opts.rs_arcsec:.2f} arcsec from transient"
                )
                out_dict[catalog] = False

        # return the info as dictionary
        return out_dict
