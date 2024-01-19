#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t2/T2Observability.py
# License:             BSD-3-Clause
# Author:              matteo.giomi@desy.de
# Date:                19.09.2018
# Last Modified Date:  05.02.2020
# Last Modified By:    valery brinnel <firstname.lastname@gmail.com>

from typing import Any

import astropy.units as u
from astropy.time import Time

from ampel.abstract.AbsLightCurveT2Unit import AbsLightCurveT2Unit
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.util.Observatory import Observatory
from ampel.view.LightCurve import LightCurve


class VisibilityConstraintModel(AmpelBaseModel):
    airmass_th: float = 2
    sun_alt_th: float = -12
    min_moon_dist: float = 30


class ObservatoryLocationModel(AmpelBaseModel):
    lat: float
    long: float
    alt: float


class ObservatoryModel(AmpelBaseModel):
    pos: ObservatoryLocationModel
    constraints: VisibilityConstraintModel = VisibilityConstraintModel()


class T2Observability(AbsLightCurveT2Unit):
    """
    cross match the position of a transient to those of sources in a set
    of catalogs and attach the required information to the transient.
    """

    # empty dict of named EarthLocation objetcs.
    observatories: dict[str, ObservatoryModel]

    # default parameters for LightCurve.get_pos method
    lc_get_pos_kwargs: dict[str, Any] = {"ret": "brightest", "filters": None}

    def post_init(self):
        self._observatories: dict[str, Observatory] = {}

    def init_observatory(self, name, pos: ObservatoryLocationModel) -> Observatory:
        """
        Return the earth location object corresponding to the desired observatory.
        Repeated requests to the same instance will not cause new duplicates.
        """

        # check if the observatory has already been init
        if not (obs := self._observatories.get(name)):
            self.logger.debug(
                f"Observatory {name} not previously instantiated. Doing it now."
            )

            # init your obs depending on the position
            self._observatories[name] = Observatory(
                name,
                logger=self.logger,
                latitude=pos.lat,
                longitude=pos.long,
                altitude=pos.alt,
            )
            return self._observatories[name]
        else:
            self.logger.debug(f"Observatory {name} already exists.")
            return obs

    def process(self, light_curve: LightCurve) -> UBson | UnitResult:
        """
        :param run_config: configuration parameter for this job.
        Eg:
        run_config = {
                'get_lc_pos_kwargs': None, # optional see ampel.view.LightCurve doc
                'observatories': {
                        'SEDm': {
                                pos: {
                                        'lat': 33.3483717		#from google maps, and for p48
                                        'long': -116.85972959
                                        'alt': 1680
                                },
                                constraints: {
                                        'airmass_th': 2,
                                        'sun_alt_th': -12,
                                        'min_moon_dist': 30
                                }
                        },
                        'SNIFS':
                                {.. },
                        ...
                        }
                }

        :returns: dict with the keys to append to each transient.
        {
                obs1: {
                                night1: {
                                        start,
                                        stop,
                                        },
                                night2: {
                                        ...
                                },
                                night3: {
                                }
                },
                obs2:
                        {
                        ...
                },
                ..
        }
        """

        self.logger.debug(
            f"getting transient position from lightcurve using args: {self.lc_get_pos_kwargs}"
        )

        if pos := light_curve.get_pos(**self.lc_get_pos_kwargs):
            transient_ra, transient_dec = pos
        else:
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        self.logger.debug(
            f"Transient position (ra, dec): {transient_ra:.4f}, {transient_dec:.4f} deg"
        )

        # initialize the catalog quer(ies). Use instance variable to aviod duplicates
        out_dict: dict[str, Any] = {}

        for name, observatory in self.observatories.items():
            my_obs = self.init_observatory(name, observatory.pos)

            for k in range(3):
                trange = [
                    (Time.now() + k * u.day).iso[:10],
                    (Time.now() + (k + 1) * u.day).iso[:10],
                ]

            ret = my_obs.compute_visibility(
                transient_ra,
                transient_dec,
                trange,
                airmass_th=observatory.constraints.airmass_th,
                sun_alt_th=observatory.constraints.sun_alt_th,
                min_moon_dist=observatory.constraints.min_moon_dist,
            )

            out_dict[name][f"night{k+1}"] = {"start": ret[0].iso, "end": ret[-1].iso}

        # return the info as dictionary
        return out_dict
