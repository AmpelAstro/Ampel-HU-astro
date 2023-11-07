#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/t3/RandomMapGenerator.py
# License:             BSD-3-Clause
# Author:              ernstand
# Date:                30.10.2023
# Last Modified Date:  30.10.2023
# Last Modified By:    <ernstand@physik.hu-berlin.de>

import numpy as np
import pandas as pd
import astropy as ap
import astropy.time as atime
from astropy import units as u
import healpy as hp
from typing import Tuple

from ampel.abstract.AbsT3PlainUnit import AbsT3PlainUnit
from ampel.struct.Resource import Resource
from ampel.struct.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson


class RandomMapGenerator(AbsT3PlainUnit):
    """
    Generate smoothed circular healpix probability values around a random coordinate.
    """

    save_dir: str = "."
    map_name: str = "simulated"

    long_range: Tuple[float, float] = [0, 360]
    lat_range: Tuple[float, float] = [0, 90]

    fwhm_range: Tuple[float, float] = [0.4, 1.2]

    nside: int = 32

    seed: int or None = None

    min_date: str = "2019-06-12"
    max_date: str = "2023-10-01"

    def process(self, t3s: T3Store) -> UBson | UnitResult:
        tmp_date = atime.Time(self.min_date)
        tmp_date.format = "gps"
        self.min_date_gps = tmp_date.value

        tmp_date = atime.Time(self.max_date)
        tmp_date.format = "gps"
        self.max_date_gps = tmp_date.value

        # generate random coordinates
        self.generate_randoms()

        # generate healpix map centered around coordinates
        self.generate_map()

        # save generated map
        self.save_map()

        resource = {
            "rand_latitude": self.latitude,
            "rand_longitude": self.longitude,
            "rand_fwhm": self.fwhm,
            "rand_seed": self.seed,
        }

        r = Resource(name="random_map", value=resource)
        t3s.add_resource(r)

        return None

    def generate_randoms(self):
        """Generate random value pair for longitude, latitude, fwhm, trigger time, distance to be used in map generation."""

        if self.seed:
            np.random.seed(self.seed)

        self.longitude = np.random.uniform(self.long_range[0], self.long_range[1])
        self.latitude = np.random.uniform(self.lat_range[0], self.lat_range[1])

        self.fwhm = np.random.uniform(self.fwhm_range[0], self.fwhm_range[1])

        self.distance = np.random.uniform(0, 2000)
        self.dist_unc = np.random.normal(240, 80)

        self.trigger_time = np.random.uniform(self.min_date_gps, self.max_date_gps)

    def generate_map(self):
        "Generate healpix probability values around coordinates"

        npix = hp.nside2npix(self.nside)

        pixels = np.zeros(npix)

        source_vector = hp.ang2vec(theta=self.longitude, phi=self.latitude, lonlat=True)
        # print(source_vector)

        radius = (3 * u.deg).to_value(u.radian)

        disc_pix = hp.query_disc(self.nside, source_vector, radius=radius, nest=True)
        # print(disc_pix)

        pixels[disc_pix] = 100

        smoothed_pix = hp.smoothing(pixels, fwhm=self.fwhm, nest=True)
        smoothed_pix /= np.sum(smoothed_pix)

        self.map_pix = smoothed_pix

    def save_map(self):
        """Save generated map with header as .fits.gz"""

        file_name = f"{self.save_dir}/{self.map_name}.fits.gz"

        print("Saving map: ", file_name)

        hdr = []
        hdr.append(("DISTMEAN", self.distance))
        hdr.append(("DISTSTD", self.dist_unc))

        trigger_time = atime.Time(self.trigger_time, format="gps")

        trigger_time.format = "fits"

        hdr.append(("DATE-OBS", trigger_time.value))
        hdr.append(("ORDERING", "NESTED"))

        hp.write_map(
            filename=file_name,
            nest=True,
            m=self.map_pix,
            extra_header=hdr,
            overwrite=True,
            dtype="float64",
        )
