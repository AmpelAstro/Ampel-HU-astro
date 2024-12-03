#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/t3/RandomMapGenerator.py
# License:             BSD-3-Clause
# Author:              ernstand
# Date:                30.10.2023
# Last Modified Date:  30.10.2023
# Last Modified By:    <ernstand@physik.hu-berlin.de>


import astropy.time as atime
import healpy as hp
import numpy as np
from astropy import units as u

from ampel.abstract.AbsT4Unit import AbsT4Unit
from ampel.types import UBson


class RandomMapGenerator(AbsT4Unit):
    """
    Generate smoothed circular healpix probability values around a random coordinate.
    """

    save_dir: str = "."
    map_name: str = "simulated"

    long_range: tuple[float, float] = (0, 360)
    lat_range: tuple[float, float] = (-30, 90)

    fwhm_range: tuple[float, float] = (0.4, 1.2)

    nside: int = 32

    seed: int | None = None

    min_date: str = "2019-06-12"
    max_date: str = "2023-10-01"

    def do(self) -> UBson:
        """Generate random coordinates, generate map around coordinates, save map and generated randoms"""

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

        return {"random_map": resource}

    def generate_randoms(self):
        """Generate random value pair for longitude, latitude, fwhm, trigger time, distance to be used in map generation."""

        if self.seed:
            np.random.seed(self.seed)
        else:
            self.seed = np.random.random_integers(0, 2147483647)
            self.logger.debug("RANDOMMAPGENERATOR", extra={"seed": self.seed})
            np.random.seed(self.seed)

        tmp_date = atime.Time(self.min_date)
        tmp_date.format = "gps"
        self.min_date_gps = tmp_date.value

        tmp_date = atime.Time(self.max_date)
        tmp_date.format = "gps"
        self.max_date_gps = tmp_date.value

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

        file_name = f"{self.save_dir}{(self.map_name).replace('.fits.gz', '')}.fits.gz"

        self.logger.debug(f"Saving map: {file_name}")

        hdr_simple = []
        hdr_simple.append(("SIMPLE", "T"))

        hdr: list[tuple[str, float | int | str | None]] = []

        hdr.append(("DISTMEAN", self.distance))
        hdr.append(("DISTSTD", self.dist_unc))
        hdr.append(("SEED", self.seed))

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
