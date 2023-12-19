#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/util/AmpelHealpix.py
# License:             BSD-3-Clause
# Author:              jakob nordin
# Date:                27.03.2023
# Last Modified Date:  27.03.2023
# Last Modified By:    jnordin


import math
import os
from base64 import b64decode, b64encode
from collections import defaultdict
from datetime import datetime
from hashlib import blake2b

import healpy as hp
import numpy as np
import requests
from astropy.time import Time
from astropy.io.fits import Header


class AmpelHealpix:
    """
    - Obtain and load (GW) Healpix map.
    - Get pixels above some probability threshold.
    - Get probability for some (ra,dec).
    - Get mean distance & distance uncertainty
    """

    # Disk storage
    save_dir: str = "."

    def __init__(
        self, map_name: str, map_url: None | str = None, save_dir: None | str = None, nside: None | int = None
    ):
        self.map_name = map_name
        self.map_url = map_url
        print(self.map_url)
        if save_dir:
            self.save_dir = save_dir
        self.nside = nside

        self._get_map()
        # Attribues
        self.credible_levels: None | list = None
        self.trigger_time: None | float = None

    def _get_map(self, clobber=False) -> int:
        path = os.path.join(self.save_dir, self.map_name)
        print(path)
        if os.path.exists(path) and not clobber:
            print("Map exists and found")
            return 1

        # Retrieve mapfile.
        map_data = requests.get(self.map_url)
        with open(path, "wb") as fh:
            fh.write(map_data.content)

        return 0

    def process_map(self) -> str:
        """
        Load map and determine prob values.
        """

        print(os.path.join(self.save_dir, self.map_name))


        # Process map
        hpx, headers = hp.read_map(
            os.path.join(self.save_dir, self.map_name), h=True, nest=True
        )
        #print("HEADERS")
        #print(headers)
        #print(type(headers))
        #print("END HEADERS")

        trigger_time = [
            datetime.fromisoformat(header[1])
            for header in headers
            if header[0] == "DATE-OBS"
        ][0]

        nside = int(hp.npix2nside(len(hpx)))

        # Downgrade resolution 
        if self.nside and self.nside<nside:
            hpx = hp.ud_grade( hpx, nside_out=self.nside,order_in='NESTED',order_out='NESTED',power=-2)
        else:
            self.nside = nside
           
        #print(headers)
        self.dist = [header[1] for header in headers if header[0] == "DISTMEAN"][0]
        #for header in headers:
            #print(header[1], header[0])
        
        #print(self.dist)
        self.dist_unc = [header[1] for header in headers if header[0] == "DISTSTD"][0]
        
        seed_list = [header[1] for header in headers if header[0] == "SEED"]
        if len(seed_list) > 0: 
            self.seed = seed_list[0]
        else:
            self.seed = self.map_name.replace(".fits.gz", "").replace(",", "_")

        # Find credible levels
        idx = np.flipud(np.argsort(hpx))

        sorted_credible_levels = np.cumsum(hpx[idx].astype(float))
        credible_levels = np.empty_like(sorted_credible_levels)
        credible_levels[idx] = sorted_credible_levels

        self.credible_levels = credible_levels
        self.trigger_time = Time(trigger_time).jd

        return b64encode(
            blake2b(sorted_credible_levels, digest_size=7).digest()
        ).decode("utf-8")

    def get_triggertime(self):

        if self.trigger_time is None:
            self.process_map()
        return self.trigger_time

    def get_pixelmask(self, pvalue_limit):
        """
        Return pixels with total probability up to some limit.
        """

        if not self.nside:
            raise ValueError("First get and process map before using.")

        # Create mask for pixel selection
        mask = np.zeros(len(self.credible_levels), int)
        mask[self.credible_levels <= pvalue_limit] = 1

        return mask.nonzero()[0].tolist()
    
    def get_maparea(self, pvalue_limit: float) -> float:
        """
        Return sky area for probability threshold in square degrees.
        """
        
        pixels = self.get_pixelmask(pvalue_limit)
        # Combine pixels when possible
        deresdict = deres(self.nside, pixels)
        healpix_regions = [
            {"nside": nside, "pixels": members} for nside, members in deresdict.items()
        ]
        # calculate relevant map area
        hp_area = 0
        for region in healpix_regions:
            npix_from_nside = 12 * region["nside"]**2
            hp_area += len(region["pixels"]) / npix_from_nside
        hp_area *= 360**2 / np.pi
        #print("USED MAP AREA: ", hp_area)

        return hp_area

    def get_cumprob(self, ra: float, dec: float) -> float:
        """
        Obtain probability for a specific coordinate based on loaded map.
        ra, dec in degrees.
        """

        if not self.nside:
            raise ValueError("First get and process map before using.")

        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        alertpix = hp.pixelfunc.ang2pix(
            hp.npix2nside(len(self.credible_levels)), theta, phi, nest=True
        )
        return self.credible_levels[alertpix]
    
    def get_mapdist(self) -> tuple[float, float]:
        """
        Obtain mean distance + std from healpix map in Mpc.
        """
        if not self.nside:
            raise ValueError("First get and process map before using.")
        
        return self.dist, self.dist_unc



def deres(nside, ipix, min_nside=1):
    """
    Decompose a set of (nested) HEALpix indices into sets of complete superpixels at lower resolutions.
    :param nside: nside of given indices
    :param ipix: pixel indices
    :min_nside: minimum nside of complete pixels

    Copied from Ampel-ZTF-archive.
    """
    remaining_pixels = set(ipix)
    decomposed = defaultdict(list)
    for log2_nside in range(int(math.log2(min_nside)), int(math.log2(nside)) + 1):
        super_nside = 2**log2_nside
        # number of base_nside pixels per nside superpixel
        scale = (nside // super_nside) ** 2
        # sort remaining base_nside pixels by superpixel
        by_superpixel = defaultdict(list)
        for pix in remaining_pixels:
            by_superpixel[pix // scale].append(pix)
        # represent sets of pixels that fill a superpixel
        # as a single superpixel, and remove from the working set
        for superpix, members in by_superpixel.items():
            if len(members) == scale:
                decomposed[super_nside].append(superpix)
                remaining_pixels.difference_update(members)

    return dict(decomposed)


def main():
    #map_name = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/new_BBH_mergers_O3a_IAS_pipeline/skymapGW190704_104834.fits.gz"
    #map_name = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/ligo-kilonova/tmp/S190408an.fits.gz,0"
    map_name = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/Ampel-HU-astro/simulated_maps/tmp_gaussmap.fits.gz"
    map_name = "/mnt/c/Users/Public/Documents/Uni/master/masterarbeit/ampel/ampel-results/ligo-kilonova/tmp/S190426c.fits.gz,1"
    ah = AmpelHealpix(
        #map_name="S191222n.fits.gz",
        map_name=map_name,
        #map_url="https://gracedb.ligo.org/api/superevents/S191222n/files/LALInference.fits.gz",
    )
    hashit = ah.process_map()
    print(hashit)
    pixels = ah.get_pixelmask(0.9)
    print(ah.get_maparea(.9))
    print(ah.get_mapdist())

    print(ah.trigger_time)

    print(ah.get_cumprob(13, 48))
    print(ah.get_cumprob(200, -14.88))
    print(ah.get_cumprob(68.58494, 33.89))


if __name__ == "__main__":
    main()
