#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU-astro/ampel/contrib/hu/t0/StellarFilter.py
# License:             BSD-3-Clause
# Author:              m. giomi <matteo.giomi@desy.de>
# Date:                06.06.2018
# Last Modified Date:  26.03.2024
# Last Modified By:    jno <jnordin@physik.hu-berlin.de.de>

from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.base.CatalogMatchUnit import CatalogMatchUnit


class StellarFilter(CatalogMatchUnit, AbsAlertFilter):
    """

    a.k.a. as the IndecentFilter, i.e. an inversion of the DecentFilter
    mainly used for finding extragalactic objects.


    todo:
    look for flare
    min mag deviation
    periodogram (from elasticc curve)

    """

    # History
    min_ndet: int  # number of previous detections
    max_ndet: int  # number of previous detections
    min_tspan: float  # minimum duration of alert detection history [days]
    max_tspan: float  # maximum duration of alert detection history [days]
    min_archive_tspan: float = 0.0  # minimum duration of alert detection history [days]
    max_archive_tspan: float = (
        10**5.0
    )  # maximum duration of alert detection history [days]

    # Brightness / Flare
    max_mag: float = 30.0
    peak_time_limit: float = (
        10.0  # Will divide lightcurve before / after this. Set to 0 to disably [days]
    )
    min_peak_diff: float = (
        1.0  # Min mag difference between peak mag before/after limit *in any band*
    )
    # Todo: select band? fit linear curve?

    # Image quality
    min_drb: float = 0.0  # deep learning real bogus score
    min_rb: float  # real bogus score
    max_fwhm: float = 5.5  # sexctrator FWHM (assume Gaussian) [pix]
    max_elong: float = 1.4  # Axis ratio of image: aimage / bimage
    max_magdiff: float = 1.0  # Difference: magap - magpsf [mag]
    max_nbad: int = 0  # number of bad pixels in a 5 x 5 pixel stamp

    # Astro
    min_sso_dist: float = 20  # distance to nearest solar system object [arcsec]
    min_gal_lat: float = (
        -1
    )  # minium distance from galactic plane. Set to negative to disable cut.
    max_gal_lat: float = 999  # maximum distance from galactic plane.

    # PS1
    require_ps_star: bool
    avoid_ps_confusion: bool = False  # Discard event if multiple nearby PS sources
    ps1_sgveto_rad: float = (
        1.0  # maximum distance to closest PS1 source for SG score veto [arcsec]
    )
    ps1_sgveto_th: float = (
        0.8  # maximum allowed SG score for PS1 source within PS1_SGVETO_RAD
    )
    ps1_confusion_rad: float = 1.0  # reject alerts if the three PS1 sources are all within this radius [arcsec]
    ps1_confusion_sg_tol: float = 0.1  # and if the SG score of all of these 3 sources is within this tolerance to 0.5

    # Gaia
    require_gaia_star: bool
    gaia_rs: float = 20.0  # search radius for GAIA DR2 matching [arcsec]
    gaia_pm_signif: float = (
        3.0  # significance of proper motion detection of GAIA counterpart [sigma]
    )
    gaia_plx_signif: float = (
        3.0  # significance of parallax detection of GAIA counterpart [sigma]
    )
    gaia_veto_gmag_min: float = (
        9.0  # min gmag for normalized distance cut of GAIA counterparts [mag]
    )
    gaia_veto_gmag_max: float = (
        20.0  # max gmag for normalized distance cut of GAIA counterparts [mag]
    )
    gaia_excessnoise_sig_max: float = 999.0  # maximum allowed noise (expressed as significance) for Gaia match to be trusted.

    def get_galactic_latitude(self, transient):
        """
        compute galactic latitude of the transient
        """
        coordinates = SkyCoord(transient["ra"], transient["dec"], unit="deg")
        return coordinates.galactic.b.deg

    def is_star_in_PS1(self, transient) -> bool:
        """
        apply combined cut on sgscore1 and distpsnr1 to reject the transient if
        there is a PS1 star-like object in it's immediate vicinity
        """

        # TODO: consider the case of alert moving wrt to the position of a star
        # maybe cut on the minimum of the distance!
        return (
            transient["distpsnr1"] < self.ps1_sgveto_rad
            and transient["sgscore1"] > self.ps1_sgveto_th
        )

    def is_confused_in_PS1(self, transient) -> bool:
        """
        check in PS1 for source confusion, which can induce subtraction artifatcs.
        These cases are selected requiring that all three PS1 cps are in the imediate
        vicinity of the transient and their sgscore to be close to 0.5 within given tolerance.
        """
        very_close = (
            max(transient["distpsnr1"], transient["distpsnr2"], transient["distpsnr3"])
            < self.ps1_confusion_rad
        )

        # Update 31.10.19: avoid costly numpy cast
        # Old:
        # In: %timeit abs(array([sg1, sg2, sg3]) - 0.5 ).max()
        # Out: 5.79 µs ± 80.5 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
        # New:
        # In: %timeit max(abs(sg1-0.5), abs(sg2-0.5), abs(sg3-0.5))
        # Out: 449 ns ± 7.01 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

        sg_confused = (
            max(
                abs(transient["sgscore1"] - 0.5),
                abs(transient["sgscore2"] - 0.5),
                abs(transient["sgscore3"] - 0.5),
            )
            < self.ps1_confusion_sg_tol
        )

        return sg_confused and very_close

    def is_star_in_gaia(self, transient: dict[str, Any]) -> bool:
        """
        match tranient position with GAIA DR2 and uses parallax
        and proper motion to evaluate star-likeliness
        returns: True (is a star) or False otehrwise.
        """

        srcs = self.cone_search_all(
            transient["ra"],
            transient["dec"],
            [
                {
                    "name": "GAIADR2",
                    "use": "catsHTM",
                    "rs_arcsec": self.gaia_rs,
                    "keys_to_append": [
                        "Mag_G",
                        "PMRA",
                        "ErrPMRA",
                        "PMDec",
                        "ErrPMDec",
                        "Plx",
                        "ErrPlx",
                        "ExcessNoiseSig",
                    ],
                }
            ],
        )[0]

        if srcs:
            gaia_tab = Table(
                [
                    {k: np.nan if v is None else v for k, v in src["body"].items()}
                    for src in srcs
                ]
            )

            # compute distance
            gaia_tab["DISTANCE"] = [src["dist_arcsec"] for src in srcs]
            gaia_tab["DISTANCE_NORM"] = (
                1.8 + 0.6 * np.exp((20 - gaia_tab["Mag_G"]) / 2.05)
                > gaia_tab["DISTANCE"]
            )
            gaia_tab["FLAG_PROX"] = [
                x["DISTANCE_NORM"]
                and self.gaia_veto_gmag_min <= x["Mag_G"] <= self.gaia_veto_gmag_max
                for x in gaia_tab
            ]

            # check for proper motion and parallax conditioned to distance
            gaia_tab["FLAG_PMRA"] = (
                abs(gaia_tab["PMRA"] / gaia_tab["ErrPMRA"]) > self.gaia_pm_signif
            )
            gaia_tab["FLAG_PMDec"] = (
                abs(gaia_tab["PMDec"] / gaia_tab["ErrPMDec"]) > self.gaia_pm_signif
            )
            gaia_tab["FLAG_Plx"] = (
                abs(gaia_tab["Plx"] / gaia_tab["ErrPlx"]) > self.gaia_plx_signif
            )

            # take into account precison of the astrometric solution via the ExcessNoise key
            gaia_tab["FLAG_Clean"] = (
                gaia_tab["ExcessNoiseSig"] < self.gaia_excessnoise_sig_max
            )

            # select just the sources that are close enough and that are not noisy
            gaia_tab = gaia_tab[gaia_tab["FLAG_PROX"]]
            gaia_tab = gaia_tab[gaia_tab["FLAG_Clean"]]

            # among the remaining sources there is anything with
            # significant proper motion or parallax measurement
            if (
                any(gaia_tab["FLAG_PMRA"] == True)  # noqa
                or any(gaia_tab["FLAG_PMDec"] == True)  # noqa
                or any(gaia_tab["FLAG_Plx"] == True)  # noqa
            ):
                return True

        return False

    # Override
    def process(self, alert: AmpelAlertProtocol) -> None | bool | int:
        """
        Mandatory implementation.
        To exclude the alert, return *None*
        To accept it, either return
        * self.on_match_t2_units
        * or a custom combination of T2 unit names
        """

        # CUT ON THE HISTORY OF THE ALERT
        #################################

        pps = [el for el in alert.datapoints if el.get("candid") is not None]
        if len(pps) < self.min_ndet or len(pps) > self.max_ndet:
            return None

        # cut on length of detection history
        detections_jds = [el["jd"] for el in pps]
        det_tspan = max(detections_jds) - min(detections_jds)
        if not (self.min_tspan <= det_tspan <= self.max_tspan):
            return None

        # IMAGE QUALITY CUTS
        ####################

        latest = alert.datapoints[0]

        if latest["isdiffpos"] == "f" or latest["isdiffpos"] == "0":
            return None

        if latest["rb"] < self.min_rb:
            return None

        if "drb" in latest and self.min_drb > 0.0 and latest["drb"] < self.min_drb:
            return None

        if latest["fwhm"] > self.max_fwhm:
            return None

        if latest["elong"] > self.max_elong:
            return None

        if abs(latest["magdiff"]) > self.max_magdiff:
            return None

        # cut on archive length
        if "jdendhist" in latest and "jdstarthist" in latest:
            archive_tspan = latest["jdendhist"] - latest["jdstarthist"]
            if not (self.min_archive_tspan < archive_tspan < self.max_archive_tspan):
                return None

        # Recent lightcurve brightness
        ###########

        if latest["magpsf"] > self.max_mag:
            return None

        pre_pp = [
            dp
            for dp in pps
            if "magpsf" in dp and (latest["jd"] - dp["jd"]) > self.peak_time_limit
        ]
        post_pp = [
            dp
            for dp in pps
            if "magpsf" in dp and (latest["jd"] - dp["jd"]) <= self.peak_time_limit
        ]
        if len(pre_pp) == 0:
            return None
        # Could also sort these for filter
        mdiff = np.mean([pp["magpsf"] for pp in pre_pp]) - np.mean(
            [pp["magpsf"] for pp in post_pp]
        )
        if mdiff < self.min_peak_diff:
            return None

        # ASTRONOMY
        ###########

        # check for closeby ss objects
        if 0 <= latest["ssdistnr"] < self.min_sso_dist:
            return None

        # cut on galactic latitude
        b = self.get_galactic_latitude(latest)
        if abs(b) < self.min_gal_lat:
            return None
        if abs(b) > self.max_gal_lat:
            return None

        # check ps1 star-galaxy score
        if self.require_ps_star and not self.is_star_in_PS1(latest):
            return None
        if self.avoid_ps_confusion and self.is_confused_in_PS1(latest):
            return None

        # check with gaia
        if self.require_gaia_star and self.is_star_in_gaia(latest):
            return None

        return True
