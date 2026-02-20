#!/usr/bin/env python
# File:                Ampel-LSST/ampel/lsst/t0/DecentVroFilter.py
# License:             BSD-3-Clause
# Author:              jno
# Date:                24.09.2025
# Last Modified Date:  24.09.2025
# Last Modified By:    jno

from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from ampel.abstract.AbsAlertFilter import AbsAlertFilter
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.base.CatalogMatchUnit import CatalogMatchUnit


class DecentVroFilter(CatalogMatchUnit, AbsAlertFilter):
    """
    General-purpose filter for VRO alerts, based on the ZTF DecentFilter.
    Goal is to provide a stream of likely extragalactic transients.

    As the VRO alerts mature and information is added, this should be
    included here. First iterations will be rough and mainly cut on
    alert history length.

    """

    # History
    min_ndet: int  # number of previous detections
    max_ndet: None | int  # number of previous detections
    min_ndet_cadence: float | None = (
        None  # time over which to coalesce previous detections when counting [days]
    )
    min_tspan: float  # minimum duration of alert detection history [days]
    max_tspan: float  # maximum duration of alert detection history [days]

    # Segments removed (for now)
    # min_archive_tspan: float = 0.0  # minimum duration of alert detection history [days]
    # max_archive_tspan: float = (
    #    10**5.0
    # )  # maximum duration of alert detection history [days]

    # Image quality, including Real-Bogus equivalents
    min_reliability: float = 0.0  # deep learning real bogus score
    # min_rb: float  # real bogus score
    # max_fwhm: float  # sexctrator FWHM (assume Gaussian) [pix]
    # max_elong: float  # Axis ratio of image: aimage / bimage
    # max_magdiff: float  # Difference: magap - magpsf [mag]
    # max_nbad: int  # number of bad pixels in a 5 x 5 pixel stamp

    # Astro
    # min_sso_dist: float  # distance to nearest solar system object [arcsec]
    min_abs_gal_lat: float  # minium absolute distance from galactic plane. Set to negative to disable cut.

    # PS1 - no PS info yet?
    # ps1_sgveto_rad: (
    #    float  # maximum distance to closest PS1 source for SG score veto [arcsec]
    # )
    # ps1_sgveto_th: (
    #    float  # maximum allowed SG score for PS1 source within PS1_SGVETO_RAD
    # )
    # ps1_confusion_rad: float  # reject alerts if the three PS1 sources are all within this radius [arcsec]
    # ps1_confusion_sg_tol: float  # and if the SG score of all of these 3 sources is within this tolerance to 0.5

    # Gaia
    gaia_rs: float  # search radius for GAIA DR2 matching [arcsec]
    gaia_pm_signif: float = (
        3.0  # significance of proper motion detection of GAIA counterpart [sigma]
    )

    gaia_plx_signif: float = (
        3.0  # significance of parallax detection of GAIA counterpart [sigma]
    )
    gaia_veto_gmag_min: float = (
        9  # min gmag for normalized distance cut of GAIA counterparts [mag]
    )
    gaia_veto_gmag_max: float = (
        20.0  # max gmag for normalized distance cut of GAIA counterparts [mag]
    )
    gaia_excessnoise_sig_max: float = 999.0  # maximum allowed noise (expressed as significance) for Gaia match to be trusted.

    def post_init(self):
        # feedback
        for k in self.__annotations__:
            self.logger.info(f"Using {k}={getattr(self, k)}")

        # To make this tenable we should create this list dynamically depending on what entries are required
        # by the filter. Now deciding not to include drb in this list, eg.
        self.keys_to_check = ("midpointMjdTai",)

    def _alert_has_keys(self, photop) -> bool:
        """
        check that given photopoint contains all the keys needed to filter
        """
        for el in self.keys_to_check:
            if el not in photop:
                self.logger.info(None, extra={"missing": el})
                return False
            if photop[el] is None:
                self.logger.info(None, extra={"isNone": el})
                return False
        return True

    def get_galactic_latitude(self, transient):
        """
        compute galactic latitude of the transient
        """
        coordinates = SkyCoord(transient["ra"], transient["dec"], unit="deg")
        return coordinates.galactic.b.deg

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
                any(gaia_tab["FLAG_PMRA"] == True)  # noqa: E712
                or any(gaia_tab["FLAG_PMDec"] == True)  # noqa: E712
                or any(gaia_tab["FLAG_Plx"] == True)  # noqa: E712
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

        pps = [el for el in alert.datapoints if el.get("diaSourceId") is not None]
        latest = pps[0]  # Most recent first - still the case?
        if not self._alert_has_keys(latest):
            return None

        #########
        ### Base cuts on lightcurve length / quality
        #########

        if self.min_reliability > 0 and (
            "reliability" not in latest
            or latest["reliability"] is None
            or latest["reliability"] < self.min_reliability
        ):
            self.logger.debug(
                None, extra={"reliability": latest.get("reliability", None)}
            )
            return None

        if self.min_ndet_cadence is not None:
            # thin detections to only count those separated by at least
            # min_ndet_cadence, to avoid counting multiple detections of the
            # same transient in a single night, for example from rapidfire
            # observations of the same field.
            thinned_pps: list[dict] = []
            for pp in sorted(pps, key=lambda x: x["midpointMjdTai"]):
                if (
                    not thinned_pps
                    or (pp["midpointMjdTai"] - thinned_pps[-1]["midpointMjdTai"])
                    >= self.min_ndet_cadence
                ):
                    thinned_pps.append(pp)
        else:
            thinned_pps = pps

        if len(thinned_pps) < self.min_ndet:
            self.logger.debug(None, extra={"nDet": len(thinned_pps)})
            return None
        if self.max_ndet and len(thinned_pps) > self.max_ndet:
            self.logger.debug(None, extra={"nDet": len(thinned_pps)})
            return None

        # cut on length of detection history
        detections_jds = [el["midpointMjdTai"] for el in pps]
        last_det = max(detections_jds)
        det_tspan = last_det - min(detections_jds)
        if not (self.min_tspan <= det_tspan <= self.max_tspan):
            self.logger.debug(None, extra={"tSpan": det_tspan})
            return None

        ## Solar System / Minor Planet

        # Currently we assume that diaObjects and ssSource are exclusie ,
        # such that solar system alerts are found here.
        # Have to revisit in case this needs to be inspected here.

        # ASTRONOMY
        ###########

        # cut on galactic latitude
        if (
            self.min_abs_gal_lat > 0
            and abs(b := self.get_galactic_latitude(latest)) < self.min_abs_gal_lat
        ):
            self.logger.debug(None, extra={"galPlane": abs(b)})
            return None

        # check ps1 star-galaxy score
        # if self.is_star_in_PS1(latest):
        #    # self.logger.debug("rejected: closest PS1 source %.2f arcsec away with sgscore of %.2f"% (latest['distpsnr1'], latest['sgscore1']))
        #    self.logger.debug(None, extra={"distpsnr1": latest["distpsnr1"]})
        #    return None

        # if self.is_confused_in_PS1(latest):
        #    # self.logger.debug("rejected: three confused PS1 sources within %.2f arcsec from alert."% (self.ps1_confusion_rad))
        #    self.logger.debug(None, extra={"ps1Confusion": True})
        #    return None

        # check with gaia
        if self.gaia_rs > 0 and self.is_star_in_gaia(latest):
            self.logger.debug(None, extra={"gaiaIsStar": True})
            return None

        self.logger.debug(
            "Alert accepted",
            extra={
                "midpointMjdTai": last_det,
                "latestdiaSourceId": latest["diaSourceId"],
            },
        )
        return True
