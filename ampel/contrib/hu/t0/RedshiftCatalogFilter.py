#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-HU_astro/ampel/contrib/hu/t0/RedshiftCatalogFilter.py
# License:             BSD-3-Clause
# Author:              Jakob Nordin <jnordin@physik.hu-berlin.de>
# Date:                02.11.2022
# Last Modified Date:  02.11.2022
# Last Modified By:    jnordin <jnordin@physik.hu-berlin.de>

from functools import partial
from typing import Literal

from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol
from ampel.ztf.t0.DecentFilter import DecentFilter


class RedshiftCatalogFilter(DecentFilter):
    """
    Filter derived from DecentFilter designed to only accept transients
    located close to a galaxy in a catalog, and within redshift bounds.

    Default parameters focused on matching with nearby catalogs.
    Current generation of filter uses a copy of the NED galaxy archive
    to z 0.02x (0.03?), downloaded ~2020.

    If any match within the radius is found, the candidate will be stat_accepted
    (even if it, for example, is not the closest match).

    Many DecentFilter parameters are set as defaults, including turning
    off matching to Gaia for stars (it is assumed that the catalog matching
    is sufficient for rejection).

    Notes: Much catalog match functions, including combining multiple catalogs
    already exist in CatalogMatchFilter. Here we wished to also make use of
    DecentFilter setup, hence the partial code repetition.
    Similarly, extcat filters could use the post_filter option for immediate
    filtering, but it is not sure how general this can be made.
    """

    # Catalog matching parameters
    catalog_name: str = "NEDz_extcats"
    catalog_type: Literal["extcats", "catsHTM"] = "extcats"
    catalog_match_radius: float = 60.0  # Match radius in arcsec
    catalog_zkey: str = "z"  # Should match redshift key in catalog
    min_z: float = 0.002  # Some stars found in e.g. NED.
    max_z: float = 0.03  # Max redshift (after catalog match)
    # Remaining decent filter parameters
    max_tspan: float  # maximum duration of alert detection history [days]
    max_archive_tspan: float  # maximum duration of alert detection history [days]
    min_drb: float = 0.0  # deep learning real bogus score
    min_rb: float  # real bogus score

    # DecentFilter default parameter overrides to match this purpose
    # Tuned based on 2021 infant SN search comparison (contact J Nordin)
    min_ndet: int = 0  # number of previous detections
    min_tspan: float = -666.0  # minimum duration of alert detection history [days]
    min_archive_tspan: float = (
        -666.0
    )  # minimum duration of alert detection history [days]
    max_fwhm: float = 4.5  # sexctrator FWHM (assume Gaussian) [pix]
    max_elong: float = 1.4  # Axis ratio of image: aimage / bimage
    max_magdiff: float = 1  # Difference: magap - magpsf [mag]
    max_nbad: int = 0  # number of bad pixels in a 5 x 5 pixel stamp
    min_sso_dist: float = 20.0  # distance to nearest solar system object [arcsec]
    min_gal_lat: float = (
        14.0  # minium distance from galactic plane. Set to negative to disable cut.
    )
    ps1_sgveto_rad: float = (
        2.0  # maximum distance to closest PS1 source for SG score veto [arcsec]
    )
    ps1_sgveto_th: float = (
        0.8  # maximum allowed SG score for PS1 source within PS1_SGVETO_RAD
    )
    ps1_confusion_rad: float = 3.0  # reject alerts if the three PS1 sources are all within this radius [arcsec]
    ps1_confusion_sg_tol: float = 0.1  # and if the SG score of all of these 3 sources is within this tolerance to 0.5
    # Gaia matching is here turned off, standard match radius is 20.
    gaia_rs: float = -999.0  # search radius for GAIA DR2 matching [arcsec]
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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        query = {
            "name": self.catalog_name,
            "use": self.catalog_type,
            "rs_arcsec": self.catalog_match_radius,
            "keys_to_append": [self.catalog_zkey],
        }
        # For exctcats, directly add the z filter
        if self.catalog_type == "extcats":
            query["post_filter"] = {
                self.catalog_zkey: {"$gte": self.min_z, "$lte": self.max_z}
            }

        self.catalog_query = partial(self.cone_search_all, catalogs=[query])

    def process(self, alert: AmpelAlertProtocol):
        """
        Run the filter on the alert. First we run the decent filter, then we match
        with the infant catalog.
        """

        # First make default DefentFilter check
        if not (decent_result := super().process(alert)):
            return decent_result

        # if the candidate has passed the decent filter, check if it is compatible
        # with the position of some nearby galaxy cluster
        latest = alert.datapoints[0]
        alert_ra = latest["ra"]
        alert_dec = latest["dec"]

        # Check for catalog matches .
        catalog_matches = self.catalog_query(alert_ra, alert_dec)[0]
        if catalog_matches is None:
            return None

        # For all catalog matches, check redshift range
        matched = None
        for catmatch in catalog_matches:
            # CatalogItem object
            if matched := (
                self.min_z <= catmatch["body"][self.catalog_zkey] <= self.max_z
            ):
                break

        return matched
