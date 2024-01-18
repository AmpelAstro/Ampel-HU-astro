#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2TNSEval.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                27.01.2021
# Last Modified Date:  17.03.2021
# Last Modified By:    jnordin@physik.hu-berlin.de

import numpy as np
from typing import Any
from collections.abc import Sequence
from astropy.coordinates import Distance, SkyCoord
from astropy.cosmology import Planck15

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit


class T2InfantCatalogEval(AbsTiedLightCurveT2Unit):
    """
    Evaluate whether a transient fulfills criteria for being a potentially
    infant (extragalactic) transient.

    This implementation requires a catalog match with redshift.
    """

    # Cuts based on T2 catalog redshifts

    # List of catalog-like output to search for redshift. It is assumed that
    # redshifts are stored as 'z'
    redshift_catalogs: list[str] = ['SDSS_spec', 'NEDz', 'GLADEv23', 'NEDz_extcats']  # Otherwise more
    # maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    max_redshift: float = 0.05 # 0.1
    # minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    min_redshift: float = 0.001
    # max abs mag through peak mag and redshift from catalog mach (require both)
    max_absmag: float = -12 # Originally -13, moved to -12 due to ZTF22aafoqrd
    # min abs mag through peak mag and redshift from catalog mach (require both)
    min_absmag: float = -20 # -17
    # arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
    min_dist: float = 1.5
    # arcsec, maximum distance
    max_dist: float = 50
    # kpc, maximum distance
    max_kpc_dist: float = 999

    # Cut on alert properties

    # A candidate need to have at least this many detections
    min_ndet: int = 1
    min_ndet_postul: int = 0  # and if it has this minimum nr of detection after the last significant (max_maglim) UL.

    # days, If a detection has an age older than this, skip (stars,age).
    max_age: float = 3.

    # Min age of detection history
    min_age: float = 0
    # range of peak magnitudes for submission
    min_peak_mag: float = 20
    max_peak_mag: float = 14
    # Reported detections in at least this many filters
    min_n_filters: int = 1
    # Require a detection in one of these filters (e.g. ZTF I-band more often spurious)
    det_filterids: list[int] = [1, 2, 3]   # default to any of them
    # Minimal galactic latitide
    min_gal_lat: float = 14
    # reject alert if ssdistnr smaller than this value for any pp
    ssdistnr_max: float = 1
    # reject alert if PS1 star for any pp
    ps1_sgveto_rad: float = 1
    ps1_sgveto_sgth: float = 0.8
    # Minimal median RB.
    rb_minmed: float = 0.3
    # Minimal median RB.
    drb_minmed: float = 0.995


    # Limiting magnitude to consider upper limits as 'significant'
    maglim_min: float = 19.5
    # A limiting magnitude max this time ago
    maglim_maxago: float = 2.5

    # Cut to apply to all the photopoints in the light curve.
    # This will affect most operations, i.e. evaluating the position,
    # computing number of detections ecc.
    lc_filters: list[dict[str, Any]] = [
        {"attribute": "sharpnr", "operator": ">=", "value": -10.15},
        {"attribute": "magfromlim", "operator": ">", "value": 0},
    ]

    def inspect_catalog(self, cat_res: dict[str, Any]) -> None | dict[str, Any]:
        """
        Check whether a redshift match can be found in matched catalogs.
        """

        # P1. # Loop through listed catalogs for redshift matches
        zmatchs = []
        info = {}
        for catname in self.redshift_catalogs:
            catinfo = cat_res.get(catname, False)
            if (
                catinfo and isinstance( catinfo.get('z',None), float) and 
                (self.min_redshift < catinfo["z"] < self.max_redshift) and
                (self.min_dist < catinfo["dist2transient"] < self.max_dist)
            ):
                self.logger.info(
                    "Found z.",
                    extra={
                        "catalog": catname,
                        "z": catinfo["z"],
                        "dist": catinfo["dist2transient"],
                    },
                )
                # Calculate physical distance
                dst_kpc = (
                    catinfo["dist2transient"] *
                    Planck15.kpc_proper_per_arcmin(catinfo["z"]).value / 60.0
                )
                if self.max_kpc_dist > 0 and dst_kpc > self.max_kpc_dist:
                    self.logger.info(
                        "Skip, physical distance too large.",
                        extra={"distance_kpc": dst_kpc},
                    )
                    continue
                zmatchs.append([catinfo["z"]])
                info[f"{catname}_z"] = catinfo["z"]
                info[f"{catname}_dist2transient"] = catinfo["dist2transient"]

        if len(zmatchs) == 0:
            return None
        info['zs'] = zmatchs

        # Special catalog searches - mark transients close to AGNs
        milliquas = cat_res.get("milliquas", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if milliquas and milliquas["redshift"] > 0:
            info["milliAGN"] = True
        if sdss_spec and sdss_spec["bptclass"] in [4, 5]:
            info["sdssAGN"] = True

        # Return collected info
        return info



    def inspect_lc(self, lc: LightCurve) -> None | dict[str, Any]:
        """
        Verify whether the transient lightcurve fulfill criteria for submission.

        """

        # apply cut on history: consider photophoints which are sharp enough
        pps = lc.get_photopoints(filters=self.lc_filters)
        assert pps is not None
        info: dict[str, Any] = {}

        # cut on number of detection
        if len(pps) < self.min_ndet:
            self.logger.info(
                'Rejected', extra={'det': len(pps)}
            )
            return None
        info["detections"] = len(pps)

        # cut on age
        jds = [pp["body"]["jd"] for pp in pps]
        most_recent_detection, first_detection = max(jds), min(jds)
        age = most_recent_detection - first_detection
        if age > self.max_age or age < self.min_age:
            self.logger.info('Rejected', extra={'age': age})
            return None
        info["age"] = age

        # cut on number of detection after last SIGNIFICANT UL
        ulims = lc.get_upperlimits(
            filters={
                "attribute": "diffmaglim",
                "operator": ">=",
                "value": self.maglim_min,
            }
        )

        if ulims and len(ulims) > 0:
            last_ulim_jd = sorted([x["body"]["jd"] for x in ulims])[-1]
            pps_after_ndet = lc.get_photopoints(
                filters=self.lc_filters + [{"attribute": "jd", "operator": ">=", "value": last_ulim_jd}]
            )
            # Check if there are enough positive detection after the last significant UL
            if (
                pps_after_ndet is not None and
                len(pps_after_ndet) < self.min_ndet_postul
            ):
                self.logger.info(
                    "not enough consecutive detections after last significant UL.",
                    extra={"NDet": len(pps), "lastUlimJD": last_ulim_jd},
                )
                return None
            # Check that there is a recent ul
            if (most_recent_detection - last_ulim_jd) > self.maglim_maxago:
                self.logger.info(
                    "No recent UL.",
                    extra={
                        "lastDet": most_recent_detection,
                        "lastUlimJD": last_ulim_jd,
                    },
                )
                return None
            info["last_UL"] = most_recent_detection - last_ulim_jd
        else:
            self.logger.info("no UL")
            return None

        # cut on number of filters
        used_filters = set([pp["body"]["fid"] for pp in pps])
        if len(used_filters) < self.min_n_filters:
            self.logger.info(
                "Rejected", extra={'nbr_filt': len(used_filters)}
            )
            return None
        # cut on which filters used
        if used_filters.isdisjoint(self.det_filterids):
            self.logger.info(
                "Rejected (wrong filter det)", extra={'det_filters': used_filters}
            )
            return None

        # cut on range of peak magnitude
        mags = [pp["body"]["magpsf"] for pp in pps]
        peak_mag = min(mags)
        if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
            self.logger.info(
                "Rejected", extra={'peak_mag': peak_mag}
            )
            return None
        info["peak_mag"] = peak_mag

        # For rapidly declining sources the latest magnitude is probably more relevant
        latest_pps = lc.get_photopoints(
            filters={
                "attribute": "jd",
                "operator": "==",
                "value": most_recent_detection,
            }
        )
        if latest_pps:
            if not len(latest_pps) == 1:
                raise ValueError("Have assumed a unique last photopoint")
            info["latest_mag"] = latest_pps[0]["body"]["magpsf"]

        # TODO: cut based on the mag rise per day (see submitRapid)

        # cut on galactic coordinates
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            raise ValueError("Light curve contains no points")
        coordinates = SkyCoord(ra, dec, unit="deg")
        b = coordinates.galactic.b.deg
        if abs(b) < self.min_gal_lat:
            self.logger.info(
                "Rejected (galactic plane)", extra={'gal_lat_b': b}
            )
            return None
        info["ra"] = ra
        info["dec"] = dec

        # cut on distance to closest solar system object
        # TODO: how to make this check: ('0.0' in list(phot["ssdistnr"])
        ssdist = np.array([pp["body"]["ssdistnr"] for pp in pps
            if "ssdistnr" in pp['body'].keys() and pp["body"]["ssdistnr"] is not None])
        close_to_sso = np.logical_and(ssdist < self.ssdistnr_max, ssdist > 0)

        # TODO: Note that this discards a transient if it was ever close to a ss object!
        if np.any(close_to_sso):
            self.logger.info(
                "Rejected (close to solar system object)",
                extra={"ssdistnr": ssdist.tolist()},
            )
            return None

        # check PS1 sg for the full alert history
        # Note that we for this check do *not* use the lightcurve filter criteria
        # TODO: Evaluate whether we should use the filters, and do a check for sufficient number of datapoints remaining
        if psdata := lc.get_tuples("distpsnr1", "sgscore1"):
            distpsnr1, sgscore1 = zip(*psdata)
            is_ps1_star = np.logical_and(
                np.array(distpsnr1) < self.ps1_sgveto_rad,
                np.array(sgscore1) > self.ps1_sgveto_sgth,
            )
            if np.any(is_ps1_star):
                self.logger.info(
                    "Rejected (PS1 SG cut)",
                    extra={"distpsnr1": distpsnr1, "sgscore1": sgscore1},
                )
                return None
        else:
            self.logger.info("No PS1 check as no data found.")

        # cut on median RB and DRB score
        rbs = [pp["body"]["rb"] for pp in pps]
        if np.median(rbs) < self.rb_minmed:
            self.logger.info(
                "Rejected (RB)",
                extra={"median_rd": np.median(rbs)},
            )
            return None
        elif (len(rbs) == 0) and self.rb_minmed > 0:
            self.logger.info("Rejected (No rb info)")
            return None
        info["rb"] = np.median(rbs)

        # drb might not exist
        drbs = [pp["body"]["drb"] for pp in pps if "drb" in pp["body"]]
        if len(drbs) > 0 and np.median(drbs) < self.drb_minmed:
            self.logger.info(
                "Rejected (dRB)",
                extra={"median_drd": np.median(drbs)},
            )
            return None
        elif (len(drbs) == 0) and self.drb_minmed > 0:
            self.logger.info("Rejected (No drb info)")
            return None

        info["drb"] = np.median(drbs)

        # Transient passed pure LC criteria
        self.logger.info("Passed T2infantCatalogEval", extra=info)
        return info


    # MANDATORY
    def process(self, light_curve: LightCurve, t2_views: Sequence[T2DocView]) -> UBson | UnitResult:
        """

        Evaluate whether a transient passes thresholds for being a nearby (young) transient.

        Parameters
        -----------
        light_curve: "ampel.view.LightCurve" instance.
        See the LightCurve docstring for more info.

        t2_views: List of T2Views (assumed to be the result of a CatalogMatch)

        Returns
        -------
        dict

        Containing transient info, and in particular the 'action' key. This will be set to true
        for transients passing all selection criteria.

        """


        # i. Check the catalog matching criteria
        # There might be multiple CatalogMatch associated with the transient
        # we here only take the first without specific origin
        t2_cat_match = t2_views[0]

        catalog_result = t2_cat_match.get_payload()
        if not isinstance(catalog_result, dict):
            return {'action': False, 'eval': 'No catlog match result'}
        transient_info = self.inspect_catalog(catalog_result)
        if not transient_info:
            return {'action': False, 'eval': 'No cat match in z-range'}


        # ii. Check whether the lightcurve passes selection criteria
        lc_info = self.inspect_lc(light_curve)

        if not lc_info:
            transient_info['action'] = False
            transient_info['eval'] = 'LC fail selection.'
            return transient_info
        transient_info.update(lc_info)


        # iii. Check absolute magnitude
        sndist = Distance(z=np.mean(transient_info['zs']), cosmology=Planck15)
        absmag = transient_info["peak_mag"] - sndist.distmod.value
        transient_info["absmag"] = absmag
        if not (self.min_absmag < absmag < self.max_absmag):
            self.logger.info("Rejected (absmag)", extra={"absmag": absmag})
            transient_info['action'] = False
            transient_info['eval'] = 'Absmag'
            return transient_info

        # Passed all criteria - ready for action
        transient_info['action'] = True
        transient_info['eval'] = 'Pass'


        return transient_info
