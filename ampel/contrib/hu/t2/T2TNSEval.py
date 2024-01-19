#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2TNSEval.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                27.01.2021
# Last Modified Date:  27.01.2021
# Last Modified By:    jnordin@physik.hu-berlin.de

from collections.abc import Sequence
from typing import Any

import numpy as np
from astropy.coordinates import SkyCoord

# T2 importing info from T3. Restructure?
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.contrib.hu.t3.ampel_tns import (
    TNSFILTERID,
)
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView


class T2TNSEval(AbsTiedLightCurveT2Unit):
    """
    Evalute whether a transient fulfills criteria for submission to TNS.
    This is done partially based on the lightcurve content and partially on the results of catalog maches.
    The T2 result will be a dictionary with the information used for TNS submission (e.g. recent upper limit).
    """

    # Lightcurve inspection parameters
    # Limiting magnitude to consider upper limits as 'significant'
    max_maglim: float = 19.5
    # Number of photometric detection we include in the TNS AT report
    nphot_submit: int = 2

    # cuts on T2 catalogs
    # reject candidates if they don't have matching in this list of T2CATALOGMATCH catalogs
    needed_catalogs: list[str] = []
    require_catalogmatch: bool = False
    # maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    max_redshift: float = 1.15
    # minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    min_redshift: float = 0
    # arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
    start_dist: float = 1.5
    # reject transient if the GAIA source brighter than this is nearby.
    max_gaia_neighbour_gmag: float = 11

    # cut on alert properties
    # A candidate need to have at least this many detections
    min_ndet: int = 2
    # and if it has this minimum nr of detection after the last significant (max_maglim) UL.
    min_ndet_postul: int = 2
    # days, If a detection has an age older than this, skip (stars,age).
    max_age: float = 50
    # Min age of detection history
    min_age: float = 0
    # range of peak magnitudes for submission
    min_peak_mag: float = 19.5
    max_peak_mag: float = 13
    # Reported detections in at least this many filters
    min_n_filters: int = 1
    # Minimal galactic latitide
    min_gal_lat: float = 14
    # reject alert if ssdistnr smaller than this value for any pp
    ssdistnr_max: float = 1
    # reject alert if PS1 star for any pp
    ps1_sgveto_rad: float = 1
    ps1_sgveto_sgth: float = 0.8
    # Minimal median RB.
    rb_minmed: float = 0.3
    drb_minmed: float = 0.95  # if drb found!
    # Try to reject likely CV through rejecting objects that quickly get very bright
    cut_fastrise: bool = True
    # Require each PP to have a magpsf lower than the diffmaglim
    require_lowerthanlim: bool = True

    # Cut to apply to all the photopoints in the light curve.
    # This will affect most operations, i.e. evaluating the position,
    # computing number of detections ecc.
    lc_filters: list[dict] = [
        {"attribute": "sharpnr", "operator": ">=", "value": -10.15},
        {"attribute": "programid", "operator": "==", "value": 1},
        {"attribute": "magfromlim", "operator": ">", "value": 0},
    ]

    # parameters for adding remarks to AT reports
    # Tag objects this close to SDSS galaxies as nuclear. Use negative to disable
    nuclear_dist: float = -1.0
    # Required distance to match with aav catalog. TODO: move?
    aav_dist: float = 1.0
    # (sigma!) if GAIA match is noisier than this, add a remark
    max_gaia_noise: float = 2.0

    # Implicitly we still assume we are submitting ZTF data
    ztf_tns_at: dict = {  # Default values to tag ZTF detections / ulims
        "flux_units": "1",
        "instrument_value": "196",
        "exptime": "30",
        "Observer": "Robot",
    }

    def inspect_catalog(self, cat_res: dict[str, Any]) -> bool:
        """
        Verify whether any catalog matching criteria prevents submission.

        """

        # Check that we got any catalogmatching results (that it was run)
        if self.require_catalogmatch and len(cat_res) == 0:
            self.logger.debug("no T2CATALOGMATCH results")
            return False

        # check that you have positive match in all of the necessary cataslogs:
        for needed_cat in self.needed_catalogs:
            if not cat_res.get(needed_cat, False):
                self.logger.debug(
                    "no T2CATALOGMATCH results for %s" % needed_cat,
                    extra={"catalog_matches": cat_res},
                )
                return False

        nedz = cat_res.get("NEDz", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if (nedz and not (self.min_redshift < nedz["z"] < self.max_redshift)) or (
            sdss_spec and not (self.min_redshift < sdss_spec["z"] < self.max_redshift)
        ):
            self.logger.debug(
                "transient z above limit.",
                extra={
                    "max_z": self.max_redshift,
                    "SDSS_spec": sdss_spec,
                    "NEDz": nedz,
                },
            )
            return False

        # another battle in the endless war against stars.
        # here we define a dict to treat each catalog in the same way
        star_filters: dict[str, dict[str, Any]] = {
            "SDSSDR10": {"class_col": "type", "star_val": 6},
            "LAMOSTDr4": {"class_col": "class", "star_val": "STAR"},
        }
        for cat_name, sfilter in star_filters.items():
            cat = cat_res.get(cat_name, False)
            cname, sval = sfilter["class_col"], sfilter["star_val"]
            if cat and cat[cname] == sval and cat["dist2transient"] < self.start_dist:
                self.logger.debug(
                    "transient matched with star in catalog.",
                    extra={"cat_name": cat_name, "cat_res": cat},
                )
                return False

        # cut matches with variable star catalog
        aavsovsx = cat_res.get("AAVSOVSX", False)
        if aavsovsx and aavsovsx["dist2transient"] < self.start_dist:
            self.logger.info("transient too close to AAVSOVSX sorce", extra=aavsovsx)
            return False

        # cut away bright stars. TODO: this considers just the closest matches...
        # as well as fast moving stars
        gaia_dr2 = cat_res.get("GAIADR2", None)
        if (
            gaia_dr2
            and gaia_dr2["Mag_G"] > 0
            and gaia_dr2["Mag_G"] < self.max_gaia_neighbour_gmag
        ):
            self.logger.debug("transient close to bright GAIA source", extra=gaia_dr2)
            return False

        # TODO: Check for fast moving stars

        # congratulation catalog, you made it!
        return True

    def inspect_lc(self, lc: LightCurve) -> bool:
        """
        Verify whether the transient lightcurve fulfill criteria for submission.

        """

        # apply cut on history: consider photophoints which are sharp enough
        if not (pps := lc.get_photopoints(filters=self.lc_filters)):
            return False

        # Current filters cannot sort two attributes
        if self.require_lowerthanlim:
            pps = [pp for pp in pps if pp["body"]["magpsf"] < pp["body"]["diffmaglim"]]

        # cut on number of detection
        if len(pps) < self.min_ndet:
            return False

        # cut on number of filters
        used_filters = set([pp["body"]["fid"] for pp in pps])
        if len(used_filters) < self.min_n_filters:
            return False

        # cut on range of peak magnitude
        mags = [pp["body"]["magpsf"] for pp in pps]
        peak_mag = min(mags)
        if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
            return False

        # cut on age
        jds = [pp["body"]["jd"] for pp in pps]
        most_recent_detection, first_detection = max(jds), min(jds)
        age = most_recent_detection - first_detection
        if age > self.max_age or age < self.min_age:
            return False

        # cut on galactic coordinates
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            return False

        coordinates = SkyCoord(ra, dec, unit="deg")
        b = coordinates.galactic.b.deg
        if abs(b) < self.min_gal_lat:
            return False

        # cut on number of detection after last SIGNIFICANT UL or r
        if ulims := lc.get_upperlimits(
            filters={
                "attribute": "diffmaglim",
                "operator": ">=",
                "value": self.max_maglim,
            }
        ):
            last_ulim = sorted(ulims, key=lambda x: x["body"]["jd"])[-1]
            pps_after_ndet = (
                lc.get_photopoints(
                    filters=self.lc_filters
                    + [
                        {
                            "attribute": "jd",
                            "operator": ">=",
                            "value": last_ulim["body"]["jd"],
                        }
                    ]
                )
                or []
            )
            # Can this work? - It does not seem like it
            # filters = self.lc_filters + [{'attribute': 'jd', 'operator': '>=', 'value': last_ulim["body"]['jd']}, {'attribute': 'magpsf', 'operator': '<', 'attribute': 'diffmaglim'}]
            # Current filters cannot sort two attributes
            if self.require_lowerthanlim:
                pps_after_ndet = [
                    pp
                    for pp in pps_after_ndet
                    if pp["body"]["magpsf"] < pp["body"]["diffmaglim"]
                ]

            if len(pps_after_ndet) < self.min_ndet_postul:
                return False

            # Check mag increase per time range for first detections
            first_pp_afterUL = sorted(pps_after_ndet, key=lambda x: x["body"]["jd"])[0]
            # This is only a relevant comparison if the obs after last significant UL is also first detection
            if self.cut_fastrise and first_pp_afterUL["body"]["jd"] == first_detection:
                delta_t = first_pp_afterUL["body"]["jd"] - last_ulim["body"]["jd"]
                delta_m = (
                    -first_pp_afterUL["body"]["magpsf"]
                    + last_ulim["body"]["diffmaglim"]
                )

                if delta_t < 3.5 and delta_m > 3:
                    self.logger.info(
                        "Likely CV", extra={"deltaT": delta_t, "deltaM": delta_m}
                    )
                    return False

        # cut on distance to closest solar system object
        # TODO: how to make this check: ('0.0' in list(phot["ssdistnr"])
        ssdist = np.array(
            [pp["body"]["ssdistnr"] for pp in pps if "ssdistnr" in pp["body"]]
        )
        ssdist[ssdist == None] = -999  # noqa: E711

        close_to_sso = np.logical_and(ssdist < self.ssdistnr_max, ssdist > 0)
        if np.any(close_to_sso):
            self.logger.info(
                "transient too close to solar system object",
                extra={"ssdistnr": ssdist.tolist()},
            )
            return False

        # check PS1 sg for the full alert history
        # Note that we for this check do *not* use the lightcurve filter criteria
        # TODO: Evaluate whether we should use the filters, and do a check for sufficient number of datapoints remaining
        # distpsnr1, sgscore1 = zip(*lc.get_tuples('distpsnr1', 'sgscore1', filters=self.lc_filters))
        if tups := lc.get_tuples("distpsnr1", "sgscore1"):
            distpsnr1, sgscore1 = zip(*tups)
        else:
            return False
        is_ps1_star = np.logical_and(
            np.array(distpsnr1) < self.ps1_sgveto_rad,
            np.array(sgscore1) > self.ps1_sgveto_sgth,
        )
        if np.any(is_ps1_star):
            self.logger.info(
                "transient below PS1 SG cut for at least one pp.",
                extra={"distpsnr1": distpsnr1, "sgscore1": sgscore1},
            )
            return False

        # cut on median RB score
        rbs = [pp["body"]["rb"] for pp in pps]
        if np.median(rbs) < self.rb_minmed:
            return False

        # cut on median dRB score
        drbs = [pp["body"]["drb"] for pp in pps if "drb" in pp["body"].keys()]
        if len(drbs) > 0 and np.median(drbs) < self.drb_minmed:
            return False

        # congratulation Lightcurve, you made it!
        return True

    def get_catalog_remarks(
        self, lc: LightCurve, cat_res: dict[str, Any]
    ) -> None | dict[str, Any]:
        """
        Look through catalogs for remarks to be added to report.
        """

        # Start building dict with remarks
        remarks: dict[str, Any] = {"remarks": ""}

        # Check redshift
        nedz = cat_res.get("NEDz", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if sdss_spec:
            remarks["remarks"] = (
                remarks["remarks"] + "SDSS spec-z %.3f. " % (sdss_spec["z"])
            )
        elif nedz:
            remarks["remarks"] = remarks["remarks"] + "NED z %.3f. " % (nedz["z"])

        # tag AGNs
        milliquas = cat_res.get("milliquas", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if (
            milliquas
            and milliquas["redshift"] is not None
            and milliquas["redshift"] > 0
        ) or (sdss_spec and sdss_spec["bptclass"] in [4, 5]):
            remarks["remarks"] = (
                remarks["remarks"] + "Known SDSS and/or MILLIQUAS QSO/AGN. "
            )
            remarks["at_type"] = 3

        # tag nuclear
        sdss_dr10 = cat_res.get("SDSSDR10", False)
        if (
            sdss_dr10
            and sdss_dr10["type"] == 3
            and sdss_dr10["dist2transient"] < self.nuclear_dist
        ):
            remarks["remarks"] = (
                remarks["remarks"] + "Close to core of SDSS DR10 galaxy. "
            )
            remarks["at_type"] = 4

        # tag noisy gaia
        if tups := lc.get_tuples("distpsnr1", "sgscore1", filters=self.lc_filters):
            distpsnr1, sgscore1 = zip(*tups)
            galaxylike_ps1 = np.logical_and(
                np.array(distpsnr1) < 1.5, np.array(sgscore1) < 0.5
            )
            gaia_dr2 = cat_res.get("GAIADR2", False)
            nedz = cat_res.get("NEDz", False)
            if (
                (
                    gaia_dr2
                    and gaia_dr2["ExcessNoise"] > self.max_gaia_noise
                    and gaia_dr2["dist2transient"] < 1
                )
                and (nedz and not (nedz["z"] > 0.01 and nedz["dist2transient"] < 1))
                and (  # if it's extragalactic
                    sdss_dr10
                    and not (sdss_dr10["type"] == 3 and sdss_dr10["dist2transient"] < 3)
                )
                and (  # and if it's not a galaxy
                    not np.any(galaxylike_ps1)
                )  # TODO: check the logic
            ):
                remarks["remarks"] = (
                    remarks["remarks"]
                    + "Significant noise in Gaia DR2 - variable star cannot be excluded. "
                )

        if len(remarks["remarks"]) == 0:
            return None
        else:
            return remarks

    def get_lightcurve_info(self, lc: LightCurve) -> None | dict[str, Any]:
        """
        Collect the data needed for the atreport. Return None in case
        you have to skip this transient for some reason.
        """

        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            return None

        # Start defining AT dict: name and position
        atdict: dict[str, Any] = {}
        # atdict.update(self.base_at_dict)
        atdict["ra"] = {"value": ra, "error": 1.0, "units": "arcsec"}
        atdict["dec"] = {"value": dec, "error": 1.0, "units": "arcsec"}

        # Add information on the latest SIGNIFICANT non detection.
        last_non_obs = 0
        if ulims := lc.get_upperlimits(
            filters={
                "attribute": "diffmaglim",
                "operator": ">=",
                "value": self.max_maglim,
            }
        ):
            last_ulim = sorted(ulims, key=lambda x: x["body"]["jd"])[-1]
            last_non_obs = last_ulim["body"]["jd"]
            filter_name = TNSFILTERID.get(last_ulim["body"]["fid"])
            atdict["non_detection"] = {
                "obsdate": last_ulim["body"]["jd"],
                "limiting_flux": last_ulim["body"]["diffmaglim"],
                "filter_value": filter_name,
            }
        else:
            atdict["non_detection"] = {
                "archiveid": "0",
                "archival_remarks": "ZTF non-detection limits not available",
            }

        atdict["non_detection"].update(self.ztf_tns_at)  # Add the default ZTF values

        # now add info on photometric detections: consider only candidates which
        # have some consecutive detection after the last ulim
        if pps := lc.get_photopoints(
            filters=self.lc_filters
            + [{"attribute": "jd", "operator": ">=", "value": last_non_obs}]
        ):
            # Lets create a few photometry points: TODO: should they be the latest or the first?
            atdict["photometry"] = {"photometry_group": {}}
            atdict["discovery_datetime"] = 10**30
            for ipp, pp in enumerate(pps[: self.nphot_submit]):
                photdict = {  # TODO: do we need to round the numerical values?
                    "obsdate": pp["body"]["jd"],
                    "flux": float("{:.2f}".format(pp["body"]["magpsf"])),
                    "flux_error": float("{:.2f}".format(pp["body"]["sigmapsf"])),
                    "limiting_flux": float("{:.2f}".format(pp["body"]["diffmaglim"])),
                    "filter_value": TNSFILTERID.get(pp["body"]["fid"]),
                }
                if pp["body"]["jd"] < atdict["discovery_datetime"]:
                    atdict["discovery_datetime"] = pp["body"]["jd"]
                photdict.update(self.ztf_tns_at)
                atdict["photometry"]["photometry_group"][str(ipp)] = photdict

        return atdict

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(
        self, light_curve: LightCurve, t2_views: Sequence[T2DocView]
    ) -> UBson | UnitResult:
        """

        Evaluate whether a transient passes thresholds for submission to TNS.
        If so, collects the required submission info.
        This unit does *not* verify whether a transient already exists in TNS.

        Parameters
        -----------
        light_curve: "ampel.view.LightCurve" instance.
        See the LightCurve docstring for more info.

        t2_records: List of T2Records

        Returns
        -------
        dict

        Warning: this dict does *not* contain the internal name (atdict['internal_name']) required.
        This will have to be added by the T3. Similarly, the base atdict entries are not included.

        In further changes, the photometry_group keys are now str instead of int (to be able to store)
        If needed by TNS, they need to be converted back at T3
        """

        self.logger.debug("starting eval from %s" % (t2_views))

        # i. Check whether the lightcurve passes selection criteria
        if not self.inspect_lc(light_curve):
            return {"tns_candidate": False, "tns_eval": "Poor lightcurve."}

        # ii. Check the catalog matching criteria
        t2_cat_match = t2_views[0]
        assert t2_cat_match.unit in (dep.unit for dep in self.t2_dependency)
        catalog_result = t2_cat_match.get_payload()
        assert isinstance(catalog_result, dict)
        if not self.inspect_catalog(catalog_result):
            return {"tns_candidate": False, "tns_eval": "Catalog match rejection."}

        # iii. Collect information for submission

        # These methods repeat a lot of calculations from above, but are kept separate
        # to easier override selection method
        atdict = self.get_lightcurve_info(light_curve)
        if atdict is None:
            return {
                "tns_candidate": False,
                "tns_eval": "Passes criteria, fails in info collection.",
            }

        catremarks = self.get_catalog_remarks(light_curve, catalog_result)
        if catremarks is not None:
            atdict.update(catremarks)

        t2result = {"tns_candidate": True, "tns_eval": "Good", "atdict": atdict}

        return t2result
