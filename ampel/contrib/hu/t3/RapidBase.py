#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : Ampel-contrib-HU/ampel/contrib/hu/t3/RapidBase.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 15.07.2019
# Last Modified Date: 06.02.2020
# Last Modified By  : vb <vbrinnel@physik.hu-berlin.de>

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from astropy.coordinates import Distance, SkyCoord
from astropy.cosmology import Planck15

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.base import abstractmethod
from ampel.model.Secret import Secret
from ampel.struct.JournalExtra import JournalExtra
from ampel.view.TransientView import TransientView
from ampel.ztf.utils import to_ampel_id, to_ztf_id


# get the science records for the catalog match
def get_catalogmatch_srecs(tran_view, logger):
    cat_res = tran_view.get_science_records(t2_class_name="CATALOGMATCH")
    if len(cat_res) == 0 or cat_res is None or cat_res[-1].get_results() is None:
        logger.info("NO CATALOG MATCH FOR THIS TRANSIENT")
        return {}
    return cat_res[-1].get_results()[-1]["output"]


class RapidBase(AbsT3Unit):
    """
    Select transients for rapid reactions. Intended as base class where the react method can be
    implemented as wished and a testreact method posts test reactions to Slack
    """

    # # weather journal will go to separate collection
    ext_journal: bool = True

    # Unless set, no full reaction will be triggered
    do_react: bool

    # # If set, will post trigger to slack
    do_testreact: bool
    slack_token: Optional[Secret]
    slack_channel: str = "#ztf_auto"
    slack_username: str = "AMPEL"

    # Cuts based on T2 catalog redshifts

    # Require a redshift max from a T2 output
    require_catalogmatch: bool = True
    # List of catalog-like output to search for redshift
    redshift_catalogs: List[str] = []
    # maximum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    max_redshift: float = 0.1
    # minimum redshift from T2 CATALOGMATCH catalogs (e.g. NEDz and SDSSspec)
    min_redshift: float = 0.001
    # max abs mag through peak mag and redshift from catalog mach (require both)
    max_absmag: float = -13
    # min abs mag through peak mag and redshift from catalog mach (require both)
    min_absmag: float = -17
    # arcsec, minimum distance to remove star matches to transient if found (eg in SDSSDR10)
    min_dist: float = 1.5
    # arcsec, maximum distance
    max_dist: float = 50
    # kpc, maximum distance
    max_kpc_dist: float = 999

    # Cut on alert properties

    # A candidate need to have at least this many detections
    min_ndet: int = 2
    min_ndet_postul: int = 2  # and if it has this minimum nr of detection after the last significant (max_maglim) UL.
    max_age: float = (
        3  # days, If a detection has an age older than this, skip (stars,age).
    )
    # Min age of detection history
    min_age: float = 0
    # range of peak magnitudes for submission
    min_peak_mag: float = 20
    max_peak_mag: float = 16
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
    # Minimal median RB.
    drb_minmed: float = 0.95
    # NOT IMPLEMENTED

    min_magrise: float = -20
    # Limiting magnitude to consider upper limits as 'significant'
    maglim_min: float = 19.5
    # A limiting magnitude max this time ago
    maglim_maxago: float = 2.5

    # Cut to apply to all the photopoints in the light curve.
    # This will affect most operations, i.e. evaluating the position,
    # computing number of detections ecc.
    lc_filters: List[Dict[str, Any]] = [
        {"attribute": "sharpnr", "operator": ">=", "value": -10.15},
        {"attribute": "magfromlim", "operator": ">", "value": 0},
    ]

    def post_init(self) -> None:

        self.name = "RapidBase"
        self.logger.info(f"Initialized T3 RapidBase instance {self.name}")

        # feedback
        for k in self.__annotations__:
            self.logger.info(f"Using {k}={getattr(self, k)}")

    def react(
        self, tran_view: TransientView, info: Dict[str, Any]
    ) -> Tuple[bool, Optional[JournalExtra]]:
        """
        Replace with react method adopted to particular facility or output
        """
        raise NotImplementedError("No real reaction implemented in RapidBase")
        return self.test_react(tran_view)

    def test_react(
        self, tran_view: TransientView, info: Dict[str, Any]
    ) -> Tuple[bool, Optional[JournalExtra]]:
        """ Trigger a test slack report """

        success = False

        if not self.slack_token:
            return False, None

        from slack import WebClient
        from slack.errors import SlackClientError
        from slack.web.slack_response import SlackResponse

        sc = WebClient(self.slack_token.get())
        assert isinstance(tran_view.id, int)
        ztf_name = to_ztf_id(tran_view.id)
        assert tran_view.lightcurve is not None
        lc = tran_view.lightcurve[-1]
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            raise ValueError("Light curve contains no points")
        msg = "Pancha says: Look up %s at RA %s DEC %s. Added info %s" % (
            ztf_name,
            ra,
            dec,
            info,
        )
        api = sc.chat_postMessage(
            channel=self.slack_channel,
            text=msg,
            username=self.slack_username,
            as_user=False,
        )
        assert isinstance(api, SlackResponse)
        if not api["ok"]:
            raise SlackClientError(api["error"])
        else:
            success = True

        description = "Sent SLACK msg"
        self.logger.info(description, extra={"channel": self.slack_channel})

        # Document what we did
        jcontent = {"t3unit": self.name, "reaction": description, "success": success}
        jup = JournalExtra(extra=jcontent)

        return success, jup

    def accept_tview(self, tran_view: TransientView) -> Optional[Dict[str, Any]]:
        """
        decide weather or not this transient is worth reacting to.

        NOTE that even if many of these cuts could defined passed directly to
        the task/job config, some of them still require relatively non trivial
        computation (e.g. 'age' of the transient). This makes this selection method
        necessary.
        """

        # We are lazy and create an info dict that can be included with the printout
        # should properly not be part of accept method
        info: Dict[str, Any] = {}

        # get the latest light curve
        assert tran_view.lightcurve is not None
        lc = tran_view.lightcurve[-1]

        # apply cut on history: consider photophoints which are sharp enough
        pps = lc.get_photopoints(filters=self.lc_filters)
        assert pps is not None
        self.logger.info("%d photop. passed filter %s" % (len(pps), self.lc_filters))
        # 		print("%d photop. passed filter %s" % (len(pps), self.lc_filters))

        # cut on number of detection
        if len(pps) < self.min_ndet:
            self.logger.info(
                "not enough detections: got %d, required %d" % (len(pps), self.min_ndet)
            )
            return None
        info["detections"] = len(pps)

        # cut on age
        jds = [pp["body"]["jd"] for pp in pps]
        most_recent_detection, first_detection = max(jds), min(jds)
        age = most_recent_detection - first_detection
        if age > self.max_age or age < self.min_age:
            self.logger.info(
                "age of %.2f days outside of range [%.2f, %.2f]"
                % (age, self.min_age, self.max_age)
            )
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
                filters=self.lc_filters
                + [{"attribute": "jd", "operator": ">=", "value": last_ulim_jd}]
            )
            # Check if there are enough positive detection after the last significant UL
            if (
                pps_after_ndet is not None
                and len(pps_after_ndet) < self.min_ndet_postul
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
                "requested detections in more than %d bands, got: %d"
                % (self.min_n_filters, len(used_filters))
            )
            return None

        # cut on range of peak magnitude
        mags = [pp["body"]["magpsf"] for pp in pps]
        peak_mag = min(mags)
        if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
            self.logger.info(
                "peak magnitude of %.2f outside of range [%.2f, %.2f]"
                % (peak_mag, self.min_peak_mag, self.max_peak_mag)
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
            if len(latest_pps) == 1:
                raise ValueError("Have assumed a unique last photopoint")
            info["latest_mag"] = latest_pps[0]["body"]["magpsf"]

        # we should here add a cut based on the mag rise per day (see submitRapid)

        # cut on galactic coordinates
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            raise ValueError("Light curve contains no points")
        coordinates = SkyCoord(ra, dec, unit="deg")
        b = coordinates.galactic.b.deg
        if abs(b) < self.min_gal_lat:
            self.logger.info(
                "transient at b=%.2f too close to galactic plane (cut at %.2f)"
                % (b, self.min_gal_lat)
            )
            return None
        info["ra"] = ra
        info["dec"] = dec

        # cut on distance to closest solar system object
        # TODO: how to make this check: ('0.0' in list(phot["ssdistnr"])
        ssdist = np.array([pp["body"]["ssdistnr"] for pp in pps])
        ssdist[ssdist is None] = -999
        # print (ssdist)

        close_to_sso = np.logical_and(ssdist < self.ssdistnr_max, ssdist > 0)
        if np.any(close_to_sso):
            self.logger.info(
                "transient too close to solar system object",
                extra={"ssdistnr": ssdist.tolist()},
            )
            return None

        # check PS1 sg for the full alert history
        # Note that we for this check do *not* use the lightcurve filter criteria
        # TODO: Evaluate whether we should use the filters, and do a check for sufficient number of datapoints remaining
        # 		print(to_ztf_id(tran_view.id))
        # 		print(lc)
        # 		print(lc.get_tuples('distpsnr1', 'sgscore1'))
        # 		print(lc.get_tuples('distpsnr1', 'sgscore1', filters=self.lc_filters))

        # 		distpsnr1, sgscore1 = zip(*lc.get_tuples('distpsnr1', 'sgscore1', filters=self.lc_filters))
        if psdata := lc.get_tuples("distpsnr1", "sgscore1"):
            distpsnr1, sgscore1 = zip(*psdata)
            is_ps1_star = np.logical_and(
                np.array(distpsnr1) < self.ps1_sgveto_rad,
                np.array(sgscore1) > self.ps1_sgveto_sgth,
            )
            if np.any(is_ps1_star):
                self.logger.info(
                    "transient below PS1 SG cut for at least one pp.",
                    extra={"distpsnr1": distpsnr1, "sgscore1": sgscore1},
                )
                return None
        else:
            self.logger.info("No PS1 check as no data found.")

        # cut on median RB and DRB score
        rbs = [pp["body"]["rb"] for pp in pps]
        if np.median(rbs) < self.rb_minmed:
            self.logger.info(
                "RB cut",
                extra={"median_rd": np.median(rbs), "rb_minmed": self.rb_minmed},
            )
            return None
        elif (len(rbs) == 0) and self.rb_minmed > 0:
            self.logger.info("No rb info for significant detection.")
            return None
        info["rb"] = np.median(rbs)

        # drb might not exist
        drbs = [pp["body"]["drb"] for pp in pps if "drb" in pp["body"]]
        if len(drbs) > 0 and np.median(drbs) < self.drb_minmed:
            self.logger.info(
                "DRB cut",
                extra={"median_drd": np.median(drbs), "drb_minmed": self.drb_minmed},
            )
            return None
        elif (len(drbs) == 0) and self.drb_minmed > 0:
            self.logger.info("No drb info for significant detection.")
            return None

        info["drb"] = np.median(drbs)

        # ---------------------------------------------------------------------#
        # CUTS ON T2 RECORDS                                                   #
        # ----------------------------------------------------------------------#

        # T2 Catalog matching
        cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)

        # check that we got any catalogmatching results (that it was run)
        if self.require_catalogmatch:

            if len(cat_res) == 0:
                self.logger.info("no T2CATALOGMATCH results")
                return None

            # Loop through listed catalogs for match
            zmatchs = []
            for catname in self.redshift_catalogs:
                catinfo = cat_res.get(catname, False)
                if (
                    catinfo
                    and (self.min_redshift < catinfo["z"] < self.max_redshift)
                    and (self.min_dist < catinfo["dist2transient"] < self.max_dist)
                ):
                    self.logger.info(
                        "z matched.",
                        extra={
                            "catalog": catname,
                            "z": catinfo["z"],
                            "dist": catinfo["dist2transient"],
                        },
                    )
                    # Calculate physical distance
                    dst_kpc = (
                        catinfo["dist2transient"]
                        * Planck15.kpc_proper_per_arcmin(catinfo["z"]).value
                        / 60.0
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
                self.logger.info("No z match.")
                return None

            # Determine physical distance

            # Determine absolute magnitue
            sndist = Distance(z=np.mean(zmatchs), cosmology=Planck15)
            absmag = info["peak_mag"] - sndist.distmod.value
            if not (self.min_absmag < absmag < self.max_absmag):
                self.logger.info("Not in absmag range.", extra={"absmag": absmag})
                # print('TEST z %.3f peakmag %.3f absmag %.3f' % (np.None(zmatchs), info['peak_mag'],absmag))
                return None
            info["absmag"] = absmag

        # tag AGNs
        milliquas = cat_res.get("milliquas", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if milliquas and milliquas["redshift"] > 0:
            info["milliAGN"] = True
        if sdss_spec and sdss_spec["bptclass"] in [4, 5]:
            info["sdssAGN"] = True

        # Potentially other checks on T2 results, eg photoz and lightcurve

        # congratulation TransientView, you made it!
        return info

    def add(self, transients):
        """
        Loop through transients and check for TNS names and/or candidates to submit
        """

        if transients is None:
            self.logger.info("no transients for this task execution")
            return []

        journal_updates = []
        # We will here loop through transients and react individually
        for tv in transients:
            matchinfo = self.accept_tview(tv)

            # Check sumission criteria
            if not matchinfo:
                continue

            self.logger.info("Passed reaction threshold", extra={"tranId": tv.id})

            # Ok, so we have a transient to react to
            if self.do_react:
                success, jup = self.react(tv, matchinfo)
                if jup is not None:
                    journal_updates.append(jup)

            # Otherwise, test
            if self.do_testreact:
                test_success, jup = self.test_react(tv, matchinfo)
                if jup is not None:
                    journal_updates.append(jup)

        return journal_updates

    def done(self):
        """ """
        # Should possibly do some accounting or verification
        self.logger.info("done running T3")
