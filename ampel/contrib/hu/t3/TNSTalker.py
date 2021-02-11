#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t3/TNSTalker.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 17.11.2018
# Last Modified Date: 04.09.2029
# Last Modified By  : Jakob van Santen <jakob.van.santen@desy.de>

import re
from itertools import islice
from typing import (
    Any,
    Dict,
    Generator,
    Iterable,
    List,
    Optional,
    Tuple,
    TYPE_CHECKING,
)

import numpy as np
from astropy.coordinates import SkyCoord

from ampel.abstract.AbsT3Unit import AbsT3Unit
from ampel.contrib.hu.t3.ampel_tns import (
    get_tnsname,
    sendTNSreports,
    TNSFILTERID,
    tnsInternal,
)
from ampel.model.Secret import Secret
from ampel.struct.JournalTweak import JournalTweak
from ampel.type import StockId
from ampel.view.TransientView import TransientView
from ampel.ztf.utils import to_ztf_id

if TYPE_CHECKING:
    from ampel.content.JournalRecord import JournalRecord
    from ampel.protocol.LoggerProtocol import LoggerProtocol


def chunks(l: Iterable, n: int) -> Generator[List, None, None]:
    source = iter(l)
    while True:
        chunk = list(islice(source, n))
        yield chunk
        if len(chunk) < n:
            break


# get the science records for the catalog match
def get_catalogmatch_srecs(
    tran_view: TransientView, logger: "LoggerProtocol"
) -> Dict[str, Any]:
    if cat_res := [
        record
        for record in (tran_view.get_t2_records(unit_id="CATALOGMATCH") or [])
        if record["body"] is not None and record["body"][-1]["result"]
    ]:
        if (body := cat_res[-1]["body"]) is not None:
            return body[-1]["result"]

    logger.info("NO CATALOG MATCH FOR THIS TRANSIENT")
    return {}


class TNSTalker(AbsT3Unit):
    """
    Get TNS name if existing, and submit selected candidates.
    If submit_tns is true, candidates fulfilling the criteria will be sent to the TNS if:
    - They are not known to the TNS OR
    - They are registered by TNS but under a non-ZTF internal name AND resubmit_tns_nonztf set to True OR
    - They are registered by TNS under a ZTF name AND resubmit_tns_ztf is set to True

    if sandbox is set to True it will try to submit candidates to the TNS sandbox, but this API has been unstable
    and might not work properly.
    """

    version = 0.1

    # TNS config

    # Bot api key frm TNS
    tns_api_key: Optional[Secret] = None
    # Check for TNS for names even if internal name is known
    get_tns_force: bool = False
    # Submit candidates passing criteria (False gives you a 'dry run')
    submit_tns: bool = True
    # Submit all candidates we have a note in the Journal that we submitted this. Overrides the resubmit entries!!
    submit_unless_journal: bool = False
    # Resubmit candidate submitted w/o the same ZTF internal ID
    resubmit_tns_nonztf: bool = True
    # Resubmit candidates even if they have been added with this name before
    resubmit_tns_ztf: bool = False

    # Submit to TNS sandbox only
    sandbox: bool = True
    # weather journal will go to separate collection.
    ext_journal: bool = True

    # AT report config
    base_at_dict: Dict = {
        "reporting_group_id": "48",
        "discovery_data_source_id": "48",
        "reporter": "J. Nordin, V. Brinnel, M. Giomi, J. van Santen (HU Berlin), A. Gal-Yam, O. Yaron, S. Schulze (Weizmann) on behalf of ZTF",
        "at_type": "1",
    }
    ztf_tns_at: Dict = {  # Default values to tag ZTF detections / ulims
        "flux_units": "1",
        "instrument_value": "196",
        "exptime": "30",
        "Observer": "Robot",
    }
    baseremark: str = "See arXiv:1904.05922 for selection criteria."
    # Limiting magnitude to consider upper limits as 'significant'
    max_maglim: float = 19.5
    # Number of photometric detection we include in the TNS AT report
    nphot_submit: int = 2

    # cuts on T2 catalogs
    # reject candidates if they don't have matching in this list of T2CATALOGMATCH catalogs
    needed_catalogs: List[str] = []
    require_catalogmatch: bool = True
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
    max_age: float = 5
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
    # Try to reject likely CV through rejecting objects that quickly get very bright
    cut_fastrise: bool = True
    # Require each PP to have a magpsf lower than the diffmaglim
    require_lowerthanlim: bool = True

    # Cut to apply to all the photopoints in the light curve.
    # This will affect most operations, i.e. evaluating the position,
    # computing number of detections ecc.
    lc_filters: List[Dict] = [
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

    slack_token: Optional[Secret] = None
    slack_channel = "#ztf_tns"
    slack_username = "AMPEL"
    # if you have more than this # of reports, send different files
    max_slackmsg_size = 200

    def search_journal_tns(
        self, tran_view: TransientView
    ) -> Tuple[Optional[str], List[str]]:
        """
        Look through the journal for a TNS name.
        Assumes journal entries came from this unit, that the TNS name is saved as "tnsName"
        and internal names as "tnsInternal"
        """
        tns_name, tns_internals = None, []

        def select(entry: "JournalRecord") -> bool:
            return bool(
                (entry["extra"] is not None and ("tnsInternal" in entry["extra"]))
                and entry["unit"]
                and entry["unit"] == self.__class__.__name__
            )

        if jentries := tran_view.get_journal_entries(tier=3, filter_func=select):
            if jentries[-1]["extra"] is not None:
                tns_name = jentries[-1]["extra"].get("tnsName", None)
            tns_internals = [
                entry["extra"].get("tnsInternal", None)
                for entry in jentries
                if entry["extra"] is not None
            ]

        self.logger.info(
            "Journal search",
            extra={
                "tranId": tran_view.id,
                "tnsName": tns_name,
                "tnsInternals": tns_internals,
            },
        )

        return tns_name, tns_internals

    def search_journal_submitted(self, tran_view: TransientView) -> bool:
        """
        Look through the journal for whether this sender submitted this to TNS.
        Assumes journal entries came from this unit, that the TNS name is saved as "tnsName"
        and tnsSender stores the api key used ('tnsSender': self.tns_api_key')
        """

        def select(entry: "JournalRecord") -> bool:
            return bool(
                (
                    entry["extra"] is not None
                    and (entry["extra"].get("tnsSender") == self.tns_api_key)
                    and "tnsSubmitResult" in entry["extra"]
                )
                and entry["unit"]
                and entry["unit"] == self.__class__.__name__
            )

        # Find the latest tns name (skipping previous)
        if tran_view.get_journal_entries(
            tier=3,
            filter_func=select,
        ):
            self.logger.info("TNS submitted", extra={"tnsSender": self.tns_api_key})
            return True
        else:
            self.logger.info("Not TNS submitted", extra={"tnsSender": self.tns_api_key})
            return False

    def _query_tns_names(self, tran_view: TransientView) -> Tuple[Optional[str], List]:
        """
        query the TNS for names and internals at the position
        of the transient.
        """
        # query the TNS for transient at this position. Note that we check the real TNS for names for compatibility...
        if tran_view.lightcurve:
            lc = tran_view.lightcurve[-1]
            if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
                ra, dec = pos
            else:
                return None, []
        else:
            return None, []

        tns_name, tns_internal = get_tnsname(
            ra=ra,
            dec=dec,
            api_key=self.tns_api_key,
            logger=self.logger,
            sandbox=False,
        )

        # Skip the AT SN prefix if present
        if tns_name is not None:
            tns_name = re.sub("^AT", "", tns_name)
            tns_name = re.sub("^SN", "", tns_name)

        # be nice and then go
        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info(
            "looking for TNS name in the TNS.",
            extra={
                "ZTFname": ztf_name,
                "ra": ra,
                "dec": dec,
                "tnsName": tns_name,
                "tnsInternals": [tns_internal],
            },
        )
        return tns_name, [tns_internal]

    def _find_tns_tran_names(
        self, tran_view: TransientView
    ) -> Tuple[Optional[str], List[str]]:
        """
        search for TNS name in tran_view.tran_names. If found,
        look in the TNS for internal names and return them
        """

        # First, look if we already registered a name
        tns_name, tns_internals = None, []
        names: List[str] = (
            [str(name) for name in (tran_view.stock["name"] or [])]
            if tran_view.stock
            else []
        )
        for tname in names:

            if "TNS" in tname and (not self.get_tns_force):
                self.logger.info(
                    "found TNS name in tran_names.",
                    extra={"TNSname": tname, "TransNames": names},
                )
                # as TNS to give you the internal names.
                # we remove the 'TNS' part of the name, since this has been
                # added by the TNSMatcher T3, plus we skip the prefix
                # We here assume that the AT/SN suffix is cut
                tns_name = tname.replace("TNS", "")
                # Not using sandbox (only checking wrt to full system). This method still only returns a single internal name of the detecting survey
                tns_internal_single, runstatus = tnsInternal(
                    tns_name, api_key=self.tns_api_key, sandbox=False
                )
                if not tns_internal_single is None:
                    # Even if this is a single name, it comes as a list
                    tns_internals.append(tns_internal_single)

        # be nice with the logging
        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info(
            "looked for TNS name in self.tran_names",
            extra={
                "ZTFname": ztf_name,
                "tnsName": tns_name,
                "tnsInternals": tns_internals,
                "TransNames": names,
            },
        )

        # if you make it till here, no match was found
        return tns_name, tns_internals

    def find_tns_name(
        self, tran_view: TransientView
    ) -> Tuple[Optional[str], List[str], Optional[JournalTweak]]:
        """
        extensive search for TNS names in:
        - tran_view.tran_names (if added by TNSMatcher)
        - the journal of tran_view (if added by this T3)
        - the TNS itself (if no name can be found with the above)

        Returns:
        --------
            tns_name, tns_internals, jup: tns_name, tns_internal, and journal update
        """

        ztf_name = to_ztf_id(tran_view.id)
        self.logger.info("looking for TNS name", extra={"ZTFname": ztf_name})

        # first we look in the journal, this is the cheapest option. If we have
        # a valid name from the journal and if you do not want to look again in
        # the TNS, we are fine. NOTE: in this case you don't return a journal update.
        tns_name, tns_internals = self.search_journal_tns(tran_view)
        self.logger.debug("Found tns name in journal: %s" % (tns_name))
        if (not tns_name is None) and (not self.get_tns_force):
            return tns_name, tns_internals, None

        # second option in case there is no TNS name in the journal: go and look in tran_names
        # and if you don't find any, go and ask TNS again.
        tns_name_new, tns_internals_new = self._find_tns_tran_names(tran_view)
        self.logger.debug(
            "Find tns names added to the ampel name list: %s internal %s"
            % (tns_name_new, tns_internals_new)
        )
        if tns_name_new is None:
            tns_name_new, tns_internals_new = self._query_tns_names(tran_view)
            self.logger.debug(
                "Proper check of tns done, found name %s" % (tns_name_new)
            )

        # now, it is possible (if you set self.get_tns_force) that the
        # new TNS name is different from the one we had in the journal. We always
        # use the most recent one. In this case we also create a JournalUpdate
        jup = None
        if not tns_name_new is None:

            # what happen if you have a new name that is different from the old one?
            if tns_name is not None and not tns_name == tns_name_new:
                self.logger.info(
                    "Adding new TNS name",
                    extra={"tnsOld": tns_name, "tnsNew": tns_name_new},
                )

                # create content of journal entry. Eventually
                # update the list with the new internal names if any are found
                jcontent = {"tnsName": tns_name_new}
                if tns_internals_new is not None:
                    tns_internals.extend(tns_internals_new)
                    for tns_int in tns_internals_new:
                        jcontent.update({"tnsInternal": tns_int})

                # create a journalUpdate and update the tns_name as well. TODO: check with JNo
                jup = JournalTweak(extra=jcontent)
                tns_name = tns_name_new

            elif tns_name is None:
                # Set the new name
                self.logger.info(
                    "Adding first TNS name", extra={"tnsNew": tns_name_new}
                )

                # create content of journal entry. Eventually
                # update the list with the new internal names if any are found
                jcontent = {"tnsName": tns_name_new}
                if tns_internals_new is not None:
                    tns_internals.extend(tns_internals_new)
                    for tns_int in tns_internals_new:
                        jcontent.update({"tnsInternal": tns_int})

                # create a journalUpdate and update the tns_name as well. TODO: check with JNo
                jup = JournalTweak(extra=jcontent)
                tns_name = tns_name_new
                # tns_internals = tns_internals_new

        # bye!
        return tns_name, tns_internals, jup

    def accept_tview(self, tran_view: TransientView) -> bool:
        """
        decide weather or not this transient is worth submitting to TNS
        or not.

        NOTE that even if many of these cuts could defined passed directly to
        the task/job config, some of them still require relatively non trivial
        computation (e.g. 'age' of the transient). This makes this selection method
        necessary.

        FOR SAKE OF SIMPLICITY during the first period of TNS submission,
        we here declare that all the transient selection logic should be
        implemented in this method, even rejection based on catalog matching
        that could be done by the select config param.
        """

        # get the latest light curve
        assert tran_view.lightcurve
        lc = tran_view.lightcurve[-1]

        # apply cut on history: consider photophoints which are sharp enough
        if not (pps := lc.get_photopoints(filters=self.lc_filters)):
            return False
        # Current filters cannot sort two attributes
        if self.require_lowerthanlim:
            pps = [pp for pp in pps if pp["body"]["magpsf"] < pp["body"]["diffmaglim"]]
        self.logger.info("%d photop. passed filter %s" % (len(pps), self.lc_filters))

        # cut on number of detection
        if len(pps) < self.min_ndet:
            self.logger.info(
                ", not enough detections: got %d, required %d"
                % (len(pps), self.min_ndet)
            )
            return False

        # cut on number of filters
        used_filters = set([pp["body"]["fid"] for pp in pps])
        if len(used_filters) < self.min_n_filters:
            self.logger.info(
                "requested detections in more than %d bands, got: %d"
                % (self.min_n_filters, len(used_filters))
            )
            return False

        # cut on range of peak magnitude
        mags = [pp["body"]["magpsf"] for pp in pps]
        peak_mag = min(mags)
        if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
            self.logger.info(
                "peak magnitude of %.2f outside of range [%.2f, %.2f]"
                % (peak_mag, self.min_peak_mag, self.max_peak_mag)
            )
            return False

        # cut on age
        jds = [pp["body"]["jd"] for pp in pps]
        most_recent_detection, first_detection = max(jds), min(jds)
        # age = Time.now().jd - min(jds)
        age = most_recent_detection - first_detection
        if age > self.max_age or age < self.min_age:
            self.logger.info(
                "age of %.2f days outside of range [%.2f, %.2f]"
                % (age, self.min_age, self.max_age)
            )
            return False

        # cut on galactic coordinates
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            return False
        coordinates = SkyCoord(ra, dec, unit="deg")
        b = coordinates.galactic.b.deg
        if abs(b) < self.min_gal_lat:
            self.logger.info(
                "transient at b=%.2f too close to galactic plane (cut at %.2f)"
                % (b, self.min_gal_lat)
            )
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
                self.logger.info(
                    "not enough consecutive detections after last significant UL.",
                    extra={"NDet": len(pps), "lastUlimJD": last_ulim["body"]["jd"]},
                )
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
        ssdist = np.array([pp["body"]["ssdistnr"] for pp in pps])
        ssdist[ssdist == None] = -999

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
            self.logger.info(
                "Median RB %below limit.",
                extra={
                    "median_rd": np.median(rbs),
                    "rb_minmed": self.rb_minmed,
                },
            )
            return False

        # -------------------------- #
        #    CUTS ON T2 RECORDS      #
        # -------------------------- #
        cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)

        # check that we got any catalogmatching results (that it was run)
        if self.require_catalogmatch and len(cat_res) == 0:
            self.logger.info("no T2CATALOGMATCH results")
            return False

        # check that you have positive match in all of the necessary cataslogs:
        for needed_cat in self.needed_catalogs:
            if not cat_res.get(needed_cat, False):
                self.logger.info(
                    "no T2CATALOGMATCH results for %s" % needed_cat,
                    extra={"catalog_matches": cat_res},
                )
                return False

        nedz = cat_res.get("NEDz", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if (nedz and not (self.min_redshift < nedz["z"] < self.max_redshift)) or (
            sdss_spec and not (self.min_redshift < sdss_spec["z"] < self.max_redshift)
        ):
            self.logger.info(
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
        star_filters: Dict[str, Dict[str, Any]] = {
            "SDSSDR10": {"class_col": "type", "star_val": 6},
            "LAMOSTDr4": {"class_col": "class", "star_val": "STAR"},
        }
        for cat_name, sfilter in star_filters.items():
            cat = cat_res.get(cat_name, False)
            cname, sval = sfilter["class_col"], sfilter["star_val"]
            if cat and cat[cname] == sval and cat["dist2transient"] < self.start_dist:
                self.logger.info(
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
        gaia_dr2 = cat_res.get("GAIADR2", None)
        if (
            gaia_dr2
            and gaia_dr2["Mag_G"] > 0
            and gaia_dr2["Mag_G"] < self.max_gaia_neighbour_gmag
        ):
            self.logger.info("transient close to bright GAIA source", extra=gaia_dr2)
            return False

        # congratulation TransientView, you made it!
        return True

    def add_atreport_remarks(
        self, tran_view: TransientView
    ) -> Optional[Dict[str, Any]]:
        """
        create additional remarks based, i.e. on catalog matching data
        """

        # TODO: check values for the cuts and and put them in the config
        # TODO: can remarks be combined? e.g. nucler + noisy?

        # get the science records for the catalog match
        cat_res = get_catalogmatch_srecs(tran_view, logger=self.logger)

        # tag AGNs
        milliquas = cat_res.get("milliquas", False)
        sdss_spec = cat_res.get("SDSS_spec", False)
        if (milliquas and milliquas["redshift"] > 0) or (
            sdss_spec and sdss_spec["bptclass"] in [4, 5]
        ):  # TODO: add distance cut?
            self.logger.info(
                "Transient is SDSS BPT or Milliquas AGN.",
                extra={
                    "tranId": tran_view.id,
                    "milliquas": milliquas,
                    "SDSS_spec": sdss_spec,
                },
            )
            return {"remarks": "Known SDSS and/or MILLIQUAS QSO/AGN. ", "at_type": 3}

        # tag nuclear
        sdss_dr10 = cat_res.get("SDSSDR10", False)
        if (
            sdss_dr10
            and sdss_dr10["type"] == 3
            and sdss_dr10["dist2transient"] < self.nuclear_dist
        ):
            self.logger.info(
                "Transient close to SDSS photometric galaxy - possibly nuclear",
                extra={"tranId": tran_view.id, "SDSSDR10": sdss_dr10},
            )
            return {"remarks": "Close to core of SDSS DR10 galaxy", "at_type": 4}

        # tag noisy gaia
        if tran_view.lightcurve and (
            tups := tran_view.lightcurve[-1].get_tuples(
                "distpsnr1", "sgscore1", filters=self.lc_filters
            )
        ):
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
                self.logger.info(
                    "Significant noise in Gaia DR2 - variable star cannot be excluded.",
                    extra={
                        "tranId": tran_view.id,
                        "GAIADR2": gaia_dr2,
                        "NEDz": nedz,
                        "SDSSDR10": sdss_dr10,
                    },
                )
                return {
                    "remarks": "Significant noise in Gaia DR2 - variable star cannot be excluded."
                }
        return None

    def create_atreport(self, tran_view: TransientView) -> Optional[Dict[str, Any]]:
        """
        Collect the data needed for the atreport. Return None in case
        you have to skip this transient for some reason.
        """

        self.logger.info("creating AT report for transient.")
        ztf_name = to_ztf_id(tran_view.id)
        if (
            tran_view.lightcurve
            and (lc := tran_view.lightcurve[-1])
            and (pos := lc.get_pos(ret="mean", filters=self.lc_filters))
        ):
            ra, dec = pos
        else:
            return None

        # Start defining AT dict: name and position
        atdict = {}
        atdict.update(self.base_at_dict)
        atdict["internal_name"] = ztf_name
        atdict["ra"] = {"value": ra, "error": 1.0, "units": "arcsec"}
        atdict["dec"] = {"value": dec, "error": 1.0, "units": "arcsec"}

        # Add information on the latest SIGNIFICANT non detection. TODO: check!
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
            atdict["discovery_datetime"] = 10 ** 30
            for ipp, pp in enumerate(pps[: self.nphot_submit]):
                photdict = {  # TODO: do we need to round the numerical values?
                    "obsdate": pp["body"]["jd"],
                    "flux": float("{0:.2f}".format(pp["body"]["magpsf"])),
                    "flux_error": float("{0:.2f}".format(pp["body"]["sigmapsf"])),
                    "limiting_flux": float("{0:.2f}".format(pp["body"]["diffmaglim"])),
                    "filter_value": TNSFILTERID.get(pp["body"]["fid"]),
                }
                if pp["body"]["jd"] < atdict["discovery_datetime"]:
                    atdict["discovery_datetime"] = pp["body"]["jd"]
                photdict.update(self.ztf_tns_at)
                atdict["photometry"]["photometry_group"][int(ipp)] = photdict

        # finally, add remarks based on catalogs adn return
        remarks = self.add_atreport_remarks(tran_view)
        if not remarks is None:
            atdict.update(remarks)
        if "remarks" in atdict.keys():
            atdict["remarks"] = "%s. %s." % (
                atdict["remarks"],
                self.baseremark,
            )
        else:
            atdict["remarks"] = self.baseremark
        return atdict

    def add(self, transients: Tuple[TransientView, ...]) -> Dict[StockId, JournalTweak]:
        """
        Loop through transients and check for TNS names and/or candidates to submit
        """

        if transients is None:
            self.logger.info("no transients for this task execution")
            return {}

        # select the transients
        transients_to_submit = [tv for tv in transients if self.accept_tview(tv)]
        self.logger.info(
            "of the %d transients presented to this task, %d passed selection criteria"
            % (len(transients), len(transients_to_submit))
        )

        # Will be saved to future journals
        journal_updates: Dict[StockId, JournalTweak] = {}
        # Reports to be sent, indexed by the transient view IDs (so that we can check in the replies)
        atreports: Dict[StockId, Dict[str, Any]] = {}
        for tran_view in transients_to_submit:

            ztf_name = to_ztf_id(tran_view.id)
            self.logger.info(
                "TNS check", extra={"tranId": tran_view.id, "ztfName": ztf_name}
            )
            self.logger.debug("TNS check for %s" % (ztf_name))

            # Simplest case to check. We wish to submit everything not noted as submitted
            if self.submit_unless_journal:

                if self.search_journal_submitted(tran_view):
                    # Note already submitted
                    self.logger.info("ztf submitted", extra={"ztfSubmitted": True})

                else:
                    # create AT report
                    if atreport := self.create_atreport(tran_view):
                        self.logger.info("Added to report list")
                        atreports[tran_view.id] = atreport

                continue

            # find the TNS name, either from the journal, from tran_names, or
            # from TNS itself. If new names are found, create a new JournalUpdate
            tns_name, tns_internals, jup = self.find_tns_name(tran_view)
            if not jup is None:
                journal_updates[tran_view.id] = jup
            self.logger.debug("TNS got %s internals %s" % (tns_name, tns_internals))

            if not tns_name is None:

                # Chech whether this ID has been submitted (note that we do not check
                # whether the same candidate was submitted as different ZTF name) and
                # depending on what's already on the TNS we can chose to submit or not
                is_ztfsubmitted = ztf_name in tns_internals
                if is_ztfsubmitted:
                    # Already registered under this name. Only submit if we explicitly configured to do this
                    if not self.resubmit_tns_ztf:
                        self.logger.info(
                            "ztf submitted",
                            extra={
                                "ztfSubmitted": is_ztfsubmitted,
                                "tnsInternals": tns_internals,
                            },
                        )
                        continue

                # Also allow for the option to not submit if someone (anyone) already did this. Not sure why this would be a good idea.
                if not is_ztfsubmitted and not self.resubmit_tns_nonztf:
                    self.logger.info(
                        "already in tns, skipping",
                        extra={
                            "ztfSubmitted": is_ztfsubmitted,
                            "tnsInternals": tns_internals,
                        },
                    )
                    continue

            # create AT report
            if atreport := self.create_atreport(tran_view):
                self.logger.info("Added to report list")
                atreports[tran_view.id] = atreport

        # TODO: we save the atreports to send them to the TNS.
        # This is just part of the tesing and will have to go away
        # 		atreports = {k: atreports[k] for k in list(atreports.keys())[:2]}
        self.atreports = atreports
        self.logger.info("collected %d AT reports to post" % len(atreports))

        # If we do not want to submit anything, or if there's nothing to submit
        if len(atreports) == 0 or (not self.submit_tns):
            self.logger.info(
                "submit_tns config parameter is False or there's nothing to submit",
                extra={
                    "n_reports": len(atreports),
                    "submit_tns": self.submit_tns,
                },
            )
            return journal_updates

        # atreports is now a dict with tran_id as keys and atreport as keys
        # what we need is a list of dicts with form {'at_report':atreport }
        # where an atreport is a dictionary with increasing integer as keys and atreports as values
        atreportlist = [
            {
                "at_report": {
                    i: report
                    for chunk in chunks(atreports.values(), 90)
                    for i, report in chunk
                }
            }
        ]
        tnsreplies = sendTNSreports(
            atreportlist,
            self.tns_api_key,
            self.logger,
            sandbox=self.sandbox,
        )

        # Now go and check and create journal updates for the cases where SN was added
        for tran_id in atreports.keys():
            ztf_name = to_ztf_id(tran_id)
            if not ztf_name in tnsreplies.keys():
                self.logger.info("No TNS add reply", extra={"tranId": tran_id})
                continue

            # Create new journal entry
            # TODO: do we want to add to the journal a failed TNS submit?
            jup = JournalTweak(
                extra={
                    "tnsName": tnsreplies[ztf_name][1]["TNSName"],
                    "tnsInternal": ztf_name,
                    "tnsSubmitresult": tnsreplies[ztf_name][0],
                    "tnsSender": self.tns_api_key,
                },
            )
            journal_updates[tran_view.id] = jup
        return journal_updates

    def done(self) -> None:
        self.logger.info("done running T3")

        if not hasattr(self, "atreports"):
            self.logger.info("No atreports collected.")
            return
        elif self.slack_token is None:
            return

        # TODO: to help debugging and verification, we post the collected atreports
        # to the slack, so that we can compare them with what JNo script is doing
        # ALL THE CONTENT OF THIS METHOD SHOULD GO AWAY AS SOON AS WE TRUST THIS T3
        self.logger.warn(
            "Posting collected ATreports to Slack. I'm still running as a test!"
        )

        import datetime
        import io
        import json

        from slack import WebClient
        from slack.errors import SlackClientError
        from slack.web.slack_reponse import SlackResponse

        sc = WebClient(token=self.slack_token.get())

        tstamp = datetime.datetime.today().strftime("%Y-%m-%d-%X")
        atlist = list(self.atreports.values())
        last = 0
        for ic, atrep in enumerate(chunks(atlist, self.max_slackmsg_size)):

            # add the atreport to a file
            self.logger.info("Posting chunk #%d" % ic)
            filename = "TNSTalker_DEBUG_%s_chunk%d.json" % (tstamp, ic)
            fbuffer = io.StringIO(filename)
            json.dump(atrep, fbuffer, indent=2)

            # upload the file with the at reports
            first = last
            last += len(atrep)
            msg = (
                "A total of %d atreports found by TNSTalker T3. Here's chunk #%d (reports from %d to %d)"
                % (len(self.atreports), ic, first, last)
            )
            api = sc.files_upload(
                channels=[self.slack_channel],
                title="TNSTalker_%s_chunk%d" % (tstamp, ic),
                initial_comment=msg,
                username=self.slack_username,
                as_user=False,
                filename=filename,
                filetype="javascript",
                file=fbuffer.getvalue(),
            )
            assert isinstance(api, SlackResponse)
            if not api["ok"]:
                raise SlackClientError(api["error"])

        self.logger.warn(
            f"DONE DEBUG Slack posting. Look at {self.slack_channel} for the results"
        )
