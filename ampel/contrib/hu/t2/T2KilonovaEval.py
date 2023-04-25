#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2KilonovaEval.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                29.03.2023
# Last Modified Date:  21.04.2023
# Last Modified By:    ernstand@physik.hu-berlin.de

import numpy as np
from typing import Any, Literal, Union
from collections.abc import Sequence
from astropy.coordinates import Distance, SkyCoord
from astropy.cosmology import Planck15

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.view.LightCurve import LightCurve
from ampel.view.T2DocView import T2DocView
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.model.StateT2Dependency import StateT2Dependency

#from nuztf.nuztf import cat_match

class T2KilonovaEval(AbsTiedLightCurveT2Unit):
    """
    Evaluate whether a transient fulfills criteria for being a potential
    kilonova-like event.

    Could include evaluations based on (if present):
    - Lightcurve directly.
    - Redshift / distance to core.
    - SNguess.
    - Healpix map probability.
    - Sncosmo fits.
    - Parsnip fits.
    - Possis fits
    
    Will combine evaluation into a "kilonovaness" grade. 
    

    """

    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal[
    	"T2DigestRedshifts", 
    	"T2RunPossis", 
    	"T2PropagateStockInfo",
        "T2CatalogMatch"
    	]]]


    # Evaluation sections
    
    # Distance
    max_redshift: float = 0.05     # Max 
    min_redshift: float = 0.001
    min_dist: float = 1.5      # Min arcsec distance (remove core variables)
    max_dist: float = 50       # Max arcsec distance 
    max_kpc_dist: float = 999  # Max distance in kpc (using redshift)
    max_redshift_uncertainty: float = 999

    # Lightcurve (using redshift, so that needs to be there)
    min_absmag: float = -20 
    max_absmag: float = -12  #  
    min_ndet: int = 1
    min_ndet_postul: int = 0  # and if it has this minimum nr of detection after the last significant (max_maglim) UL.
    min_age: float = 0.01 # require 2 detections separated by 15 minutes
    max_age: float = 3.
    # Min age of detection history
    # range of peak magnitudes for submission
    min_peak_mag: float = 20
    max_peak_mag: float = 16
    # Reported detections in at least this many filters
    min_n_filters: int = 1
    # Require a detection in one of these filters (e.g. ZTF I-band more often spurious)
    det_filterids: list[int] = [1, 2, 3]   # default to any of them
    # Below are copied from filter - not sure is needed, but kept for legacy
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


    # Catalogs 
    # catalogs =["GLADEv23",
    #             "NEDz",
    #             "NEDz_extcats",
    #             "SDSS_spec",
    #             "LSPhotoZZou",
    #             "twoMPZ",
    #             "wiseScosPhotoz",
    #             "PS1_photoz",
    #             "CRTS_DR1",
    #             "milliquas",
    #             "GAIADR2",
    #             "SDSSDR10",
    #             "wise_color",
    #             "TNS"
    #             ]

    def inspect_ampelz(self, t2res: dict[str, Any]) -> None | dict[str, Any]:
        """
        Check whether Ampel Z data (from T2DigestRedshifts) fulfill criteria.
        """
        info = {'pass':0, 'rejects': []}
        
        if not t2res.get('ampel_z'):
            # No match
            criterium_name = "no ampelz match"
            info["rejects"].append(criterium_name)
            info["pass"] -= 10 # TODO maybe find a better way to punish this
            return info
        #info["pass"] += 1
            
        info['ampel_z'] = t2res['ampel_z']
        info['ampel_z_precision'] = t2res['group_z_precision']
        info['ampel_dist'] = t2res['ampel_dist']
        # Add score
        if (self.min_redshift < info["ampel_z"] < self.max_redshift):
            info['pass'] += 1
        else:
            criterium_name = "redshift"
            info['rejects'].append(criterium_name)
        if (self.min_dist < info["ampel_dist"] < self.max_dist):
            info['pass'] += 1
        else:
            criterium_name = "ampel_dist"
            info['rejects'].append(criterium_name)

        # Calculate physical distance
        info['dst_kpc'] = (
            info["ampel_dist"] *
                    Planck15.kpc_proper_per_arcmin(info["ampel_z"]).value / 60.0
                )
        if info['dst_kpc'] < self.max_kpc_dist:
            info['pass'] += 1
        else:
            criterium_name = "dst_kpc"
            info['rejects'].append(criterium_name)
        
        # Return collected info
        return info

    def inspect_possis(self, t2res: dict[str, Any]) -> None | dict[str, Any]:
        """
        Check whether a fit to T2RunPossis models look good.
        """
        info = {'pass':0, 'model':t2res['model_name'], 'rejects': []}
        if not t2res['z'] or not t2res['sncosmo_result']['success']:
            criterium_name = "no possis fits"
            info['rejects'].append(criterium_name)
            info["pass"] -= 10 # TODO maybe find a better way to punish this
            return info	# doesnt make sense to continue analysis if no values available
        #info["pass"]
        
        info['possis_abspeak'] = t2res['fit_metrics']['restpeak_model_absmag_B']
        info['possis_obspeak'] = t2res['fit_metrics']['obspeak_model_B']
        info['possis_chisq'] = t2res['sncosmo_result']['chisq']
        info['possis_ndof'] = t2res['sncosmo_result']['ndof']
        
        if info['possis_ndof']<0:
            criterium_name = "possis_ndof"
            info['rejects'].append(criterium_name)
            info["pass"] -= 10 # TODO maybe find a better way to punish this
            return info
        info["pass"] += 1
        
        if (self.min_absmag < info["possis_abspeak"] < self.max_absmag):
            info['pass'] += 1
        else:
            criterium_name = "absmag"
            info['rejects'].append(criterium_name)
        
        return info


    def inspect_lc(self, lc: LightCurve) -> None | dict[str, Any]:
        """
        Verify whether the transient lightcurve fulfill criteria for submission.

        """

        # apply cut on history: consider photophoints which are sharp enough
        pps = lc.get_photopoints(filters=self.lc_filters) # pps: photopoints from lightcurve w filters

        assert pps is not None
        info: dict[str, Any] = {'pass': 0, "rejects": []}
        

        # cut on number of detection
        # if len(pps) < self.min_ndet:
        #     self.logger.info(
        #         'Rejected', extra={'det': len(pps)}
        #     )
        #     return None
        # info["detections"] = len(pps)

        # pass on number of detections
        if len(pps) >= self.min_ndet:
            info['pass'] += 1
        elif len(pps) == 0:
            criterium_name = "no pps"
            info["rejects"].append(criterium_name)
            info["pass"] -= 10 # TODO maybe find a better way to punish this
            return info
        else:
            criterium_name = "pps"
            info['rejects'].append(criterium_name)
        info['detections'] = len(pps)

        # cut on age
        # jds = [pp["body"]["jd"] for pp in pps]
        # most_recent_detection, first_detection = max(jds), min(jds)
        # age = most_recent_detection - first_detection
        # if age > self.max_age or age < self.min_age:
        #     self.logger.info('Rejected', extra={'age': age})
        #     return None
        # info["age"] = age

        # pass on age
        jds = [pp["body"]["jd"] for pp in pps]
        most_recent_detection, first_detection = max(jds), min(jds)
        age = most_recent_detection - first_detection
        if age <= self.max_age and age >= self.min_age:
            info['pass'] += 1
        else:
            criterium_name = "age"
            info['rejects'].append(criterium_name)
        info["age"] = age


        # cut on number of detection after last SIGNIFICANT UL
        ulims = lc.get_upperlimits(
            filters={
                "attribute": "diffmaglim",
                "operator": ">=",
                "value": self.maglim_min,
            }
        ) # upper limits of lightcurve w filter: diffmaglim higher than maglim_min

        if ulims and len(ulims) > 0:
            last_ulim_jd = sorted([x["body"]["jd"] for x in ulims])[-1]
            pps_after_ndet = lc.get_photopoints(
                filters=self.lc_filters + [{"attribute": "jd", "operator": ">=", "value": last_ulim_jd}]
            )
            info["pass"] += 1
            # Check if there are enough positive detection after the last significant UL
            if (
                pps_after_ndet is not None and
                len(pps_after_ndet) < self.min_ndet_postul
            ):
                self.logger.info(
                    "not enough consecutive detections after last significant UL.",
                    extra={"NDet": len(pps), "lastUlimJD": last_ulim_jd},
                )

                criterium_name = "min_ndet_postul"
                info['rejects'].append(criterium_name)
                #return None
            else:
                info["pass"] += 1

            # Check that there is a recent ul
            if (most_recent_detection - last_ulim_jd) > self.maglim_maxago:
                self.logger.info(
                    "No recent UL.",
                    extra={
                        "lastDet": most_recent_detection,
                        "lastUlimJD": last_ulim_jd,
                    },
                )
                criterium_name = "maglim_maxago"
                info['rejects'].append(criterium_name)
                #return None
            else:
                info["pass"] += 1
            info["last_UL"] = most_recent_detection - last_ulim_jd
        else:
            self.logger.info("no UL")
            criterium_name = "no_UL"
            info['rejects'].append(criterium_name)
            #return None

        # cut on number of filters
        used_filters = set([pp["body"]["fid"] for pp in pps])
        if len(used_filters) < self.min_n_filters:
            self.logger.info(
                "Rejected", extra={'nbr_filt': len(used_filters)}
            )
            info["rejects"].append("min_n_filters")
            #return None
        else:
            info["pass"] += 1
        # cut on which filters used
        if used_filters.isdisjoint(self.det_filterids):
            self.logger.info(
                "Rejected (wrong filter det)", extra={'det_filters': used_filters}
            )
            criterium_name = "wrong_filter"
            info['rejects'].append(criterium_name)
            #return None
        info["pass"] += 1

        # cut on range of peak magnitude
        mags = [pp["body"]["magpsf"] for pp in pps]
        peak_mag = min(mags)
        if peak_mag > self.min_peak_mag or peak_mag < self.max_peak_mag:
            self.logger.info(
                "Rejected", extra={'peak_mag': peak_mag}
            )
            criterium_name = "peak_mag"
            info['rejects'].append(criterium_name)
            #return None
        else:
            info["pass"] += 1
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
                criterium_name = "unique latest_pps"
                info["rejects"].append(criterium_name)
                info["pass"] -= 5
                return info
                #raise ValueError("Have assumed a unique last photopoint")
            info["latest_mag"] = latest_pps[0]["body"]["magpsf"]

        # TODO: cut based on the mag rise per day (see submitRapid)

        # cut on galactic coordinates
        if pos := lc.get_pos(ret="mean", filters=self.lc_filters):
            ra, dec = pos
        else:
            criterium_name = "lc no points"
            info["rejects"].append(criterium_name)
            info["pass"] -= 5
            return info
            #raise ValueError("Light curve contains no points")
        coordinates = SkyCoord(ra, dec, unit="deg")
        b = coordinates.galactic.b.deg
        if abs(b) < self.min_gal_lat:
            self.logger.info(
                "Rejected (galactic plane)", extra={'gal_lat_b': b}
            )
            criterium_name = "min_gal_lat"
            info['rejects'].append(criterium_name)
            #return None
        else:
            info["pass"] += 1
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
            criterium_name = "close_to_sso"
            info['rejects'].append(criterium_name)
            #return None
        else:
            info["pass"] += 1

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
                criterium_name = "is_ps1_star"
                info['rejects'].append(criterium_name)
                info["pass"] -= 5
                #return None
            else:
                info["pass"] += 0 # TODO: think about whether optional checks should be rewarded more
        else:
            self.logger.info("No PS1 check as no data found.")

        # cut on median RB and DRB score
        rbs = [pp["body"]["rb"] for pp in pps]
        if np.median(rbs) < self.rb_minmed:
            self.logger.info(
                "Rejected (RB)",
                extra={"median_rb": np.median(rbs)},
            )
            criterium_name = "rb_minmed"
            info['rejects'].append(criterium_name)
            #return None
        elif (len(rbs) == 0) and self.rb_minmed > 0:
            self.logger.info("Rejected (No rb info)")
            criterium_name = "no rb info"
            info['rejects'].append(criterium_name)
            #return None
        else: 
            info["pass"] += 1
        info["rb"] = np.median(rbs)
        

        # drb might not exist
        drbs = [pp["body"]["drb"] for pp in pps if "drb" in pp["body"]]
        if len(drbs) > 0 and np.median(drbs) < self.drb_minmed:
            self.logger.info(
                "Rejected (dRB)",
                extra={"median_drd": np.median(drbs)},
            )
            criterium_name = "drb_minmed"
            info['rejects'].append(criterium_name)
            #return None
        elif (len(drbs) == 0) and self.drb_minmed > 0:
            self.logger.info("Rejected (No drb info)")
            criterium_name = "no drb info"
            info['rejects'].append(criterium_name)
            #return None
        else:
            info["pass"] += 1

        info["drb"] = np.median(drbs)

        # Transient passed pure LC criteria
        self.logger.info("Passed T2infantCatalogEval", extra=info)
        return info
    

    def inspect_catmatch(self, t2res: dict[str, Any]) -> None | dict[str, Any]:
        """
        Check wether any catalog has a match for the transit.
        """
        catalogKeys = t2res.keys()
        # print(catalogKeys)

        info = {'pass':0, 'rejects': []}

        for cat in catalogKeys:
            if (t2res[cat] is not None):
                info["pass"] -= 5
                info["rejects"].append(cat)
                info[cat] = t2res[cat]
                
            else:
                pass
                #info["pass"] -= 5

        #print(t2res.keys())
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

        kilonovaness: int = 0
        z_kilonovaness: int = 0
        lc_kilonovaness: int = 0
        possis_kilonovaness: int = 0
        cat_kilonovaness: int = 0
        info = {'possis':[]}
        rejects = []
        
        # Check t2 ouputs
        for t2_view in t2_views:
            self.logger.info('Parsing t2 results from {}'.format(t2_view.unit))
            t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res

            # Redshift
            if t2_view.unit == 'T2DigestRedshifts':
                zinfo = self.inspect_ampelz(t2_res)
                if len(zinfo["rejects"]) > 0:
                    rejects.extend(zinfo["rejects"])
                info.update(zinfo)
                kilonovaness += zinfo['pass']
                z_kilonovaness = zinfo["pass"]
            # Fit to kilonova model
            if t2_view.unit == 'T2RunPossis':
                pinfo = self.inspect_possis(t2_res)
                info['possis'].extend(pinfo)   # Could be multiple possis fits
                if len(pinfo["rejects"]) > 0:
                    rejects.extend(pinfo["rejects"])
                kilonovaness += pinfo['pass']
                possis_kilonovaness = pinfo["pass"]
                if 'possis_abspeak' in info.keys():
                    if info['possis_chisq']>zinfo['possis_chisq']:
                        info.update( pinfo )
                else:
                    info.update( pinfo )

            if t2_view.unit == "T2CatalogMatch":
                #print(info("Catalog match {}".format(t2_res.keys)))
                cinfo = self.inspect_catmatch(t2_res)
                if len(cinfo["rejects"]) > 0:
                    rejects.extend(cinfo["rejects"])
                info.update(cinfo)
                kilonovaness += cinfo['pass']
                cat_kilonovaness = cinfo["pass"]
                #print("T2CatalogMatch in kilonovaeval")

            # Propagate map info
            if t2_view.unit == 'T2PropagateStockInfo':
                info.update( t2_res )   # Could there be multiple maps associated? E.g. after updates? TODO



        # Check whether the lightcurve passes selection criteria
        # TODO: add kilonovaness criteria
        lc_info = self.inspect_lc(light_curve)
        if lc_info:
           info.update(lc_info)
           kilonovaness += lc_info['pass']
           lc_kilonovaness = lc_info["pass"]
           if len(lc_info["rejects"]) > 0:
            rejects.extend(lc_info["rejects"])


        # iii. Check absolute magnitude - again (but directly from lightcurve)
        if (z := info.get('ampel_z')) and (obsmag := info.get("peak_mag")):
            sndist = Distance(z=z, cosmology=Planck15)
            info["absmag"] = obsmag - sndist.distmod.value
            if (self.min_absmag < info['absmag'] < self.max_absmag):
                kilonovaness += 1
                lc_kilonovaness += 1
            else:
                criterium_name = "absmag from lc"
                rejects.append(criterium_name)

        # Categorize
        if (kilonovaness < 1):
            rank_decimal = 0
        else:
            rank_decimal = kilonovaness/(kilonovaness + len(rejects))
        
        if rank_decimal > .9: #TODO arbitrary rn
            info['is_gold'] = True
        info['kilonovaness'] = kilonovaness
        info["kilonovaness_dec"] = rank_decimal
        info["z_kilonovaness"] = z_kilonovaness
        info["lc_kilonovaness"] = lc_kilonovaness
        info["possis_kilonovaness"] = possis_kilonovaness
        info["cat_kilonovaness"] = cat_kilonovaness
        info["rejects"] = rejects

        return info
