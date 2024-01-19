#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                Ampel-contrib-HU/ampel/contrib/hu/t0/TransientInClusterFilter.py
# License:             BSD-3-Clause
# Author:              m. giomi <matteo.giomi@desy.de>
# Date:                28.08.2018
# Last Modified Date:  24.11.2021
# Last Modified By:    Jakob van Santen <jakob.van.santen@desy.de>

import numpy as np
from functools import partial
from ampel.ztf.t0.DecentFilter import DecentFilter
from ampel.protocol.AmpelAlertProtocol import AmpelAlertProtocol


class TransientInClusterFilter(DecentFilter):
    """
    Filter derived from the DecentFilter, in addition selecting candidates
    with position compatible with that of nearby galaxy clusters.
    """

    big_search_radius_arcmin: float  # conservative search radius around cluster position. Max in RASSEBCS is 16.a arcmin
    cluserter_rad_multiplier: float  # if we want to enlarge the search region around each cluster.

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # feedback
        for k in self.__annotations__:
            self.logger.info(f"Using {k}={getattr(self, k)}")

        self.rassebcs_query = partial(
            self.cone_search_all,
            catalogs=[
                {
                    "name": "RASSEBCS",
                    "use": "extcats",
                    "rs_arcsec": self.big_search_radius_arcmin
                    * 60
                    * self.cluserter_rad_multiplier,
                    "keys_to_append": ["ANGULAR_RADIUS"],
                }
            ],
        )

    def process(self, alert: AmpelAlertProtocol):
        """
        run the filter on the alert. First we run the decent filter, then we match
        with the cluster catalog.
        """

        # if the candidate has passed the decent filter, check if it is compatible
        # with the position of some nearby galaxy cluster
        latest = alert.datapoints[0]
        alert_ra = latest["ra"]
        alert_dec = latest["dec"]

        # A) find position of all the nearby clusters. If none is found, reject alert.
        nearby_clusters = self.rassebcs_query(alert_ra, alert_dec)[0]

        if nearby_clusters is None:
            self.logger.debug(
                f"rejected: no cluster from RASSEBCS within {self.big_search_radius_arcmin*60:.2f} arcsec from alert position"
            )
            return None

        # B) for all the nearby clusters, compute their distances to the alert
        # (extcats works in arcsec but ANGULAR_RADIUS is in arcmin)
        clust_dists_2_alert = (
            np.array([d["dist_arcsec"] for d in nearby_clusters]) / 60.0
        )
        clust_radius = np.array([d["body"]["ANGULAR_RADIUS"] for d in nearby_clusters])

        # C) if, for any of the nearby clusters, the distance to the alert is smaller
        # than the cluster size, count this as a match
        if not any(clust_radius > clust_dists_2_alert):
            self.logger.debug(
                "rejected: distance to alert is larger than the cluster size for all the nearby clusters"
            )

            for ii in range(len(nearby_clusters)):
                self.logger.debug(
                    f"Angular radius: {clust_radius[ii]:.2f}, "
                    f"distance to alert: {clust_dists_2_alert[ii]:.2f} (both arcmin)"
                )

            return None

        for ii in range(len(nearby_clusters)):
            self.logger.debug(
                f"Angular radius: {clust_radius[ii]:.2f}, "
                f"distance to alert: {clust_dists_2_alert[ii]:.2f} (both arcmin)"
            )

        # now run the decent filter (faster to do it afterwards ;-)
        return super().process(alert)
