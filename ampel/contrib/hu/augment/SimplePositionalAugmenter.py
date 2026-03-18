from typing import Any

import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.augment.AbsAugmenter import AbsAugmenter
from ampel.types import Tag

ARCSEC_IN_RAD = np.pi / 180 / 3600
SQDEG_IN_SR = (np.pi / 180) ** 2


def mean_position(
    ra: list[float], dec: list[float]
) -> tuple[float, float, float, float]:
    ra = np.radians(ra)
    dec = np.radians(dec)
    ref_ra = ra[0]
    ref_dec = dec[0]
    dx = (ra - ref_ra) * np.cos(ref_dec)
    dy = dec - ref_dec
    mean_dx = np.mean(dx)
    mean_dy = np.mean(dy)
    mean_ra = ref_ra + mean_dx / np.cos(ref_dec)
    mean_dec = ref_dec + mean_dy
    dra = ref_ra + dx / np.cos(ref_dec)
    ddec = ref_dec + dy
    std_ra = np.std(dra)
    std_dec = np.std(ddec)
    return (
        np.degrees(mean_ra),
        np.degrees(mean_dec),
        np.degrees(std_ra),
        np.degrees(std_dec),
    )


class SimplePositionalAugmenter(AbsAugmenter, abstract=True):
    """
    Augment based on the position of alerts.
    Based on https://www.overleaf.com/read/hpbjsjrrxpym#b8d122
    """

    # the minimum posterior probability for any match
    min_posterior: float = 0.9

    # The density of sources in the primary alert stream in
    # 1/deg^2 used as the prior.
    # In the future this could be replaced with a
    # general PriorModel that can be evaluated per sky position.
    nu1: float

    # the resolutions of the two alert streams in arcseconds
    sigma1: float
    sigma2: float

    # tag to add to the alert
    tag: Tag | None | list[Tag] = None

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        sigma_sq_rad = (self.sigma1**2 + self.sigma2**2) * ARCSEC_IN_RAD**2
        self._rho1 = 4 * np.pi * self.nu1 / SQDEG_IN_SR
        self._radius_arcsec = (
            np.sqrt(
                np.log((1 / self.min_posterior - 1) / (self._rho1 * sigma_sq_rad) * 2)
                * 2
                * sigma_sq_rad
            )
            / ARCSEC_IN_RAD
        )

    @abstractmethod
    def name_from_alert(self, dp: dict) -> str: ...
    @abstractmethod
    def ra_from_alert(self, dp: dict) -> float: ...
    @abstractmethod
    def dec_from_alert(self, dp: dict) -> float: ...
    @abstractmethod
    def jd_from_alert(self, dp: dict) -> float: ...

    @abstractmethod
    def cone_search(
        self,
        ra: float,
        dec: float,
        radius_arcsec: float,
        jd_start: float,
        jd_end: float,
    ) -> list[dict]: ...

    @abstractmethod
    def shape_alert_dict(self, dp: dict) -> AmpelAlert: ...

    def select_closest_source(
        self, alert_dicts: list[dict[str, Any]], ra_ref: float, dec_ref: float
    ) -> tuple[list[dict[str, Any]], float, tuple[float, float, float, float]]:
        # sort the alerts per source name
        source_names = np.unique([self.name_from_alert(ad) for ad in alert_dicts])
        sorted_alerts_dict = {
            zn: [ad for ad in alert_dicts if self.name_from_alert(ad) == zn]
            for zn in source_names
        }

        # loop over all sources and find closest one
        min_distance = np.inf
        closest_source_name = None
        associated_position = None
        for zn, ads in sorted_alerts_dict.items():
            mp = mean_position(
                [self.ra_from_alert(ad) for ad in ads],
                [self.dec_from_alert(ad) for ad in ads],
            )
            distance = np.sqrt(
                (np.radians(mp[0]) - np.radians(ra_ref)) ** 2
                * np.cos(np.radians(dec_ref)) ** 2
                + (np.radians(mp[1]) - np.radians(dec_ref)) ** 2
            )
            if distance < min_distance:
                min_distance = distance
                closest_source_name = zn
                associated_position = mp

        if closest_source_name is None:
            # We checked earlier already that there are some alerts so this should not happen
            raise ValueError("No source found!")

        # selecting the alerts for the closest source
        return (
            sorted(sorted_alerts_dict[closest_source_name], key=self.jd_from_alert),
            min_distance,
            associated_position,
        )

    def augment(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert | None:
        # only use detection of the primary stream
        detections = [
            dp for dp in dps if "ra" in dp["body"]
        ]  # TODO: refine detection definition

        # calculate mean position and standard deviation of alerts
        ra, dec, dra, ddec = mean_position(
            [dp["body"]["ra"] for dp in detections],
            [dp["body"]["dec"] for dp in detections],
        )

        # get alerts within radius from secondary stream
        cone_search_res = self.cone_search(
            ra, dec, self._radius_arcsec, jd_center - time_pre, jd_center + time_post
        )
        if len(cone_search_res) == 0:
            return None

        # select closest source, calculate distance, position and standard deviation
        ztf_alerts_from_closest_object, distance_to_closest_ztf_source, pos = (
            self.select_closest_source(cone_search_res, ra, dec)
        )

        # calculate posterior association probability
        sigma1_sq = dra**2 + ddec**2
        sigma2_sq = pos[2] ** 2 + pos[3] ** 2
        sigma_sq_rad = (sigma1_sq + sigma2_sq) * SQDEG_IN_SR
        posterior = 1 / (
            self._rho1
            * sigma_sq_rad
            / 2
            * np.exp(distance_to_closest_ztf_source**2 / (2 * sigma_sq_rad))
            + 1
        )
        # TODO: save association probability somewhere (database specified in resources?)

        return self.shape_alert_dict(ztf_alerts_from_closest_object)
