from typing import Any

import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.ingest.AbsArchiveAugmenter import AbsArchiveAugmenter
from ampel.types import Tag
from ampel.view.ReadOnlyDict import ReadOnlyDict
from ampel.ztf.base.ArchiveUnit import ArchiveUnit
from ampel.ztf.util.ZTFIdMapper import to_ampel_id


def mean_position(ra: list[float], dec: list[float]) -> tuple[float, float]:
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
    return np.degrees(mean_ra), np.degrees(mean_dec)


class ZiArchiveAugmenter(AbsArchiveAugmenter, ArchiveUnit):
    """
    Add datapoints from archived ZTF alerts.
    """

    # the radius in which to look for ZTF counterparts
    radius_arcsec: float = 1

    # tag to add to the alert
    tag: Tag | None | list[Tag] = None

    @staticmethod
    def select_closest_ztf_source(
        alert_dicts: list[dict[str, Any]], ra_ref: float, dec_ref: float
    ) -> tuple[list[dict[str, Any]], float]:
        ztf_names = np.unique([ad["objectId"] for ad in alert_dicts])
        sorted_alerts_dict = {
            zn: [ad for ad in alert_dicts if ad["objectId"] == zn] for zn in ztf_names
        }
        min_distance = np.inf
        closest_ztf_name = None
        for zn, ads in sorted_alerts_dict.items():
            mp = mean_position(
                [float(ad["candidate"]["ra"]) for ad in ads],
                [float(ad["candidate"]["dec"]) for ad in ads],
            )
            distance = np.sqrt(
                (np.radians(mp[0]) - np.radians(ra_ref)) ** 2
                * np.cos(np.radians(dec_ref)) ** 2
                + (np.radians(mp[1]) - np.radians(dec_ref)) ** 2
            )
            if distance < min_distance:
                min_distance = distance
                closest_ztf_name = zn

        if closest_ztf_name is None:
            # We checked earlier already that there are some alerts so this should not happen
            raise ValueError("No ZTF source found!")

        # selecting the alerts for the closest ZTF source
        return sorted(
            sorted_alerts_dict[closest_ztf_name], key=lambda x: x["candidate"]["jd"]
        ), min_distance

    def get_alerts(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert | None:
        detections = [
            dp for dp in dps if "ra" in dp["body"]
        ]  # TODO: refine detection definition
        ra, dec = mean_position(
            [dp["body"]["ra"] for dp in detections],
            [dp["body"]["dec"] for dp in detections],
        )
        url = (
            f"/alerts/cone_search?"
            f"ra={ra}&dec={dec}&radius={self.radius_arcsec / 3600}&"
            f"jd_start={jd_center - time_pre}&jd_end={jd_center + time_post}"
        )
        queried_alerts = self.session.get(url)
        queried_detections = [
            a for a in queried_alerts if "ra" in a["candidate"]
        ]  # TODO: refine detection definition
        if len(queried_detections) == 0:
            return None
        ztf_alerts_from_closest_object, distance_to_closest_ztf_source = (
            self.select_closest_ztf_source(queried_detections, ra, dec)
        )

        latest_alert_from_closest_source = ztf_alerts_from_closest_object[-1]

        return AmpelAlert(
            id=latest_alert_from_closest_source["candid"],
            stock=to_ampel_id(latest_alert_from_closest_source["objectId"]),
            datapoints=tuple(d["candidate"] for d in ztf_alerts_from_closest_object),
            extra=ReadOnlyDict({"name": latest_alert_from_closest_source["objectId"]}),
            tag=self.tag,
        )
