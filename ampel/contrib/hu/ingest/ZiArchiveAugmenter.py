from typing import Any

import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.ingest.AbsArchiveAugmenter import AbsArchiveAugmenter
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.base.ArchiveUnit import ArchiveUnit


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

    @staticmethod
    def select_closest_ztf_source(
        alert_dicts: list[dict[str, Any]], ra_ref: float, dec_ref: float
    ) -> dict[str, Any]:
        ztf_names = np.unique([ad["objectId"] for ad in alert_dicts])
        sorted_alerts_dict = {
            zn: [ad["candidate"] for ad in alert_dicts if ad["objectId"] == zn]
            for zn in ztf_names
        }
        mp = np.inf
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
            if distance < mp:
                mp = distance
                closest_ztf_name = zn

        if closest_ztf_name is None:
            # We checked earlier already that there are some alerts so this should not happen
            raise ValueError("No ZTF source found!")

        # selecting the alerts for the closest ZTF source
        selected_alerts = sorted(
            sorted_alerts_dict[closest_ztf_name], key=lambda x: x["candidate"]["jd"]
        )

        # mimicking the structure of a photopoints query result, so we can use the ZiAlertSupplier later
        pps_dict = selected_alerts[-1]
        pps_dict["prv_candidates"] = [ad for ad in selected_alerts[:-1]]
        return pps_dict

    def get_alerts(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert | None:
        ra, dec = mean_position(
            [dp["body"]["ra"] for dp in dps], [dp["body"]["dec"] for dp in dps]
        )
        url = f"/alerts/cone_search?ra={ra}&dec={dec}&radius={self.radius_arcsec / 3600}&jd_start={jd_center - time_pre}&jd_end={jd_center + time_post}"
        response = self.session.get(url)
        response.raise_for_status()
        queried_alerts = response.json()
        if len(queried_alerts) == 0:
            return None
        closest_ztf_alerts = self.select_closest_ztf_source(response.json(), ra, dec)
        return ZiAlertSupplier.shape_alert_dict(closest_ztf_alerts)
