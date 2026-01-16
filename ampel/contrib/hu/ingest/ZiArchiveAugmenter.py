from typing import Any

import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.ingest.AbsArchiveAugmenter import AbsArchiveAugmenter
from ampel.ztf.alert.ZiAlertSupplier import ZiAlertSupplier
from ampel.ztf.base.ArchiveUnit import ArchiveUnit


class ZiArchiveAugmenter(AbsArchiveAugmenter, ArchiveUnit):
    """
    Add datapoints from archived ZTF alerts.
    """

    # the radius in which to look for ZTF counterparts
    radius_arcsec: float = 1

    def select_closest_ztf_source(
        self, alert_dicts: dict[str, Any]
    ) -> dict[str, Any]: ...

    def get_alerts(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert:
        ra = float(np.median([dp["body"]["ra"] for dp in dps]))
        dec = float(np.median([dp["body"]["dec"] for dp in dps]))
        url = f"/alerts/cone_search?ra={ra}&dec={dec}&radius={self.radius_arcsec}&jd_start={jd_center - time_pre}&jd_end={jd_center + time_post}"
        response = self.session.get(url)
        response.raise_for_status()
        closest_ztf_alerts = self.select_closest_ztf_source(response.json())
        return ZiAlertSupplier.shape_alert_dict(closest_ztf_alerts)
