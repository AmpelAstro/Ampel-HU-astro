import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.contrib.hu.t0.SimpleConeAdder import SimpleConeAdder
from ampel.view.ReadOnlyDict import ReadOnlyDict
from ampel.ztf.base.ArchiveUnit import ArchiveUnit
from ampel.ztf.util.ZTFIdMapper import to_ampel_id


class ZiArchiveAdder(SimpleConeAdder, ArchiveUnit):
    """
    Add datapoints from archived ZTF alerts.
    """

    def cone_search(
        self,
        ra: float,
        dec: float,
        radius_arcsec: float,
        jd_start: float,
        jd_end: float,
    ) -> list[dict]:
        url = (
            f"/alerts/cone_search?"
            f"ra={ra}&dec={dec}&radius={radius_arcsec / 3600}&"
            f"jd_start={jd_start}&jd_end={jd_end}"
        )
        queried_alerts = self.session.get(url)
        queried_detections = [
            a for a in queried_alerts if "ra" in a["candidate"]
        ]  # TODO: refine detection definition
        return queried_detections

    def shape_alert_dict(self, dp: list[dict]) -> list[AmpelAlert]:
        objids = np.unique([d["objectId"] for d in dp])

        alerts = []
        for objid in objids:
            idp = sorted(
                filter(
                    lambda d: d["objectId"] == objid if "objectId" in d else False, dp
                ),
                key=lambda x: x["candidate"]["jd"],
                reverse=True,
            )
            latest_dp = idp[0]
            alert = AmpelAlert(
                id=latest_dp["candid"],
                stock=to_ampel_id(latest_dp["objectId"]),
                datapoints=tuple(d["candidate"] for d in idp),
                extra=ReadOnlyDict({"name": latest_dp["objectId"]}),
            )
            alerts.append(alert)

        return alerts
