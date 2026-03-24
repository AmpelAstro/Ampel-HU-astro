import numpy as np

from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t0.AbsAdder import AbsAdder
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


class SimpleConeAdder(AbsAdder, abstract=True):
    """
    Simply get all alerts from secondary alert source within a radius
    """

    radius_arcsec: float

    # tag to add to the alert
    tag: Tag | None | list[Tag] = None

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

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

    def add(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert | None:
        # only use detection of the primary stream
        detections = [
            dp for dp in dps if "ra" in dp["body"]
        ]  # TODO: refine detection definition

        # calculate mean position and standard deviation of alerts
        ra, dec, _, _ = mean_position(
            [dp["body"]["ra"] for dp in detections],
            [dp["body"]["dec"] for dp in detections],
        )

        # get alerts within radius from secondary stream
        cone_search_res = self.cone_search(
            ra, dec, self.radius_arcsec, jd_center - time_pre, jd_center + time_post
        )
        if len(cone_search_res) == 0:
            return None

        return self.shape_alert_dict(cone_search_res)
