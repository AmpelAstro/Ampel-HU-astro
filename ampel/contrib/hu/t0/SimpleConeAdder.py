from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t0.AbsAdder import AbsAdder
from ampel.contrib.hu.util.meanpos import mean_position
from ampel.types import Tag


class SimpleConeAdder(AbsAdder, abstract=True):
    """
    Simply get all alerts from secondary alert source within a radius
    """

    radius_arcsec: float

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
