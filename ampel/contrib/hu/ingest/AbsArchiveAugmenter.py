from ampel.abstract.AbsT0Muxer import AbsT0Muxer
from ampel.abstract.AbsT0Unit import AbsT0Unit
from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.model.UnitModel import UnitModel
from ampel.types import StockId, Tag


class AbsArchiveAugmenter(AbsT0Muxer):
    """
    Add datapoints from archived alerts.
    """

    #: Number of days of history to add, relative to the earliest point in the
    #: t0 collection
    history_days: float = 0
    #: Number of days of looking ahead based on the JD contained in alert
    #: Only effect if original alert obtained through time-dependent search
    #:Warning: could interfer if further alserts added ex through ZiMongoMuxer
    future_days: float = 0

    shaper: UnitModel | str

    tag: Tag

    # Standard projection used when checking DB for existing PPS/ULS
    projection: dict[str, int] = {
        "_id": 1,
        "tag": 1,
        "excl": 1,
        "body.jd": 1,
        "body.fid": 1,
        "body.rcid": 1,
        "body.magpsf": 1,
    }

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        self._shaper = self.context.loader.new_logical_unit(
            model=UnitModel(unit=self.shaper)
            if isinstance(self.shaper, str)
            else self.shaper,
            logger=self.logger,
            sub_type=AbsT0Unit,
        )

        self._t0_col = self.context.db.get_collection("t0", "w")

    @abstractmethod
    def get_alerts(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> AmpelAlert | None: ...

    def process(
        self, dps: list[DataPoint], stock_id: None | StockId = None
    ) -> tuple[None | list[DataPoint], None | list[DataPoint]]:
        """
        :param dps: datapoints from alert
        :param stock_id: stock id from alert
        Attempt to determine which pps/uls should be inserted into the t0 collection,
        and which one should be marked as superseded.
        """

        if not stock_id:
            # no new points to add; use input points for combination
            return dps, dps

        # Alert jd, assumed to be latest dp
        alert_jd = max(
            dp["body"]["jd"] for dp in dps if dp["id"] > 0 and "ZTF" in dp["tag"]
        )

        # Obtain archive alert
        archive_alert = self.get_alerts(
            dps, alert_jd, self.history_days, self.future_days
        )
        if not archive_alert:
            # nothing found in archive
            return dps, dps
        archive_dps = self._shaper.process(archive_alert.datapoints, stock_id)

        # Create combined state of alert and archive
        # Add all dps because the archive has data from a different instrument
        jds_alert = [dp["body"]["jd"] for dp in dps]
        extended_dps = sorted(
            dps + [dp for dp in archive_dps if dp["body"]["jd"] not in jds_alert],
            key=lambda d: d["body"]["jd"],
        )

        # TODO:
        #  check the DB content, both to extend
        #  the state and avoid duplicate inserts. Would then also need to
        #  check for superseeded dps.

        return extended_dps, extended_dps
