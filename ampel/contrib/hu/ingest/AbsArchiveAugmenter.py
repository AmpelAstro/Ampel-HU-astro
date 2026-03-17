from ampel.abstract.AbsT0Muxer import AbsT0Muxer
from ampel.abstract.AbsT0Unit import AbsT0Unit
from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.model.UnitModel import UnitModel
from ampel.types import StockId, Tag


class AbsArchiveAugmenter(AbsT0Muxer, abstract=True):
    """
    Augment alerts with data from another alert source. This explicitly assumes
    that the source of the primary stream is deeper than the source of the
    augmenting data. This means that all objects detected in the primary stream
    should also be detected by the augmenting stream.
    TODO: Ideally the matches would be recorded in the database. Maybe a resource should
        be defined where a suitable mongo db can be specified?
    """

    #: Number of days of history to add, relative to the earliest point in the
    #: t0 collection
    history_days: float = 0
    #: Number of days of looking ahead based on the JD contained in alert
    #: Only effect if original alert obtained through time-dependent search
    #:Warning: could interfer if further alserts added ex through ZiMongoMuxer
    future_days: float = 0

    augmenting_shaper: UnitModel | str

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

        self._augmenting_shaper = self.context.loader.new_logical_unit(
            model=UnitModel(unit=self.augmenting_shaper)
            if isinstance(self.augmenting_shaper, str)
            else self.augmenting_shaper,
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
        # TODO: check if these ztf alerts are already in the database.
        #   In that case we'd have to check if they are associated to any other alert
        #   from the primary survey and if that association is better or worse.

        if not archive_alert:
            # nothing found in archive
            return dps, dps

        archive_dps = self._augmenting_shaper.process(
            archive_alert.datapoints, stock_id
        )

        # Create combined state of alert and archive
        # Add all dps because the archive has data from a different instrument
        extended_dps = sorted(dps + archive_dps, key=lambda d: d["body"]["jd"])

        return extended_dps, extended_dps
