from ampel.abstract.AbsT0Muxer import AbsT0Muxer
from ampel.abstract.AbsT0Unit import AbsT0Unit
from ampel.alert.AmpelAlert import AmpelAlert
from ampel.base.AmpelABC import AmpelABC
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.core.ContextUnit import ContextUnit
from ampel.log import AmpelLogger
from ampel.model.UnitModel import UnitModel
from ampel.types import StockId


class AbsAdder(AmpelABC, ContextUnit, abstract=True):
    """
    Add alerts from another alert source. This explicitly assumes
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

    # process augmented datapoints
    mux: UnitModel | str
    augmenting_shaper: UnitModel | str

    logger: AmpelLogger

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        self._adding_shaper = self.context.loader.new_logical_unit(
            model=UnitModel(unit=self.augmenting_shaper)
            if isinstance(self.augmenting_shaper, str)
            else self.augmenting_shaper,
            logger=self.logger,
            sub_type=AbsT0Unit,
        )
        self._muxer = self.context.loader.new_context_unit(
            model=UnitModel(unit=self.mux) if isinstance(self.mux, str) else self.mux,
            logger=self.logger,
            context=self.context,
            sub_type=AbsT0Muxer,
        )

        self._t0_col = self.context.db.get_collection("t0", "w")

    @abstractmethod
    def add(
        self, dps: list[DataPoint], jd_center: float, time_pre: float, time_post: float
    ) -> list[AmpelAlert] | None: ...

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

        # Obtain augment alert
        add_alerts = self.add(dps, alert_jd, self.history_days, self.future_days)

        if not add_alerts:
            # nothing found in archive
            return [], []

        # process datapoints for individual alerts independently
        add_insert = []
        add_combine = []
        for a in add_alerts:
            add_dps = self._adding_shaper.process(a.datapoints, a.stock)

            # the muxer should check the database for already inserted datapoints
            mux_res = self._muxer.process(add_dps, a.stock)

            add_insert.extend(mux_res[0])
            add_combine.extend(mux_res[1])

        return add_insert, add_combine
