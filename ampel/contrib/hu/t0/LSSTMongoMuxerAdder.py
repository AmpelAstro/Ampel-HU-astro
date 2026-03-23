from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.t0.AbsAdder import AbsAdder
from ampel.lsst.ingest.LSSTMongoMuxer import LSSTMongoMuxer
from ampel.model.UnitModel import UnitModel
from ampel.types import StockId


class LSSTMongoMuxerAdder(LSSTMongoMuxer):
    add: UnitModel

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._adder = self.context.loader.new_context_unit(
            model=self.add,
            context=self.context,
            logger=self.logger,
            sub_type=AbsAdder,
        )

    def process(
        self, dps_al: list[DataPoint], stock_id: None | StockId = None
    ) -> tuple[None | list[DataPoint], None | list[DataPoint]]:
        dps_insert, dps_combine = super().process(dps_al, stock_id)
        aug_insert, aug_combine = self._adder.process(dps_combine, stock_id)
        return dps_insert + aug_insert, dps_combine + aug_combine
