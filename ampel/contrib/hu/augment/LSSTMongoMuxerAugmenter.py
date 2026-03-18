from ampel.content.DataPoint import DataPoint
from ampel.contrib.hu.augment.AbsAugmenter import AbsAugmenter
from ampel.lsst.ingest.LSSTMongoMuxer import LSSTMongoMuxer
from ampel.model.UnitModel import UnitModel
from ampel.types import StockId


class LSSTMongoMuxerAugmenter(LSSTMongoMuxer):
    augment: UnitModel

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self._augment = self.context.loader.new_context_unit(
            model=self.augment,
            context=self.context,
            logger=self.logger,
            sub_type=AbsAugmenter,
        )

    def process(
        self, dps_al: list[DataPoint], stock_id: None | StockId = None
    ) -> tuple[None | list[DataPoint], None | list[DataPoint]]:
        dps_insert, dps_combine = super().process(dps_al, stock_id)
        aug_insert, aug_combine = self._augment.process(dps_combine, stock_id)
        return dps_insert + aug_insert, dps_combine + aug_combine
