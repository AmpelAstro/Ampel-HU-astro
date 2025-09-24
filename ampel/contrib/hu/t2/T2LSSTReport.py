from collections.abc import Sequence

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.model.LSSTReport import (
    Classification,
    LSSTReport,
    ModelClassification,
    Object,
    PhotometricPoint,
)
from ampel.contrib.hu.model.ParsnipRiseDeclineResult import (
    ParsnipResult,
    ParsnipRiseDeclineResult,
)
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.UnitModel import UnitModel
from ampel.struct.UnitResult import UnitResult
from ampel.view.T2DocView import T2DocView


class T2LSSTReport(AbsTiedStateT2Unit, AbsTabulatedT2Unit):
    """
    Create an LSST report from T2RunParsnipRiseDecline output and additional information
    """

    result_adapter: UnitModel | None = None

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UnitResult:
        for view in t2_views:
            if view.unit == "T2RunParsnipRiseDecline" and (
                payload := view.get_payload(code=DocumentCode.OK)
            ):
                parsnip_result = ParsnipRiseDeclineResult(**payload)
                break
        else:
            raise ValueError("T2RunParsnipRiseDecline result not found in T2 views")
        stock = compound["stock"]
        if not isinstance(stock, int):
            raise TypeError(f"Stock ID should be an integer, got {stock}")

        photometry = [
            PhotometricPoint(
                time=time,
                flux=flux,
                fluxerr=fluxerr,
                band=band,
                zp=zp,
                zpsys=zpsys,
            )
            for time, flux, fluxerr, band, zp, zpsys in self.get_flux_table(
                datapoints
            ).iterrows("time", "flux", "fluxerr", "band", "zp", "zpsys")
        ]
        # Fill object record from latest LSST_OBJ datapoint (diaObject)
        for dp in sorted(
            datapoints, key=lambda x: x["meta"][-1].get("ts", 0), reverse=True
        ):
            if "LSST_OBJ" in dp.get("tag", {}):
                obj = Object(
                    id=stock,
                    ra=float(dp["body"]["ra"]),
                    ra_err=float(dp["body"]["raErr"]),
                    dec=float(dp["body"]["dec"]),
                    dec_err=float(dp["body"]["decErr"]),
                    ra_dec_cov=float(dp["body"]["ra_dec_Cov"]),
                    # FIXME: add redshift if available
                )
                break
        else:
            raise ValueError("No LSST_OBJ found in datapoints")
        classifications = [
            Classification(
                name="ParsnipRiseDecline",
                version="1",
                models=[
                    ModelClassification(
                        model=model.model, probabilities=model.classification
                    )
                    for model in parsnip_result.classifications[0].parsnip
                    if isinstance(model, ParsnipResult)
                ],
            )
        ]
        report = LSSTReport(
            object=obj, photometry=photometry, classification=classifications
        )
        return UnitResult(body=report.model_dump(), adapter=self.result_adapter)
