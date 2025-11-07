from collections.abc import Mapping, Sequence
from typing import Any

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.model.LSSTReport import (
    Classification,
    Feature,
    Host,
    LSSTReport,
    Object,
    PhotometricPoint,
)
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.UnitModel import UnitModel
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.struct.UnitResult import UnitResult
from ampel.view.T2DocView import T2DocView


class T2MissingDependency(RuntimeError):
    pass


class T2LSSTReport(AbsTiedStateT2Unit, AbsTabulatedT2Unit):
    """
    Create an LSST report for subsequent distribution. Base class connects current photometric information.

    Subclasses can allow filter and information from T2s to be added.

    """

    tabulator: Sequence[UnitModel] = [
        UnitModel(unit="LSSTT2Tabulator", config={"zp": 27.5})
    ]

    result_adapter: UnitModel | None = None

    def process_t2s(
        self, report_views: dict[str, Mapping[str, Any]]
    ) -> Sequence[Classification | Host | Feature] | None:
        """
        Process T2 views to extract information to be propagated.
        This is a placeholder for any processing needed on the T2 views.

        Return None if submission criteria not met, in which case no report is propagated.

        Can generate one of three kinds of information:
        - classification (Model name and probabilities, info)
        - host (name, redshift, source, info)
        - features (name, dict)

        """
        # Placeholder implementation
        return None

    def _get_payloads(
        self, t2_views: Sequence[T2DocView]
    ) -> dict[str, Mapping[str, Any]]:
        units = {dep.unit for dep in self.t2_dependency}
        payloads = {}
        for view in t2_views:
            assert isinstance(view.unit, str)
            if view.unit in units and (
                payload := view.get_payload(code=DocumentCode.OK)
            ):
                units.remove(view.unit)
                payloads[view.unit] = payload
            if not units:
                break
        if self.t2_dependency and units:
            raise T2MissingDependency(f"Missing T2 dependencies: {units}")
        return payloads

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UnitResult:
        # Inspect and collect T2 data if requested
        try:
            unit_reports = self.process_t2s(self._get_payloads(t2_views)) or []
        except T2MissingDependency:
            return UnitResult(code=DocumentCode.T2_MISSING_DEPENDENCY)
        if self.t2_dependency and not unit_reports:
            return UnitResult(
                code=DocumentCode.OK,
                journal=JournalAttributes(extra={"skipped": True}),
            )

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

        report = LSSTReport(
            object=obj,
            state=compound["link"],
            photometry=photometry,
            classification=[v for v in unit_reports if isinstance(v, Classification)],
            host=[v for v in unit_reports if isinstance(v, Host)],
            features=[v for v in unit_reports if isinstance(v, Feature)],
        )
        return UnitResult(body=report.model_dump(), adapter=self.result_adapter)
