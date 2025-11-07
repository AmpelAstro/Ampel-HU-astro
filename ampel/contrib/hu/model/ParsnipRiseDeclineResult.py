from collections.abc import Mapping
from functools import cached_property
from typing import Any

from pydantic import BaseModel, computed_field
from scipy.stats import chi2


class FitDataInfo(BaseModel):
    z: list[float]
    z_source: str
    z_weights: Any | None
    jdstart: float
    jdend: float


Classification = dict[str, float]
RiseDeclineFeatures = dict[str, Any]


class ParsnipPrediction(BaseModel):
    ra: float | None
    dec: float | None
    type: str
    redshift: float
    parsnip_reference_time: float
    parsnip_scale: float
    reference_time: float
    reference_time_error: float
    color: float
    color_error: float
    amplitude: float
    amplitude_error: float
    s1: float
    s1_error: float
    s2: float
    s2_error: float
    s3: float
    s3_error: float
    total_s2n: float
    count: int
    count_s2n_3: int
    count_s2n_5: int
    count_s2n_3_pre: int
    count_s2n_3_rise: int
    count_s2n_3_fall: int
    count_s2n_3_post: int
    model_chisq: float
    model_dof: int
    luminosity: float
    luminosity_error: float

    @computed_field  # type: ignore[prop-decorator]
    @cached_property
    def chi2pdf(self) -> float:
        return float(
            chi2.pdf(self.model_chisq, self.model_dof) if self.model_dof > 0 else 0.0
        )


class PredictionEntry(BaseModel):
    z: float
    prediction: ParsnipPrediction


class ClassificationEntry(BaseModel):
    z: float
    classification: Classification


class ParsnipResult(BaseModel):
    model: str
    predictions: list[PredictionEntry]
    classifications: list[ClassificationEntry]
    marginal_lc_classifications: Classification
    z_at_minchi: float
    prediction: ParsnipPrediction
    classification: Classification
    plot: Mapping[str, Any] | None = None


class ParsnipFailure(BaseModel):
    model: str
    failed: str


class Classifications(BaseModel):
    name: str
    version: str
    parsnip: list[ParsnipResult | ParsnipFailure]
    xgbbinary: dict[str, Any]
    xgbmulti: dict[str, Any]
    features: RiseDeclineFeatures | None


class ParsnipRiseDeclineResult(BaseModel):
    fitdatainfo: FitDataInfo
    risedeclinefeatures: RiseDeclineFeatures
    classifications: list[Classifications]
