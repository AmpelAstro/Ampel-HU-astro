from collections.abc import Sequence

from ampel.base.AmpelBaseModel import AmpelBaseModel


class PhotometricPoint(AmpelBaseModel):
    time: float
    flux: float
    fluxerr: float
    band: str
    zp: float
    zpsys: str


class Object(AmpelBaseModel):
    id: int
    ra: float
    ra_err: float
    dec: float
    dec_err: float
    ra_dec_cov: float
    # redshift: float
    # redshift_err: float


class ModelClassification(AmpelBaseModel):
    model: str
    probabilities: dict[str, float]


class Classification(AmpelBaseModel):
    name: str
    version: str
    info: str | None = None
    models: list[ModelClassification]


class Host(AmpelBaseModel):
    name: str
    source: str
    redshift: float
    redshift_error: float | None = None
    distance: float
    info: str | None = None


class Feature(AmpelBaseModel):
    name: str
    version: str
    info: str | None = None
    features: dict[str, float]


class LSSTReport(AmpelBaseModel):
    object: Object
    photometry: Sequence[PhotometricPoint]
    classification: list[Classification] = []
    host: list[Host] = []
    features: list[Feature] = []
