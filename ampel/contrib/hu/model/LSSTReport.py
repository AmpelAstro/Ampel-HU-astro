from collections.abc import Sequence
from typing import Annotated

from pydantic import Field

from ampel.base.AmpelBaseModel import AmpelBaseModel


class PhotometricPoint(AmpelBaseModel):
    """Observed photometric point in an LSST alert report"""

    time: Annotated[float, Field(description="epoch of observation")]
    flux: Annotated[float, Field(description="observed flux")]
    fluxerr: Annotated[float, Field(description="flux uncertainty")]
    band: Annotated[str, Field(description="photometric band")]
    zp: Annotated[float, Field(description="zero point")]
    zpsys: Annotated[str, Field(description="zero point system")]


class Object(AmpelBaseModel):
    """Object associated with an LSST alert report (diaObject)"""

    id: Annotated[int, Field(description="diaObjectId")]
    ra: Annotated[float, Field(description="right ascension (deg)")]
    ra_err: Annotated[float, Field(description="right ascension uncertainty (deg)")]
    dec: Annotated[float, Field(description="declination")]
    dec_err: Annotated[float, Field(description="declination uncertainty")]
    ra_dec_cov: Annotated[
        float, Field(description="right ascension/declination covariance")
    ]
    # redshift: float
    # redshift_err: float


class ModelClassification(AmpelBaseModel):
    model: Annotated[str, Field(description="name of the model")]
    probabilities: Annotated[
        dict[str, float], Field(description="probabilities for each class")
    ]


class Classification(AmpelBaseModel):
    name: Annotated[str, Field(description="classifier name")]
    version: Annotated[str, Field(description="classifier version")]
    info: Annotated[str | None, Field(description="additional information about the classifier")]
    models: Annotated[
        list[ModelClassification], Field(description="list of model classifications")
    ]


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
    """
    Data model for LSST alert reports from Ampel.
    """

    object: Object
    photometry: Sequence[PhotometricPoint]
    classification: list[Classification] = []
    host: list[Host] = []
    features: list[Feature] = []
