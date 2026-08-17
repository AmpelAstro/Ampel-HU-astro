import warnings

from ampel.lsst.t0.DecentVroFilter import RUBIN_ALERT_FLAGS, DecentVroFilter

__all__ = ["RUBIN_ALERT_FLAGS", "DecentVroFilter"]

warnings.warn(
    "ampel.contrib.hu.t0.DecentVroFilter is deprecated and will be removed in a future release. Please use ampel.lsst.t0.DecentVroFilter instead",
    stacklevel=2,
)
