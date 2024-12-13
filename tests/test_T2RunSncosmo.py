from unittest.mock import MagicMock, Mock, patch

import numpy as np
import pytest

from ampel.base.AuxUnitRegister import AuxUnitRegister
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo
from ampel.log.AmpelLogger import AmpelLogger
from ampel.model.UnitModel import UnitModel
from ampel.view.T2DocView import T2DocView
from ampel.ztf.view.ZTFT2Tabulator import ZTFT2Tabulator


@pytest.fixture
def mock_t2runsncosmo():
    AuxUnitRegister._dyn["ZTFT2Tabulator"] = ZTFT2Tabulator  # noqa: SLF001

    t2 = T2RunSncosmo(
        sncosmo_model_name="salt2",
        redshift_kind=None,
        fixed_z=None,
        backup_z=None,
        scale_z=None,
        sncosmo_bounds={},
        apply_mwcorrection=False,
        phaseselect_kind=None,
        noisified=False,
        plot_db=False,
        plot_props=None,
        plot_matplotlib_suffix=None,
        plot_matplotlib_dir=".",
        t2_dependency=[],
        tabulator=[UnitModel(unit="ZTFT2Tabulator")],
        logger=AmpelLogger.get_logger(),
    )
    t2.post_init()
    return t2


def inputs():
    compound = MagicMock(spec=T1Document)
    base = {"tag": ["ZTF"]}
    datapoints = [
        base | {"body": d}
        for d in (
            {"jd": 2450000, "mag": 20, "magerr": 0.1},
            {"jd": 2450001, "mag": 21, "magerr": 0.1},
            {"jd": 2450002, "mag": 22, "magerr": 0.1},
        )
    ]
    t2_views = [MagicMock(spec=T2DocView)]
    return compound, datapoints, t2_views


def test_T2RunSncosmo_no_redshift(mock_t2runsncosmo):
    compound, datapoints, t2_views = inputs()

    result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["z_source"] == "Fitted"


def test_T2RunSncosmo_no_phaselimit(mock_t2runsncosmo):
    compound = MagicMock(spec=T1Document)
    datapoints = [MagicMock(spec=DataPoint)]
    t2_views = [MagicMock(spec=T2DocView)]

    with patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["jdstart"] == -np.inf


def test_T2RunSncosmo_fit_error(mock_t2runsncosmo):
    compound, datapoints, t2_views = inputs()

    with (
        patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")),
        patch.object(mock_t2runsncosmo, "_get_phaselimit", return_value=(0, 1000)),
        patch("sncosmo.fit_lc", side_effect=RuntimeError("fit error")),
    ):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["run_error"] is True


def test_T2RunSncosmo_success(mock_t2runsncosmo):
    compound, datapoints, t2_views = inputs()

    mock_fit_result = {
        "parameters": np.array([0.1]),
        "data_mask": np.array([True]),
        "covariance": None,
        "param_names": ["z"],
        "chisq": 1.0,
        "ndof": 1,
    }
    mock_fitted_model = Mock()

    with (
        patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")),
        patch.object(mock_t2runsncosmo, "_get_phaselimit", return_value=(0, 1000)),
        patch("sncosmo.fit_lc", return_value=(mock_fit_result, mock_fitted_model)),
        patch.object(
            mock_t2runsncosmo, "_get_fit_metrics", return_value={"metric": 1.0}
        ),
    ):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert "sncosmo_result" in result
    assert result["fit_metrics"] == {"metric": 1.0}
