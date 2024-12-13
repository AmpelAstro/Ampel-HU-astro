from unittest.mock import Mock, patch

import pytest

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo
from ampel.view.T2DocView import T2DocView


@pytest.fixture
def mock_t2runsncosmo():
    return T2RunSncosmo(
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
    )


def test_T2RunSncosmo_no_redshift(mock_t2runsncosmo):
    compound = Mock(spec=T1Document)
    datapoints = [Mock(spec=DataPoint)]
    t2_views = [Mock(spec=T2DocView)]

    result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["z_source"] is None


def test_T2RunSncosmo_no_phaselimit(mock_t2runsncosmo):
    compound = Mock(spec=T1Document)
    datapoints = [Mock(spec=DataPoint)]
    t2_views = [Mock(spec=T2DocView)]

    with patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["jdstart"] is None


def test_T2RunSncosmo_fit_error(mock_t2runsncosmo):
    compound = Mock(spec=T1Document)
    datapoints = [Mock(spec=DataPoint)]
    t2_views = [Mock(spec=T2DocView)]

    with (
        patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")),
        patch.object(mock_t2runsncosmo, "_get_phaselimit", return_value=(0, 1000)),
        patch.object(mock_t2runsncosmo, "get_flux_table", return_value=Mock()),
        patch("sncosmo.fit_lc", side_effect=RuntimeError("fit error")),
    ):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert result["run_error"] is True


def test_T2RunSncosmo_success(mock_t2runsncosmo):
    compound = Mock(spec=T1Document)
    datapoints = [Mock(spec=DataPoint)]
    t2_views = [Mock(spec=T2DocView)]

    mock_flux_table = Mock()
    mock_fit_result = {
        "parameters": [0.1],
        "data_mask": [True],
        "covariance": None,
        "param_names": ["z"],
        "chisq": 1.0,
        "ndof": 1,
    }
    mock_fitted_model = Mock()

    with (
        patch.object(mock_t2runsncosmo, "_get_redshift", return_value=(0.1, "Fixed")),
        patch.object(mock_t2runsncosmo, "_get_phaselimit", return_value=(0, 1000)),
        patch.object(mock_t2runsncosmo, "get_flux_table", return_value=mock_flux_table),
        patch("sncosmo.fit_lc", return_value=(mock_fit_result, mock_fitted_model)),
        patch.object(
            mock_t2runsncosmo, "_get_fit_metrics", return_value={"metric": 1.0}
        ),
    ):
        result = mock_t2runsncosmo.process(compound, datapoints, t2_views)

    assert isinstance(result, dict)
    assert "sncosmo_result" in result
    assert result["fit_metrics"] == {"metric": 1.0}
