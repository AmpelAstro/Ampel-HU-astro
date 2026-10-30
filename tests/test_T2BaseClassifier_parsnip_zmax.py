"""Regression tests for clamping parsnip redshifts to the model band-weight grid.

Parsnip interpolates precomputed per-band weights on a fixed grid and reads the
shifted index both at ``int_locs`` and ``int_locs + 1``. Above the grid end this
raises an IndexError that ``run_parsnip_zsample`` does not catch (it only handles
ValueError), which aborts the T2 worker. ``get_parsnip_class`` therefore clamps
the sampled redshifts to the ceiling reported by ``parsnip_model_zmax``.

These tests exercise the geometry helper and the clamping contract without
loading a real parsnip model, Torch, or any network resources.
"""

from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from ampel.contrib.hu.model.ParsnipRiseDeclineResult import ParsnipFailure
from ampel.contrib.hu.t2 import T2BaseClassifier as mod
from ampel.contrib.hu.t2.T2BaseClassifier import T2BaseClassifier, parsnip_model_zmax

# Geometry of the model in use (sn+2ulens+dwarfs / parsnipModelAug10Pct30.h5).
# Above this grid the band-weight interpolation index runs off the tensor:
#   n=25485, loc_max=15249.0, spacing=6.829252e-05  ->  zmax=3.99926
REAL_N = 25485
REAL_LOC_MAX = 15249.0
REAL_SPACING = 6.829252e-05


def make_model(n: int, loc_max: float, spacing: float) -> SimpleNamespace:
    """Stand-in parsnip model exposing only the band-weight geometry."""
    return SimpleNamespace(
        band_interpolate_weights=np.zeros((6, n)),
        band_interpolate_locations=np.array([0.0, loc_max]),
        band_interpolate_spacing=spacing,
    )


def make_classifier(model: SimpleNamespace) -> T2BaseClassifier:
    """Bare T2BaseClassifier with only the attributes get_parsnip_class uses."""
    inst = T2BaseClassifier.__new__(T2BaseClassifier)
    inst._class_parsnip = {"m": {"model": model, "classifier": object()}}
    inst.parsnipplot_properties = None
    inst.parsnip_zeropoint_offset = 0.0
    return inst


# ---------------------------------------------------------------------------
# parsnip_model_zmax: pure tensor geometry
# ---------------------------------------------------------------------------


def test_zmax_matches_real_model_geometry():
    """Reproduce the ceiling measured on the affected model."""
    zmax = parsnip_model_zmax(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    assert zmax == pytest.approx(3.99926, abs=1e-4)


def test_zmax_brackets_observed_crash_boundary():
    """The ceiling sits just below 4.00; z=4.407 (the reported alert) is above it.

    Empirically z=4.00 was OK and z=4.01 crashed; the exact ceiling (3.99926)
    lies in between, so it must reject 4.01/4.407 while not reaching 4.00 by
    more than the grid resolution.
    """
    zmax = parsnip_model_zmax(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    assert 3.99 < zmax < 4.01
    assert zmax < 4.01  # the observed crash at z=4.01 is rejected
    assert zmax < 4.407  # so is the reported alert


def test_zmax_leaves_room_for_plus_one_neighbour():
    """At zmax the shifted index plus the +1 neighbour stays in bounds."""
    n, loc_max, spacing = REAL_N, REAL_LOC_MAX, REAL_SPACING
    zmax = parsnip_model_zmax(make_model(n, loc_max, spacing))
    # Worst-case shifted index, mirroring parsnip's own interpolation maths.
    int_locs = int(loc_max + np.log10(1 + zmax) / spacing)
    assert int_locs + 1 <= n - 1


def test_zmax_degenerate_inputs_return_zero():
    """A grid with no headroom (or bad spacing) yields 0.0, not a negative z."""
    assert parsnip_model_zmax(make_model(10, 10.0, 1e-4)) == 0.0  # budget <= 0
    assert parsnip_model_zmax(make_model(100, 1.0, 0.0)) == 0.0  # spacing <= 0
    assert parsnip_model_zmax(make_model(100, 1.0, -1e-4)) == 0.0  # spacing < 0


def test_zmax_scales_with_grid_size():
    """A larger grid permits a higher redshift ceiling."""
    small = parsnip_model_zmax(make_model(1000, 100.0, REAL_SPACING))
    large = parsnip_model_zmax(make_model(25485, REAL_LOC_MAX, REAL_SPACING))
    assert 0.0 < small < large


# ---------------------------------------------------------------------------
# get_parsnip_class: clamping contract
# ---------------------------------------------------------------------------


def test_get_parsnip_class_drops_out_of_range_samples():
    """Out-of-range z are filtered and zsample/zweights stay in sync."""
    inst = make_classifier(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    sentinel = object()
    with patch.object(mod, "run_parsnip_zsample", return_value=sentinel) as rpz:
        out = inst.get_parsnip_class(
            "m",
            None,
            zsample=[0.5, 4.407, 1.2],
            zweights=[0.2, 0.3, 0.5],
        )

    assert out is sentinel
    zs = rpz.call_args.args[1]
    zw = rpz.call_args.args[2]
    assert zs == [0.5, 1.2]
    assert zw == [0.2, 0.5]
    assert rpz.call_args.kwargs["model"] is inst._class_parsnip["m"]["model"]


def test_get_parsnip_class_all_out_of_range_returns_failure():
    """When nothing survives the clamp, a ParsnipFailure is returned instead."""
    inst = make_classifier(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    with patch.object(mod, "run_parsnip_zsample") as rpz:
        out = inst.get_parsnip_class("m", None, zsample=[4.407, 5.0, 6.0])

    rpz.assert_not_called()
    assert isinstance(out, ParsnipFailure)
    assert out.model == "m"
    assert out.failed == "RedshiftOutOfRange"


def test_get_parsnip_class_unweighted_defaults_applied():
    """Without zweights, surviving samples get uniform 1/len(zsample) weights."""
    inst = make_classifier(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    with patch.object(mod, "run_parsnip_zsample", return_value=object()) as rpz:
        inst.get_parsnip_class("m", None, zsample=[0.5, 1.2, 4.407])

    zw = rpz.call_args.args[2]
    # 1/len(original zsample) == 1/3, not 1/len(surviving) == 1/2.
    assert zw == [pytest.approx(1 / 3), pytest.approx(1 / 3)]


def test_get_parsnip_class_in_range_passes_everything_through():
    """All-in-range input is forwarded untouched (no spurious clamping)."""
    inst = make_classifier(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    with patch.object(mod, "run_parsnip_zsample", return_value=object()) as rpz:
        inst.get_parsnip_class("m", None, zsample=[0.1, 1.0, 2.5])

    assert rpz.call_args.args[1] == [0.1, 1.0, 2.5]


def test_get_parsnip_class_mismatched_weights_raise():
    """zip(..., strict=True) surfaces a length mismatch rather than truncating."""
    inst = make_classifier(make_model(REAL_N, REAL_LOC_MAX, REAL_SPACING))
    with pytest.raises(ValueError):
        inst.get_parsnip_class("m", None, zsample=[0.5, 1.0], zweights=[1.0])
