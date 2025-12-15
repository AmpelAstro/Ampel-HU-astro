from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.core.AmpelContext import AmpelContext
from ampel.model.UnitModel import UnitModel
from ampel.struct.T3Store import T3Store
from ampel.view.TransientView import TransientView


def test_TransientTablePublisher(
    mock_context: AmpelContext,
    t3_transient_views: list[TransientView],
    ampel_logger,
):
    unit = mock_context.loader.new_logical_unit(
        UnitModel(
            unit="TransientTablePublisher",
            config={
                "transient_table_schema": {
                    "T2LightCurveSummary": {
                        "len": "len",
                    },
                    "T2RunSNCosmo": {
                        "dof": "sncosmo_info.ndof",
                        "t0": "fit_results.t0",
                    },
                },
                "convert_stock_to": "ztf",
            },
        ),
        sub_type=AbsPhotoT3Unit,
        logger=ampel_logger,
    )

    result = unit.process(iter(t3_transient_views), T3Store())  # type: ignore[arg-type]
    assert result is None
