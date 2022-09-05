
from ampel.config.builder.DistConfigBuilder import DistConfigBuilder
from ampel.config.builder.DisplayOptions import DisplayOptions

def test_build_config() -> None:
    cb = DistConfigBuilder(DisplayOptions(verbose=True, debug=True))
    cb.load_distributions()
    config = cb.build_config(
        stop_on_errors=2,
        config_validator="ConfigValidator",
        get_unit_env=False,
    )