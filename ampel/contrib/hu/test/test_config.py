from ampel.config.builder.DisplayOptions import DisplayOptions
from ampel.config.builder.DistConfigBuilder import DistConfigBuilder


def test_build_config() -> None:
    cb = DistConfigBuilder(DisplayOptions(verbose=True, debug=True))
    cb.load_distributions(
        prefixes=[
            "ampel-interface",
            "ampel-core",
            "ampel-alerts",
            "ampel-photometry",
            "ampel-ztf",
            "ampel-hu-astro",
        ]
    )
    config = cb.build_config(
        stop_on_errors=2,
        config_validator="ConfigValidator",
        get_unit_env=False,
    )
