from pathlib import Path

import pytest
import yaml

from ampel.util.template import apply_templates


@pytest.fixture
def data_path():
    return Path(__file__).parent / "data" / "template-job"


def dump_without_anchors(data: dict):
    class NoAnchorDumper(yaml.SafeDumper):
        def ignore_aliases(self, data) -> bool:
            return True

    return yaml.dump(data, Dumper=NoAnchorDumper)


@pytest.fixture
def original_dict(data_path):
    with (data_path / "original.yaml").open() as f:
        return yaml.safe_load(f)


@pytest.fixture(params=["template", "template-multi"])
def template_dict(request, data_path):
    with (data_path / f"{request.param}.yaml").open() as f:
        return yaml.safe_load(f)


def test_consumer_template(mock_context, original_dict, template_dict, ampel_logger):
    transformed = apply_templates(
        mock_context, ["ingest_lsst_alerts"], target=template_dict, logger=ampel_logger
    )
    transformed.pop("template")  # remove template key for comparison
    assert dump_without_anchors(transformed) == dump_without_anchors(original_dict), (
        "template resolves to the original verbose config"
    )
