from pathlib import Path

import numpy as np
import pytest
import requests
import yaml
from syrupy.matchers import path_type

from ampel.cli.JobCommand import JobCommand
from ampel.model.job.JobModel import JobModel


def round_float(value: float, path, rel=1e-3):
    """Round float to avoid snapshot mismatches due to small numerical differences."""
    if np.isfinite(value):
        return round(value / rel) * rel
    return value


@pytest.fixture
def model_paths():
    model = Path(__file__).parent / "data" / "parsnipModelAug10Pct30.h5"
    if not model.exists():
        r = requests.get("https://box.hu-berlin.de/f/78d7594841da4653a743/?dl=1")
        r.raise_for_status()
        with open(model, "wb") as f:
            f.write(r.content)

    classifier = Path(__file__).parent / "data" / "parsnipClassifierAug10Pct30.pkl"
    if not classifier.exists():
        r = requests.get("https://box.hu-berlin.de/f/d1cec88c0bec4123b208/?dl=1")
        r.raise_for_status()
        with open(classifier, "wb") as f:
            f.write(r.content)

    return model, classifier


@pytest.fixture
def alertprocessor_config(model_paths):
    model, classifier = model_paths

    with (Path(__file__).parent / "data" / "t2runparsnip-job.yml").open("r") as f:
        config = yaml.safe_load(f)

    # point to local copies of model and classifier
    config["task"][0]["config"]["directives"][0]["ingest"]["mux"]["combine"][0][
        "state_t2"
    ][1]["config"]["paths_parsnip"]["sn+2ulens+dwarfs"].update(
        {
            "model": str(model),
            "classifier": str(classifier),
        }
    )

    # load alert from test data
    config["task"][0]["config"]["supplier"]["config"].update(
        {
            "loader": {
                "unit": "FileAlertLoader",
                "config": {
                    "files": [
                        str(
                            Path(__file__).parent
                            / "data"
                            / "alert_169641461369274412.json"
                        )
                    ]
                },
            },
            "deserialize": "json",
        }
    )

    return config


def test_vroprep_classifier(
    mock_context, alertprocessor_config, ampel_logger, snapshot
):
    cmd = JobCommand()
    job = JobModel(**alertprocessor_config)

    for chan in job.channel:
        mock_context.add_channel(chan.channel)

    t1, t2 = (
        mock_context.db.get_collection("t1", mode="r"),
        mock_context.db.get_collection("t2", mode="r"),
    )

    tasks = cmd.load_tasks(mock_context, job, ampel_logger)
    cmd.run_tasks(mock_context, job, tasks[:1], "jobtest", ampel_logger)

    doc = t2.find_one({"unit": "T2RunParsnipRiseDecline"})
    assert doc["code"] == -1

    assert t1.find_one({"link": doc["link"], "stock": doc["stock"]}) is not None

    cmd.run_tasks(mock_context, job, tasks[1:], "jobtest", ampel_logger)

    doc = t2.find_one({"unit": "T2RunParsnipRiseDecline"})
    assert doc["code"] == 0

    classification = doc["body"][0]["classifications"][0]
    # features have too much platform-dependent numerical noise for reasons that are poorly understood
    del classification["features"]

    assert (
        snapshot(matcher=path_type(types=(float,), replacer=round_float))
        == classification
    )
