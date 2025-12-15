import csv
from io import StringIO

import requests
from slack_sdk import WebClient

from ampel.contrib.hu.t3.SlackSummaryPublisher import SlackSummaryPublisher
from ampel.log.AmpelLogger import AmpelLogger
from ampel.struct.T3Store import T3Store
from ampel.view.TransientView import TransientView


def test_slacksummary(t3_transient_views: list[TransientView], mocker):
    unit = SlackSummaryPublisher(
        **{
            "cols": [
                "ztf_name",
                "ra",
                "dec",
                "magpsf",
                "sgscore1",
                "rb",
                "most_recent_detection",
                "first_detection",
                "n_detections",
                "distnr",
                "distpsnr1",
                "isdiffpos",
                "_id",
            ],
            "excitement": {"Low": 50, "Mid": 200, "High": 400},
            "slack_token": dict(label="tokeytoke", value="xoxoxox"),
            "slack_channel": "#ampel-live",
            "full_photometry": True,
            "logger": AmpelLogger.get_logger(),
        }
    )

    assert len(t3_transient_views) < unit.excitement["Low"], "Small number passed"

    # intercept Slack API calls
    mocker.patch("requests.post")
    mocker.patch("slack_sdk.WebClient.api_call")

    unit.process(iter(t3_transient_views), T3Store())  # type: ignore[arg-type]

    api_call = WebClient.api_call
    api_call.assert_called_once()  # type: ignore[attr-defined]
    assert (
        "MEH!" in api_call.call_args[1]["json"]["text"]  # type: ignore[attr-defined]
    ), "Text matches number of transients selected"

    requests.post.assert_called()  # type: ignore[attr-defined]
    assert len(requests.post.call_args_list) == 2, "2 explicit requests issued"  # type: ignore[attr-defined]

    # Verify summary
    content = requests.post.call_args_list[0][1]["files"]["file"]  # type: ignore[attr-defined]
    with StringIO(content) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # verify that T2 information is in summary
    t2s = set(t2.unit for tv in t3_transient_views for t2 in (tv.t2 or []))
    assert len(t2s) > 0

    # Verify that nested t2 results were extrected
    # This tests assumes that sncosmo was run on the test data
    assert reader.fieldnames
    for key in ["T2-model", "T2-sncosmo_info_ndof", "T2-fit_results_t0"]:
        assert key in reader.fieldnames

    assert len(rows) == len(t3_transient_views), "1 row per transient"

    # Verify photometry dump
    content = requests.post.call_args_list[1][1]["files"]["file"]  # type: ignore[attr-defined]
    with StringIO(content) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == sum(len(v.t0 or []) for v in t3_transient_views), (
        "1 row per photopoint"
    )
