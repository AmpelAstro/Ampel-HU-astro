import pytest
from ampel.log.AmpelLogger import AmpelLogger
from ampel.contrib.hu.t3.SlackSummaryPublisher import SlackSummaryPublisher
from ampel.secret.NamedSecret import NamedSecret

from ampel.content.StockDocument import StockDocument
from ampel.content.DataPoint import DataPoint
from ampel.content.T2Document import T2Document
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ampel_id

import requests
from slack import WebClient
import csv
from io import StringIO


@pytest.fixture
def t3_transient_views() -> list[TransientView]:
    return [
        TransientView(
            id=(stock_id := to_ampel_id(name)),
            stock=StockDocument(
                {"stock": stock_id, "channel": ["CHANNYCHAN"]},
            ),
            t0=[DataPoint(id=i, stock=stock_id, body={}) for i in range(10)],
            t2=[
                T2Document(stock=stock_id, unit="T2LightCurveSummary", body=[{"len": 10}]),
                T2Document(
                    stock=stock_id,
                    unit="T2SNCosmo",
                    body=[
                        {
                            "model": "salt2",
                            "sncosmo_info": {"ndof": 42},
                            "fit_results": {"t0": 13480238.2309423},
                        }
                    ],
                ),
            ],
        )
        for name in [
            "ZTF18actmutj",
            "ZTF20abthfto",
            "ZTF21aabltvh",
        ]
    ]


def test_slacksummary(t3_transient_views, mocker):

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
            "slack_token": NamedSecret(label="tokeytoke", value="xoxoxox"),
            "slack_channel": "#ampel-live",
            "full_photometry": True,
            "logger": AmpelLogger.get_logger(),
        }
    )

    assert len(t3_transient_views) < unit.excitement["Low"], "Small number passed"

    # intercept Slack API calls
    mocker.patch("requests.post")
    mocker.patch("slack.WebClient.api_call")

    unit.process(iter(t3_transient_views))

    api_call = WebClient.api_call
    api_call.assert_called_once()
    assert (
        "MEH!" in api_call.call_args[1]["json"]["text"]
    ), "Text matches number of transients selected"

    requests.post.assert_called()
    assert len(requests.post.call_args_list) == 2, "2 explicit requests issued"

    # Verify summary
    content = requests.post.call_args_list[0][1]["files"]["file"]
    with StringIO(content) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # verify that T2 information is in summary
    t2s = set(t2["unit"] for tv in t3_transient_views for t2 in tv.t2)
    assert len(t2s) > 0

    # Verify that nested t2 results were extrected
    # This tests assumes that sncosmo was run on the test data
    for key in ["T2-model", "T2-sncosmo_info_ndof", "T2-fit_results_t0"]:
        assert key in reader.fieldnames

    assert len(rows) == len(t3_transient_views), "1 row per transient"

    # Verify photometry dump
    content = requests.post.call_args_list[1][1]["files"]["file"]
    with StringIO(content) as f:
        rows = list(csv.reader(f))
    assert (
        len(rows) == sum(len(v.t0) for v in t3_transient_views) + 1
    ), "1 row per photopoint"
