import pytest

from ampel.content.DataPoint import DataPoint
from ampel.content.StockDocument import StockDocument
from ampel.enum.DocumentCode import DocumentCode
from ampel.view.T2DocView import TYPE_STATE_T2, T2DocView
from ampel.view.TransientView import TransientView
from ampel.ztf.util.ZTFIdMapper import to_ampel_id


@pytest.fixture()
def t3_transient_views() -> list[TransientView]:
    return [
        TransientView(
            id=(stock_id := to_ampel_id(name)),
            stock=StockDocument(
                {"stock": stock_id, "channel": ["CHANNYCHAN"]},
            ),
            t0=[
                DataPoint(id=i, stock=stock_id, body={}, channel=[], meta=[])
                for i in range(10)
            ],
            t2=[
                T2DocView(
                    stock=stock_id,
                    link=0,
                    tag=[],
                    code=DocumentCode.OK,
                    meta=[{"code": DocumentCode.OK, "tier": 2}],
                    t2_type=TYPE_STATE_T2,
                    unit="T2LightCurveSummary",
                    confid=None,
                    body=[{"len": 10}],
                ),
                T2DocView(
                    stock=stock_id,
                    link=0,
                    tag=[],
                    code=DocumentCode.OK,
                    meta=[{"code": DocumentCode.OK, "tier": 2}],
                    t2_type=TYPE_STATE_T2,
                    unit="T2SNCosmo",
                    confid=None,
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
