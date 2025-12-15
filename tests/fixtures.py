import pytest

from ampel.content.DataPoint import DataPoint
from ampel.content.StockDocument import StockDocument
from ampel.content.T1Document import T1Document
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
                DataPoint(
                    id=i,
                    stock=stock_id,
                    body={"ra": 298.8424232, "dec": 44.4},
                    channel=[],
                    meta=[],
                )
                for i in range(10)
            ],
            t1=[
                T1Document(
                    unit="T1LightCurveBuilder",
                    config=0,
                    stock=stock_id,
                    link=0,
                    dps=list(range(10)),
                    channel=["CHANNYCHAN"],
                    tag=[],
                    code=DocumentCode.OK,
                    meta=[{"code": DocumentCode.OK, "tier": 1}],
                )
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
