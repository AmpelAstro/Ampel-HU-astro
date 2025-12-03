from collections.abc import Mapping, Sequence
from typing import Any, Literal

from ampel.contrib.hu.model.LSSTReport import (
    Classification,
    Feature,
    Host,
    ModelClassification,
)
from ampel.contrib.hu.t2.T2LSSTReport import T2LSSTReport


def get_parsnip_summary(full_classification: dict, summary_mode: str) -> dict:
    """
    Collapse a full ELAsTiCC classification matrix until a subset
    of (possibly combined) probabilities.
    """

    outclass = {}
    if summary_mode == "elasticc_simple":
        # Simple ones
        for singleclass in ["KN", "SLSN-I", "TDE"]:
            outclass["P(" + singleclass + ")"] = full_classification[singleclass]
        # Combined SNIa prop
        outclass["P(SNIa)"] = sum(
            [full_classification[s] for s in ["SNIa", "SNIa91bg", "SNIax"]]
        )
        # Combined SN prop
        outclass["P(SN)"] = sum(
            [full_classification[s] for s in ["SLSN-I", "PISN", "SNII", "SNibc"]]
        )
        outclass["P(SN)"] += outclass["P(SNIa)"]
    elif summary_mode == "full":
        outclass = {
            'P('+class_name+')': prob for class_name, prob in full_classification.items()
        }

    return outclass


class T2ClassificationReport(T2LSSTReport):
    """
    Require and propagate classification information from ParsnipRiseDecline.
    """

    report_t2s: Sequence[str] = ["T2RunParsnipRiseDecline"]
    summary_mode: Literal["elasticc_simple", "full"] = "full"

    # Only transients with features within these limits will be reported
    risedecline_select_map: dict[str, Sequence[float]] = {
        "t_lc": [5, 400],
        "ampel_dist": [0.0, 10.0],
        "redshift": [0, 0.5],
        "parsnip_modelchidof": [0.0, 4],
        "model_dof": [3, 99],
    }

    def process_t2s(
        self, report_views: dict[str, Mapping[str, Any]]
    ) -> Sequence[Classification | Host | Feature] | None:
        """
        Process T2 views to extract information to be propagated.
        Will parse results from T2ParsnipRiseDecline and:
        - Check whether minimal quality criteria are met.
        - Collect sparse feature and host information.
        - Look for Parsnip classification results and propagate a summary.
        """

        body = report_views["T2RunParsnipRiseDecline"]
        # Check whether lightcurve quality limits are met
        if "risedeclinefeatures" not in body:
            return None
        features = {}
        for field, limits in self.risedecline_select_map.items():
            if (
                field in body["risedeclinefeatures"]
                and limits[0] < body["risedeclinefeatures"][field] < limits[1]
            ):
                features[field] = body["risedeclinefeatures"][field]
                continue
            return None

        # To be propagated, collect information.
        unitresults: list[Classification | Host | Feature] = []

        # Add selection features. Todo: add a selection of other interesting ones.
        unitresults.append(
            Feature(name="RiseDeclineFeatures", version="1", features=features)
        )

        # At the moment, only the base parsnip is run, so only looking for this
        if "parsnip" in body["classifications"][0]:
            unitresults.append(
                Classification(
                    name="ParsnipRiseDecline",
                    version="1",
                    models=[
                        ModelClassification(
                            model=modelname,
                            probabilities=get_parsnip_summary(
                                classdoc["classification"], summary_mode=self.summary_mode, 
                            ),
                        )
                        for modelname, classdoc in body["classifications"][0][
                            "parsnip"
                        ].items()
                    ],
                )
            )

        # Add host information if available
        # (should be unless the classifier is run in a non-redshift mode)
        if "fitdatainfo" in body and "redshift" in features:
            unitresults.append(
                Host(
                    name="AmpelDigestRedshifts",
                    redshift=features["redshift"],
                    distance=features["ampel_dist"],
                    source="AmpelZ",
                    info=body["fitdatainfo"].get("z_source", None),
                )
            )

        return unitresults
