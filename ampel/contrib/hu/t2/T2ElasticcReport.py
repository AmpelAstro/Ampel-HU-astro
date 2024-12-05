#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2ElasticcReport.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                25.08.2022
# Last Modified Date:  13.12.2022
# Last Modified By:    jnordin@physik.hu-berlin.de

from collections import defaultdict
from collections.abc import Sequence
from typing import Literal

import numpy as np

from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.util import get_payload
from ampel.enum.DocumentCode import DocumentCode
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.struct.JournalAttributes import JournalAttributes
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView

# Tying classification output classes with ELAsTICC taxonomy classes
parsnip_taxonomy = {
    # extracted from https://github.com/LSSTDESC/elasticc/blob/bc0de488c5276ce61b650117db19e93634b10815/taxonomy/taxonomy.ipynb
    #    'SLSN':  2242,   # never used?
    "SNII": 2224,
    "SNIa": 2222,
    "SNibc": 2223,
    "SNIbc": 2223,
    "TDE": 2243,
    "CART": 2245,
    "ILOT": 2244,
    "Mdwarf-flare": 2233,
    "PISN": 2246,
    "KN": 2232,
    "SLSN-I": 2242,
    "SNIa91bg": 2226,
    "SNIax": 2225,
    "dwarf-nova": 2234,
    # mssing
    "uLens": 2235,
}

# T2XgbClassifier can also yield direct evaluations for some cases
# when the main run fails. We here map these to elasticc codes
direct_evaluations = {"AGN": 2330, "uLens": 2235, "EB": 2320}


# Prior section - can be applied as demo, but not matching the simulations.

# Elasticc redshift priors.
# Warning: Based on normalized distributions of training sample, no rate info taken into account.
# Bins based on z 0.1 bins starting at -0.09 and extending to 3.4
# Probabilities assumed to add up to 1000 for each redshift bin.
zmap = {
    "CART": {
        0: 9,
        1: 124,
        2: 135,
        3: 143,
        4: 143,
        5: 127,
        6: 97,
        7: 62,
        8: 33,
        9: 16,
        10: 8,
        11: 6,
        12: 6,
        13: 5,
        14: 5,
        15: 6,
        16: 7,
        17: 6,
        18: 3,
        19: 1,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "uLens": {
        0: 304,
        1: 0,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,
        7: 0,
        8: 0,
        9: 0,
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SLSN-I": {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 5,
        5: 11,
        6: 20,
        7: 35,
        8: 55,
        9: 80,
        10: 112,
        11: 153,
        12: 215,
        13: 307,
        14: 421,
        15: 536,
        16: 645,
        17: 746,
        18: 827,
        19: 879,
        20: 911,
        21: 936,
        22: 958,
        23: 976,
        24: 988,
        25: 995,
        26: 998,
        27: 999,
        28: 999,
        29: 999,
        30: 999,
        31: 999,
        32: 999,
        33: 999,
        34: 999,
        35: 999,
    },
    "dwarf-nova": {
        0: 304,
        1: 0,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,
        7: 0,
        8: 0,
        9: 0,
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SNIa91bg": {
        0: 5,
        1: 78,
        2: 100,
        3: 127,
        4: 151,
        5: 156,
        6: 132,
        7: 88,
        8: 45,
        9: 17,
        10: 6,
        11: 3,
        12: 2,
        13: 2,
        14: 1,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SNIa": {
        0: 1,
        1: 19,
        2: 26,
        3: 39,
        4: 60,
        5: 92,
        6: 134,
        7: 181,
        8: 217,
        9: 225,
        10: 200,
        11: 150,
        12: 93,
        13: 48,
        14: 23,
        15: 11,
        16: 5,
        17: 2,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "ILOT": {
        0: 24,
        1: 257,
        2: 223,
        3: 173,
        4: 114,
        5: 59,
        6: 23,
        7: 6,
        8: 1,
        9: 0,
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SNIax": {
        0: 4,
        1: 64,
        2: 77,
        3: 95,
        4: 116,
        5: 136,
        6: 145,
        7: 137,
        8: 111,
        9: 77,
        10: 46,
        11: 23,
        12: 9,
        13: 3,
        14: 1,
        15: 1,
        16: 1,
        17: 2,
        18: 1,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "TDE": {
        0: 3,
        1: 46,
        2: 56,
        3: 69,
        4: 84,
        5: 100,
        6: 113,
        7: 123,
        8: 133,
        9: 141,
        10: 145,
        11: 142,
        12: 130,
        13: 108,
        14: 83,
        15: 63,
        16: 50,
        17: 43,
        18: 38,
        19: 34,
        20: 30,
        21: 25,
        22: 20,
        23: 14,
        24: 9,
        25: 4,
        26: 1,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "KN": {
        0: 29,
        1: 293,
        2: 234,
        3: 166,
        4: 98,
        5: 45,
        6: 15,
        7: 3,
        8: 0,
        9: 0,
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SNII": {
        0: 3,
        1: 47,
        2: 59,
        3: 75,
        4: 94,
        5: 112,
        6: 126,
        7: 131,
        8: 127,
        9: 117,
        10: 105,
        11: 95,
        12: 86,
        13: 77,
        14: 70,
        15: 66,
        16: 63,
        17: 59,
        18: 48,
        19: 33,
        20: 21,
        21: 11,
        22: 5,
        23: 2,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "SNibc": {
        0: 4,
        1: 60,
        2: 73,
        3: 89,
        4: 105,
        5: 119,
        6: 127,
        7: 128,
        8: 121,
        9: 107,
        10: 90,
        11: 72,
        12: 53,
        13: 36,
        14: 25,
        15: 18,
        16: 13,
        17: 8,
        18: 3,
        19: 1,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "Mdwarf-flare": {
        0: 304,
        1: 0,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,
        7: 0,
        8: 0,
        9: 0,
        10: 0,
        11: 0,
        12: 0,
        13: 0,
        14: 0,
        15: 0,
        16: 0,
        17: 0,
        18: 0,
        19: 0,
        20: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 0,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
    "PISN": {
        0: 0,
        1: 7,
        2: 10,
        3: 16,
        4: 24,
        5: 38,
        6: 62,
        7: 100,
        8: 152,
        9: 215,
        10: 284,
        11: 352,
        12: 401,
        13: 409,
        14: 368,
        15: 296,
        16: 211,
        17: 131,
        18: 76,
        19: 49,
        20: 36,
        21: 26,
        22: 15,
        23: 6,
        24: 1,
        25: 0,
        26: 0,
        27: 0,
        28: 0,
        29: 0,
        30: 0,
        31: 0,
        32: 0,
        33: 0,
        34: 0,
        35: 0,
    },
}
# Parsnip type priors
# Based on BTS _observed_ rate distributions, i.e. no knowledge of elasticc
# or LSST simulation properties.
# Probabilities assumed to add up to 1000
btsmap = {
    "SNIa": 663,
    "uLens": 10,
    "SLSN-I": 14,
    "dwarf-nova": 10,
    "SNIa91bg": 10,
    "ILOT": 10,
    "SNIax": 10,
    "TDE": 10,
    "KN": 10,
    "SNII": 168,
    "SNIbc": 58,
    "Mdwarf-flare": 10,
    "PISN": 10,
    "CART": 10,
}

# Priors based on the host galaxy u-g color (from elasticc_galcol notebook)
galcol_prior: dict[str, list[float]] = {
    "CART": [1.1681864134021618, 0.5451610139053239],
    "ILOT": [1.1458356003495815, 0.5093055923613805],
    "KN": [1.3089799059907905, 0.3303618787483199],
    "PISN": [0.3426406749881123, 0.4465663913955291],
    #'SLSN': [0.37962520970022234, 0.3723189361658543],
    "SLSN-I": [0.37962520970022234, 0.3723189361658543],
    "SNIa91bg": [1.8926892171924607, 0.4252978932935133],
    "SNIa": [0.6102884945155177, 0.39269590167347995],
    "SNIax": [0.6945950080255028, 0.36977359393920906],
    "SNIbc": [0.9490882652856778, 0.5687292278916529],
    "SNII": [0.8752024059665746, 0.5297978836976114],
    "TDE": [0.8753566579007913, 0.7194828408852172],
    "uLens": [0, 0],
    "Mdwarf-flare": [0, 0],
    "dwarf-nova": [0, 0],
}
# Correlation between galaxy color and mwebv
galcol_mwebv_slope = 1.1582821


class T2ElasticcReport(AbsTiedStateT2Unit):
    """

    Parse a series of T2 results from T2RunParsnip and T2XgbClassifier, and
    create combined classifications according to the taxonomy of
    https://github.com/LSSTDESC/elasticc/blob/main/taxonomy/taxonomy.ipynb.

    This unit assumes there are three XGB and one Parsnip classifier to run.

    """

    # XGB binary models
    # Order taken to mean that "prob0" is the probability that the first
    # term is correct.
    tree_1v2: str = "xgb_v6_simmod_tree12_"
    tree_21v22: str = "xgb_v6_simmod_tree2122_"
    tree_1113v12: str = "xgb_v6_simmod_tree121113_"

    # Setting for report to construct
    broker_name: str = "AMPEL"
    broker_version: str = "v0.4"
    classifier_name: str = "ElasticcLive"
    classifier_version: str = "XGBUnified+Parsnip05"

    # Combine multiple classifiers in report
    multiple_classifiers: bool = False

    # Use a redshift prior from elasticc alerts and obs rate prior from BTS.
    use_priors: bool = False

    # Which units should this be changed to
    t2_dependency: Sequence[
        StateT2Dependency[Literal["T2RunParsnip", "T2MultiXgbClassifier"]]
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def _one_report_per_classifier(report: dict):
        """
        reformat a v0.9 brokerClassification for v0.9.1
        see: https://raw.githubusercontent.com/LSSTDESC/elasticc/5d7b314b537197c99086acf019e6e2c1dc4aa267/alert_schema/elasticc.v0_9_1.brokerClassification.avsc
        """
        by_classifier = defaultdict(list)
        for c in report["classifications"]:
            by_classifier[(c["classifierName"], c["classifierParams"])].append(
                {k: c[k] for k in ("classId", "probability")}
            )
        for (name, params), classifications in by_classifier.items():
            new_report = report.copy()
            new_report["classifications"] = [
                {k: c[k] for k in ("classId", "probability")} for c in classifications
            ]
            new_report["classifierName"] = name
            new_report["classifierParams"] = params
            yield new_report

    def add_zprior(self, parsnip_prob: dict, z: float):
        """
        Adjust probabilities based on redshift prior.
        """
        # Out of range
        if z > 3.409 or z < 0:
            return parsnip_prob

        zbin = int((z + 0.09) / 0.1)

        scaled_prob = {}
        for model_name, zdist in zmap.items():
            scaled_prob[parsnip_taxonomy[model_name]] = (
                parsnip_prob[parsnip_taxonomy[model_name]] * float(zdist[zbin]) / 1000
            )

        # Reweight probabilities
        p = sum([v for v in scaled_prob.values()])
        if not p > 0:
            p = 1.0

        return {k: v / p for k, v in scaled_prob.items()}

    def add_rateprior(self, parsnip_prob: dict):
        """
        Modify fit probabilities based on observed rates.
        """

        scaled_prob = {}
        for model_name, model_prob in btsmap.items():
            scaled_prob[parsnip_taxonomy[model_name]] = (
                parsnip_prob[parsnip_taxonomy[model_name]] * float(model_prob) / 1000
            )

        # Reweight probabilities
        p = sum([v for v in scaled_prob.values()])
        if not p > 0:
            p = 1.0

        return {k: v / p for k, v in scaled_prob.items()}

    def add_hostprior(self, parsnip_prob: dict, host_ug: float | None, mwebv: float):
        """
        Modify fit probabilities based on u-g color, if available.
        """

        if host_ug is None:
            return parsnip_prob

        # Correct for mwebv dependence
        host_ug = host_ug - galcol_mwebv_slope * mwebv

        scaled_prob = {}
        for model_name, model_prob in galcol_prior.items():
            if model_prob[1] > 0:
                scaled_prob[parsnip_taxonomy[model_name]] = parsnip_prob[
                    parsnip_taxonomy[model_name]
                ] * np.exp(
                    -((host_ug - model_prob[0]) ** 2) / (2.0 * model_prob[1] ** 2)
                )
            else:
                scaled_prob[parsnip_taxonomy[model_name]] = 0.0

        # Reweight probabilities (or rescale host properties first? write it out to see...)
        p = sum([v for v in scaled_prob.values()])

        return {k: v / p for k, v in scaled_prob.items()}

    def get_hostcol(self, dia_object: dict, z: float):
        """
        Extract the most relevant u-g color, if present.
        """

        z1 = dia_object.get("hostgal_zphot_q050")
        z2 = dia_object.get("hostgal2_zphot_q050")
        # No z information
        if z1 is None or z2 is None:
            return None

        # Which "mid"fix, '' or '2'?
        midfix = ""
        if abs(z - z2) < abs(z - z1):
            midfix = "2"
        u = dia_object["hostgal" + midfix + "_mag_u"]
        g = dia_object["hostgal" + midfix + "_mag_g"]
        if not u < 40:
            return None
        if not g < 40:
            return None
        return u - g

    def make_unit_result(
        self, compound: T1Document, class_report: dict[str, UBson]
    ) -> UnitResult:
        return UnitResult(
            body={
                "report": class_report,
            },
            code=DocumentCode.OK,
            journal=JournalAttributes(extra={"link": compound["link"]}),
        )

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UnitResult:
        """

        Extract and combine results.

        Parameters
        -----------

        t2_records: List of T2Records from (here) TabulatorRiseDecline.

        Returns
        -------
        dict
        """

        # Construct the base reply. We have some magic here to fill
        class_report: dict[str, UBson] = {
            "brokerName": self.broker_name,
            "brokerVersion": self.broker_version,
        }

        # use an alias variable to inform mypy that classifications is always a list
        class_report["classifications"] = classifications = []

        # Get alert info from T1Document if present
        for metarecord in compound["meta"]:
            if isinstance(alert_id := metarecord.get("alert"), int):
                class_report["alertId"] = alert_id
                class_report["brokerIngestTimestamp"] = metarecord["ts"] * 1000
            if isinstance(alert_ts := metarecord.get("alert_ts"), int):
                class_report["elasticcPublishTimestamp"] = alert_ts

        # Get diaSourceId
        # Get the last diaSource id in datapoints
        for dp in reversed(datapoints):
            if "diaSourceId" in dp["body"]:
                class_report["diaSourceId"] = dp["body"]["diaSourceId"]
                break

        # 1. Obtain base classifications from XGBClassifier and Parsnip
        is1, is21, is1113, parsnip_class, parsnip_z, direct_eval = (
            None,
            None,
            None,
            None,
            None,
            None,
        )

        # Parse t2views - should not be more than one.
        for t2_view in t2_views:
            self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
            # Xgb results either from multiple instances of T2XgbClassifier...
            if t2_view.unit == "T2XgbClassifier":
                t2_res = get_payload(t2_view)
                if "prob0" in t2_res:
                    if t2_res["model"] == self.tree_1v2:
                        is1 = t2_res["prob0"]
                    elif t2_res["model"] == self.tree_21v22:
                        is21 = t2_res["prob0"]
                    elif t2_res["model"] == self.tree_1113v12:
                        is1113 = t2_res["prob0"]
                else:
                    # Direct evaluation available even though XGB did not run
                    direct_eval = t2_res.get("direct_eval", None)
            # ... or all from T2MultiXgbClassifier
            elif t2_view.unit == "T2MultiXgbClassifier":
                t2_res = get_payload(t2_view)
                if t2_res["model"] == "multiXgb":
                    is1 = t2_res["classifications"][self.tree_1v2]["prob0"]
                    is21 = t2_res["classifications"][self.tree_21v22]["prob0"]
                    is1113 = t2_res["classifications"][self.tree_1113v12]["prob0"]
                elif t2_res["model"] == "directEval":
                    direct_eval = t2_res.get("direct_eval", None)
            elif t2_view.unit == "T2RunParsnip":
                t2_res = get_payload(t2_view)

                if "classification" in t2_res:
                    parsnip_class = {
                        parsnip_taxonomy[parsnip_type]: prob
                        for parsnip_type, prob in t2_res["classification"].items()
                    }
                parsnip_z = t2_res.get("z_at_minchi", None)

        # Did XGBClassifier run?
        if is1 is None or is21 is None or is1113 is None:
            classId = 0  # Default for non-reply

            # Replace with  direct evaluation from XGB if present
            if direct_eval is not None:
                classId = direct_evaluations[direct_eval]

            self.logger.debug("No XGB result, file none report")
            # This is meant to say that we really do not know.
            classifications.append(
                {
                    "classifierName": self.classifier_name + "SNGuess",
                    "classifierParams": self.classifier_version,
                    "classId": classId,
                    "probability": 1.0,
                }
            )
            return self.make_unit_result(compound, class_report)

        # Create first set of probabilities
        prob1 = is1
        prob21 = (1 - is1) * is21
        prob22 = (1 - is1) * (1 - is21)

        # We can now finalize probabilities for Reucurring events (2)
        classifications.append(
            {
                "classifierName": self.classifier_name + "SNGuess",
                "classifierParams": self.classifier_version,
                "classId": 21,
                "probability": prob21,
            }
        )
        classifications.append(
            {
                "classifierName": self.classifier_name + "SNGuess",
                "classifierParams": self.classifier_version,
                "classId": 22,
                "probability": prob22,
            }
        )
        classifications.append(
            {
                "classifierName": self.classifier_name + "SNGuess",
                "classifierParams": self.classifier_version,
                "classId": 1,
                "probability": prob1,
            }
        )

        # If Parsnip did not run we are done here
        if parsnip_class is None:
            self.logger.debug("No Parsnip result, file simple report")
            return self.make_unit_result(compound, class_report)

        # Create a new series of classifications including base Parsnip
        parsnip_classifications = [
            dict(d, classifierName=self.classifier_name + "SNGuess" + "Parsnip")
            for d in classifications
            if d["classId"] != 1
        ]
        for klass, parsnip_prob in parsnip_class.items():
            parsnip_classifications.append(
                {
                    "classifierName": self.classifier_name + "SNGuess" + "Parsnip",
                    "classifierParams": self.classifier_version,
                    "classId": klass,
                    "probability": parsnip_prob * prob1,
                }
            )

        # Scale parsnip with redshift prior if requested
        if self.use_priors and parsnip_z is not None:
            # Modify the parsnip probabilities based on the priors

            # Need to identify the host galaxy_color
            host_ug, mwebv = None, 0
            for dp in datapoints:
                if "mwebv" in dp["body"]:
                    host_ug = self.get_hostcol(dp["body"], parsnip_z)
                    mwebv = dp["body"]["mwebv"]
                    break

            # This requires a redshift to be present, but should also be the
            # case if parsnip has completed (?)
            parsnip_class = self.add_hostprior(
                self.add_rateprior(self.add_zprior(parsnip_class, parsnip_z)),
                host_ug,
                mwebv,
            )
            # Create a new series of classifications including priored  Parsnip
            pprior_classifications = [
                dict(
                    d,
                    classifierName=self.classifier_name
                    + "SNGuess"
                    + "Parsnip"
                    + "Prior",
                )
                for d in classifications
                if d["classId"] != 1
            ]
            for klass, parsnip_prob in (parsnip_class or {}).items():
                pprior_classifications.append(
                    {
                        "classifierName": self.classifier_name
                        + "SNGuess"
                        + "Parsnip"
                        + "Prior",
                        "classifierParams": self.classifier_version,
                        "classId": klass,
                        "probability": parsnip_prob * prob1,
                    }
                )

            # Done with prob collection, decide what to submit
            if self.multiple_classifiers:
                classifications.extend(parsnip_classifications)
                classifications.extend(pprior_classifications)
            else:
                # Only submit prior version here
                class_report["classifications"] = pprior_classifications
        # Done with prob collection, decide what to submit
        elif self.multiple_classifiers:
            classifications.extend(parsnip_classifications)
        else:
            # Only submit parsnip
            class_report["classifications"] = parsnip_classifications

        # Return report
        return self.make_unit_result(compound, class_report)
