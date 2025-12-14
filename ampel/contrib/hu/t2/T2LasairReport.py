#!/usr/bin/env python
# File:                Ampel-HU-astro/ampel/contrib/hu/t2/T2LasairReport.py
# License:             BSD-3-Clause
# Author:              jno
# Date:                7.12.2025
# Last Modified Date:  7.12.2025
# Last Modified By:    jno

from collections.abc import Sequence

from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2ClassificationReport import T2ClassificationReport
from ampel.contrib.hu.util.LasairAnnotator import LasairAnnotator
from ampel.struct.UnitResult import UnitResult
from ampel.view.T2DocView import T2DocView


class T2LasairReport(T2ClassificationReport, LasairAnnotator):
    """
    Propagate classification results to Lasair.
    """

    # Which classification results should go into Lasair annotations
    # For results from T2RunParsnipRiseDecline, which requires the classifier (eg ParsnipRiseDecline) and model (eg snlong)
    # to be specified. as lasair_classification_source = ['ParsnipRiseDecline', 'snlong']
    lasair_classification_source: None | Sequence[str] = [
        "ParsnipRiseDecline",
        "snlong",
    ]

    # Min prob for classification to be reported
    lasair_min_classification_prob: float = 0.5

    # Info propagated to lasair
    explanation: str = "AMPEL classification"

    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UnitResult:
        # First, run the parent classification report processing
        res = super().process(compound, datapoints, t2_views)
        t2output = {"T2ClassificationReport": res.body, "lasair_annotations": {}}
        if res.body is None:
            return UnitResult(
                body=t2output,
            )

        # For now we only have annotations based on T2RunParsnipRiseDecline results
        if (
            self.lasair_classification_source is not None
            and isinstance(res.body, dict)
            and res.body["classification"] is not None
        ):
            # Check that we have classification results to propagate
            for classif in res.body["classification"]:
                if classif["name"] != self.lasair_classification_source[0]:
                    continue
                for model in classif["models"]:
                    if model["model"] != self.lasair_classification_source[1]:
                        continue
                    # Finally can gather the data
                    classdict = model["probabilities"]
                    most_likely_class = max(classdict, key=classdict.get)
                    if (
                        classdict[most_likely_class]
                        > self.lasair_min_classification_prob
                    ):
                        classification = most_likely_class
                    else:
                        classification = "Unclear"

                    # Classifications in format P(class), so remove the P( and )
                    if (
                        classification is not None
                        and classification.startswith("P(")
                        and classification.endswith(")")
                    ):
                        classification = classification[2:-1]

                    # Then, propagate to Lasair
                    t2output["lasair_annotations"] = {
                        f"{self.lasair_classification_source[0]}_{self.lasair_classification_source[1]}": self.annotate(
                            res.body["object"]["external_id"],
                            classification,
                            classdict,
                            self.explanation,
                            classif["version"],
                        )
                    }

        return UnitResult(
            body=t2output,
        )
