#!/usr/bin/env python
# File:                Ampel-HU-astro/ampel/contrib/hu/util/LasairAnnotator.py
# License:             BSD-3-Clause
# Author:              jno
# Date:                9.10.2025
# Last Modified Date:  7.12.2025
# Last Modified By:    jno

from typing import Literal

import lasair  # type: ignore[import]

from ampel.secret.NamedSecret import NamedSecret


class LasairAnnotator:
    """
    Send an annotation to the Lasair system.
    """

    # Lasair annotator name, used for searching their.
    # It is assumed that each topic configures rules for explanation, url, report formating and classifiaction selection.
    lasair_topic: Literal["AMPEL"]
    # Which lasair version to use
    lasair_version: Literal["ztf", "lsst"]
    # API authentication for group owner
    lasair_api_token: NamedSecret[str]
    # Check existence
    check_existence: bool = True

    def post_init(self) -> None:
        self.endpoint = f"https://lasair-{self.lasair_version}.lsst.ac.uk/api"
        self.lasairclient = lasair.lasair_client(
            self.lasair_api_token.value, endpoint=self.endpoint
        )

    def annotate(
        self,
        objectId: str,
        classification: str,
        classdict: dict,
        explanation: str,
        version: str,
    ) -> bool:
        """
        Submit annotation to Lasair.:
        objectId: str - The Lasair object id to annotate (e.g. ZTF ID)
        classification: str - Single classification
        classdict: dict - Full class probability
        explanation: str - "Natural language" explanation
        """

        if self.check_existence:
            try:
                self.lasairclient.object(objectId, lite=True)
            except Exception:
                # Object not find in Lasair, return.
                return False

        lasairout = self.lasairclient.annotate(
            self.lasair_topic,
            objectId,
            classification,
            version=version,
            explanation=explanation,
            classdict=classdict,
            url="",  # Could be added later
        )

        return lasairout.get("status") == "success"
