#!/usr/bin/env python
# File:                Ampel-HU-astro/ampel/contrib/hu/util/AmpelLasair.py
# License:             BSD-3-Clause
# Author:              jno
# Date:                9.10.2025
# Last Modified Date:  9.10.2025
# Last Modified By:    jno

from ampel.base.AmpelBaseModel import AmpelBaseModel

from typing import Literal
import json, sys, settings
import lasair


class LasairAnnotator():
    """
    Send an annotation to the Lasair system.
    """

    # Lasair annotator name, used for searching their. 
    # It is assumed that each topic configures rules for explanation, url, report formating and classifiaction selection. 
    lasair_topic: Literal[ "AMPEL" ]
    # Which lasair version to use
    lasair_version: Literal["ztf","lsst"]
    # API authentication for group owner 
    lasair_api_token: NamedSecret[str]
    # Check existence 
    check_existence: bool = True
    

    def post_init(self) -> None:
        self.endpoint = "https://lasair-{}.lsst.ac.uk/api".format(self.lasair_version)
        self.lasairclient = lasair.lasair_client(self.lasair_api_token, endpoint=self.endpoint)

    def annotate(self, objectId: str, transientdict: dict) -> bool:    

        # Shape info for Lasair, according to topic
        if self.lasair_topic == "AMPEL":
            (classification, version, explanation, classdict, url) = self.get_dummy_annotation_info(transientdict)
        else:
            raise ValueError("Unknown Lasair topic %s" % self.lasair_topic)

        if self.check_existence:
            try:
                objectInfo = self.lasairclient.object( objectId, lite=True )
            except Exception as e:
                # Object not find in Lasair, return.
                print("Object %s not found in Lasair, not annotating." % objectId) 
                return 0

        self.lasairclient.annotate(
            self.lasair_topic, 
            objectId, 
            classification,
            version=version, 
            explanation=explanation, 
            classdict=classdict,
            url=url
        )

        return 1

    def get_dummy_annotation_info(self, transientdict: dict) -> tuple[str, float, str, dict, str]:
        """
        Dummy example of how to extract annotation info from a transient dict.
        """

        classification = "unknown"
        version = 1
        explanation = "No explanation available"
        classdict = {}
        url = "http://example.com"

        if "class" in transientdict:
            classification = transientdict["class"]
            explanation = f"Classified as {classification}"
            classdict = {classification: 0.9}
            url = f"http://example.com/{transientdict['id']}"

        return (classification, version, explanation, classdict, url)