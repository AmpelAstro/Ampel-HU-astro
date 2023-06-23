#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2Elasticc2Classification.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                25.08.2022
# Last Modified Date:  13.12.2022
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Literal, Union
from collections.abc import Sequence
import numpy as np

import joblib

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.content.T1Document import T1Document
from ampel.content.DataPoint import DataPoint
from ampel.view.T2DocView import T2DocView

from ampel.contrib.hu.t2.T2ElasticcRedshiftSampler import get_elasticc_redshift_samples
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import T2TabulatorRiseDeclineBase, FitFailed, BaseLightCurveFeatures


class T2Elasticc2Classifier(AbsStateT2Unit, AbsTabulatedT2Unit, T2TabulatorRiseDeclineBase, BaseLightCurveFeatures):
    """
    Carry out operations which yield one (or more) Elasticc Taxonomy classifications.

    Base feature extraction structure:
    - 1a. Extract alert ID.
    - 1b. Extract host props (a la the RedshiftSampler)
    - 1c. Extract RiseDecline features.

    A series of base binary trees.
    - 2a. Run galactic vs nongalactic tree.
    - 2b. Run agn vs [kn+sn+long] IF significant nongalactic prob, otherwise set all of these to zero.
    - 2c. Run ulens vs [periodic+mdwarf+dwarf] IF significant galactic prob, otherwise set all of these to zero.
    - 2d. Run periodic star subclassification IF significant galactic prob, otherwise set all of these to zero.

    If there are few detection + short lightcurve:
    - 3a. Run kn vs [sn+long]
    - 3b. Run [mdwarf, dwarf-nova] x [periodic+mdwarf/dwarf-nova]
    [kn,mdwarf and dwarf-nova prob set to zero if the condition is not fulfilled.]

    If there are multiple detections:
    - 4a. Run parsnip for sn+long
    (- 4b. Periodic star subclassification could be here if we do not always run it.
    Note that the limits for sec 3 and 4 might not be the same. There are KN with up to 10 detections (even if very rare) and the parsnip model can work with less data (they use an internal algorithm).

    Adding priors:
    - 5. Accept a classification dict. Return this with rescaled probabilities based on some priors.
    # Lingering questions:

    Classifiers are expected to be one of:
    - xgboost binary classifier dumped jusing joblib, together with columns to use:
        {'columns':['ndet','t_lc',...],'model':...}
    - xgboost multi classifier...
    - parsnip
    Each are configured by label:path

    - The order of how classes are separate can be changed. For example could ulens events be picked out later, inside the low/high ndet forks. Same thing for AGN. Would this improve anything?
    - Should ulens be considered part of the periodic star multiclasser?
    - Should we also try to subclassify periodic stars with very few detections? Or just return the group id? A: At least for some types it works well with few det, so try.
    - What about the special "inclusive" KN channel?
    """


    # Setting for report to construct
    broker_name: str = 'AMPEL'
    broker_version: str = 'v0.5'
    classifier_name: str = 'ElasticcLive'
    classifier_version: str = 'v230622'

    ## Paths to classifiers to load.
    paths_xgbbinary: dict[str,str] = {
        'kn_vs_nonrec': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring'
    }
    paths_xgbmulti: dict[str,str] = {
    }
    paths_parsnip: dict[str,str] = {

    }


    # Parameters controlling host galaxy information
    # How many redshifts samples to return ? Implemented are 1,3 or 5.
    nbr_samples: int = 3
    # Check for final/spectroscopic redshifts, and use if present
    use_final_z: bool = False
    # Which bands should be used to estimate the galaxy color?
    use_galcol: list[list[str]] = [['u','i'], ['u','g']]

    # Parameters controlling feature extraction
    max_ndet: int = 200000


    def read_class_models(self) -> None:
        self.class_xgbbinary = {
            label : joblib.load(path)
                for label, path in self.paths_xgbbinary.items()
        }

    def post_init(self) -> None:
        self.init_lightcurve_extractor()
        self.read_class_models()

    def get_binary_class(self, classlabel, features):
        """
        Return the classification for the labeled xgb model
        Classify based on features ordered according to columns
        """

        return self.class_xgbbinary[classlabel]['model'].predict_proba(
            np.array( [ features.get(col, np.nan)
                    for col in self.class_xgbbinary[classlabel]['columns'] ] ).reshape(1,-1)
        )


    def submit(self, class_reports: list[dict]) -> list[dict]:
        """
        Submit reports if so configured, record response.
        """

        return [
                {"class_report": report, "submitted": False}
                for report in class_reports
            ]


    def classify(self, base_class_dict, features) -> list[dict]:
        """
        Based on the provided features, extend the base_class dict into
        one or more elasticc2 classification reports.

        Meant to be replaced in other classes.
        """

        testlabel = 'kn_vs_nonrec'
        gotclass = self.get_binary_class(testlabel, features)[0]

        my_class = {k:v for k,v in base_class_dict.items()}
        my_class['classifications'].extend( [
            {'classId': 2232, 'probability': float(gotclass[0])},
            {'classId': 2220, 'probability': float(gotclass[1])},
            ] )

        return [my_class]

    def process(self,
        compound: T1Document,
        datapoints: Sequence[DataPoint]
        ) -> Union[UBson, UnitResult]:
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
        class_report: dict[str,UBson] = {
            'brokerName': self.broker_name,
            'brokerVersion': self.broker_version,
        }
        # use an alias variable to inform mypy that classifications is always a list
        class_report['classifications'] = classifications = []


        ## 1a+b: extract alert information

        # Get alert info from T1Document if present
        for metarecord in compound['meta']:
            if isinstance(alert_id := metarecord.get('alert'), int):
                class_report['alertId'] = alert_id
                class_report['brokerIngestTimestamp'] = metarecord['ts']
            if isinstance(alert_ts := metarecord.get('alert_ts'), int):
                class_report['elasticcPublishTimestamp'] = alert_ts

        # Extract information from the diaObject
        for dp in reversed(datapoints):
            if 'diaSourceId' in dp['body'].keys():
                class_report['diaSourceId'] = dp['body']['diaSourceId']
                break
        # Get the last diaSource id in datapoints
        features = {}
        for dp in reversed(datapoints):
            if 'mwebv' in dp['body'].keys():
                features.update( get_elasticc_redshift_samples(dp['body'],
                        self.nbr_samples, self.use_galcol, self.use_final_z) )
                features['mwebv'] = dp['body']['mwebv']
                break


        ### 1c: Extract features
        # Convert input datapoints to standardized Astropy Table
        # Using standard tabulators
        flux_table = self.get_flux_table(datapoints)

        # Cut the flux table if requested
        if self.max_ndet>0 and len(flux_table)>self.max_ndet:
            flux_table = self.cut_flux_table(flux_table)

        try:
            # Calculate RiseDecline
            features.update( self.compute_stats(flux_table) )
        except FitFailed:
            # Continue, but likely cant say much without these
            pass
        try:
            # Calculate light_curve features
            features.update(self.extract_lightcurve_features(flux_table))
        except ValueError:
            self.logger.info('lightcurve extraction fail')
            pass
        # Create averaged values
        avgfeat = self.average_filtervalues(features)
        features.update(avgfeat)

        # Fix galaxy color mismatch: xgboost assumes galcol corresponds to
        # u-i color, while galcol here is a dict with potentially differnet
        # colors. Fix this.
        if features.get('galaxy_color') is not None:
            for col, val in features['galaxy_color'].items():
                features['galaxy_color_'+col] = val
            features['galaxy_color'] = features.get('galaxy_color_u-i',None)
        # Similarly we need homogeneous treatment of redshifts
        zkind = {'HOSTGAL_ZQUANT':1,'HOSTGAL2_ZQUANT':2,'HOSTGAL_ZSPEC':3,'HOSTGAL2_ZSPEC':3,'default':0}
        if features['z_source'] == 'default':
            features['z_source'] = 0
            features['z'] = 0
        else:
            # Assuming three samples - this is what was used in feature training
            features['z'] = features['z_samples'][1]
            features['z_source'] = zkind[features['z_source']]

        # TODO: Check for minimal data to attempt classifcation?
        # Otherwise, prepare and report classification: {'classId':300, probability: 1}

        # Run the classifiers, tie together to one or more reports
        classification_reports = self.classify(class_report, features)

        # Potentially, immediately submit reports through kafka
        replies = self.submit( classification_reports )

        # Prepare t2 document
        return {
            "reports_replies": replies,
            "features": features
        }
