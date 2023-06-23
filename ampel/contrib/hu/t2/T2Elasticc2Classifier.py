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
from ampel.contrib.hu.t2.T2RunParsnip import run_parsnip_zsample
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import T2TabulatorRiseDeclineBase, FitFailed, BaseLightCurveFeatures

# Parsnip models - avoid importing these if running in the light mode?
import parsnip
import lcdata



class BaseElasticc2Classifier(AbsStateT2Unit, AbsTabulatedT2Unit, T2TabulatorRiseDeclineBase, BaseLightCurveFeatures):
    """
    Base class for carrying out operations which yield one (or more) Elasticc Taxonomy classifications.
    Classifications carried out by classifier method, which is assumed to be
    implemented in child.

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
    broker_version: str = 'v0.8'
    classifier_name: str
    classifier_version: str

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

    # Parameters controlling any parsnip fit
    # The parsnip model zeropoint will by default be set to that of
    # the input lightcurve. This can be adjusted by zp offset.
    # Offset can e.g. appear if training performed with wrong zeropoint...
    parsnip_zeropoint_offset: float = 0


    def read_class_models(self) -> None:
        self.class_xgbbinary = {
            label : joblib.load(path)
                for label, path in self.paths_xgbbinary.items()
        }
        self.class_xgbmulti = {
            label : joblib.load(path)
                for label, path in self.paths_xgbmulti.items()
        }
        self.class_parsnip = {
            label: {
                'classifier': parsnip.Classifier.load(paths['classifier']),
                'model':     parsnip.load_model(paths['model'], threads=1) }
                            for label, paths in self.paths_parsnip.items()
        }

    def post_init(self) -> None:
        self.init_lightcurve_extractor()
        self.read_class_models()

    def get_xgb_class(self, classlabel, features):
        """
        Return the classification for the labeled xgb model
        Classify based on features ordered according to columns
        """

        return self.class_xgbbinary[classlabel]['model'].predict_proba(
            np.array( [ features.get(col, np.nan)
                    for col in self.class_xgbbinary[classlabel]['columns'] ] ).reshape(1,-1)
        )


    def get_parsnip_class(self, classlabel, features, lctable)->dict | None:
        """
        Return the classification for the labeled parsnip model
        Classify based on features ordered according to columns
        """
        pdict = run_parsnip_zsample(lctable,
            {'label':classlabel, 'z':features['z_samples'], 'z_weights':features['z_weights']},
            self.class_parsnip[classlabel]['model'],
            self.class_parsnip[classlabel]['classifier'],
            self.parsnip_zeropoint_offset )
        if pdict.get('Failed',None) is not None:
            return pdict
        return pdict



    def submit(self, class_reports: list[dict]) -> list[dict]:
        """
        Submit reports if so configured, record response.
        """

        return [
                {"class_report": report, "submitted": False}
                for report in class_reports
            ]


    def classify(self, base_class_dict, features, lc_table) -> (list[dict], dict):
        """
        Based on the provided features, extend the base_class dict into
        one or more elasticc2 classification reports.

        Meant to be replaced in other classes.
        """

        testlabel = 'kn_vs_nonrec'
        gotclass = self.get_xgb_class(testlabel, features)[0]

        my_class = {k:v for k,v in base_class_dict.items()}
        my_class['classifications'].extend( [
            {'classId': 2232, 'probability': float(gotclass[0])},
            {'classId': 2220, 'probability': float(gotclass[1])},
            ] )

        return ([my_class], {'binary':[float(gotclass[0])]})

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
            'classifierName': self.classifier_name,
            'classifierParams': self.classifier_version,
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
        (classification_reports, probabilities) = self.classify(class_report, features, flux_table)

        # Potentially, immediately submit reports through kafka
        replies = self.submit( classification_reports )

        # Prepare t2 document
        return {
            "reports_replies": replies,
            "features": features,
            "probabilities": probabilities
        }



class T2Elasticc2Classifier(BaseElasticc2Classifier):
    """

    Feature extraction, submission in Base.

    Classification steps implemented here:
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


    Q: What should the min ndet be for running parsnip? Remember that these are not counted
    the same way. Need to run and compare...

    """

    # Setting for report to construct
    classifier_name: str = 'ElasticcMonster'
    classifier_version: str = 'v230622'

    # Classifier mode?
    # fast: do not run parsnip (or multixgb if that is timeconsuming)
    # extragalactic: run parsnip if extragalctic prob>1%
    classifier_mode: Literal["full", "fast", "extragalactic"] = "extragalactic"

    ## Paths to classifiers to load.
    paths_xgbbinary: dict[str,str] = {
        'gal_vs_nongal': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring',
        'agn_vs_knsnlong': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring',
        'kn_vs_snlong': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring',
        'mdwarf_vs_galactic': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring',
        'dflare_vs_galactic': '/home/jnordin/Downloads/model_kn_vs_other_non_recurring',
    }
    paths_xgbmulti: dict[str,str] = {
    }
    paths_parsnip: dict[str,dict[str,str]] = {
        'snlong': {
'model':'/home/jnordin/github/ampel83elasticc2/ml_workshop/parsnip/trained_models/model-elasticc2_run20_SNLONG.h5',
'classifier':'/home/jnordin/github/ampel83elasticc2/ml_workshop/parsnip/trained_models_classifiers/elasticc2_run20_SNLONG_pred_classifier.pkl'
        }
    }

    def classify(self, base_class_dict, features, flux_table) -> (list[dict], dict):
        """
        Base:
        - galactic vs nongalactic

        TODO: Could improve speed by not running leafs if node
        below prob threshold.
        """


        pgal = self.get_xgb_class('gal_vs_nongal', features)[0]
        pagn = self.get_xgb_class('agn_vs_knsnlong', features)[0]
        pkn = self.get_xgb_class('kn_vs_snlong', features)[0]
        pmd = self.get_xgb_class('mdwarf_vs_galactic', features)[0]
        pdf = self.get_xgb_class('dflare_vs_galactic', features)[0]

        if self.classifier_mode in ['full','extragalactic']:
            if features['ndet']<4:
                print('NOT ENOUGH NDET; NOT RUNNING PARSNIP')
                parsnip_class = {'fail':'few det'}
            elif self.classifier_mode=='extragalactic' and pgal[1]<0.01:
                print('LOOKS GALACTIC')
                parsnip_class = {'fail':'galactic'}
            else:
                # attempt to run parsnip
                parsnip_class = self.get_parsnip_class('snlong', features, flux_table)
                print('parsnip', parsnip_class)


        my_class = {k:v for k,v in base_class_dict.items()}
        # Prob some more elegant way using the tree structure to add these
        my_class['classifications'].extend( [
            {'classId': 2232, 'probability': float(pgal[1]*pagn[1]*pkn[0])},
            {'classId': 2220, 'probability': float(pgal[1]*pagn[1]*pkn[1])},
            {'classId': 2332, 'probability': float(pgal[1]*pagn[0])},
            {'classId': 2320, 'probability': float(pgal[0])},
            ] )

        prob = {'binary':[float(pgal[0]), float(pagn[0]), float(pkn[0]), float(pmd[0]), float(pdf[0]) ],
                'multiple': [],
                'parsnip': parsnip_class
        }

        return ([my_class], prob)
