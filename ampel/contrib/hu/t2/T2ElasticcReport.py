#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2ElasticcReport.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                25.08.2022
# Last Modified Date:  25.08.2022
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Literal, Union
from collections.abc import Sequence

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.content.T1Document import T1Document
from ampel.content.DataPoint import DataPoint
from ampel.view.T2DocView import T2DocView

# Tying classification output classes with ELAsTICC taxonomy classes
parsnip_taxonomy = {
    # Correct Elasticc outputs but they changed?
    # https://github.com/plasticc/taxonomy/blob/main/taxonomy.ipynb ?
    'SLSN':  131,
    'SNII':  113,
    'SNIa':  111,
    'SNibc': 112,
    'SNIbc': 112,
    'TDE':  132,
    'CART': 134,
    'ILOT': 133,
    'Mdwarf-flare': 122,
    'PISN': 135,
    'KN': 121,
    'SLSN-I': 131,
    'SNIa91bg': 115,
    'SNIax': 114,
    'dwarf-nova': 123,
    # mssing
    'uLens': 124,
    }

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
    tree_1v2: str = 'xgb_v6_simmod_tree12_'
    tree_21v22: str = 'xgb_v6_simmod_tree2122_'
    tree_1113v12: str = 'xgb_v6_simmod_tree121113_'

    # Setting for report to construct
    broker_name: str = 'AMPEL'
    broker_version: str = 'v0.2'
    classifier_name: str = 'SNguessParsnip'
    classifier_version: str = 'XGBUnified+Parsnip04'


    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal["T2RunParsnip", "T2XgbClassifier"]]]

    def submit(self, report: dict) -> str:
        """
        Placeholder for actually doing a quick T2 submit.
        """

        # For now we do a quick probability eval
        psum = 0
        for klass in report['classifications']:
            psum += klass['probability']

        if not psum==1:
            print('WARNING; SOMETHING NOT RIGHT.')

        return 'Not submitted.'


    def process(self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView]
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
        class_report = {
            'alertId': 'alertId',
            'diaSourceId': 'diaSourceId',
            'elasticcPublishTimestamp': 'elasticcPublishTimestamp',
            'brokerIngestTimestamp': 'brokerIngestTimestamp',
            'brokerName': self.broker_name,
            'brokerVersion': self.broker_version,
            'classifications': []
        }

        # Get alert info from T1Document if present
        for metarecord in compound['meta']:
            if 'alert' in metarecord.keys():
                class_report['alertId'] = metarecord['alert']
                class_report['brokerIngestTimestamp'] = metarecord['ts']
            if 'alert_ts' in metarecord.keys():
                class_report['elasticcPublishTimestamp'] = metarecord['alert_ts']

        # Get diaSourceId
        # Get the last diaSource id in datapoints
        for dp in reversed(datapoints):
            if 'diaSourceId' in dp['body'].keys():
                class_report['diaSourceId'] = dp['body']['diaSourceId']
                break

        # 1. Obtain base classifications from XGBClassifier and Parsnip
        is1, is21, is1113, parsnip_class = None, None, None, None

        # Parse t2views - should not be more than one.
        for t2_view in t2_views:
            self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
            # So far only knows how to parse phases from T2TabulatorRiseDecline
            if t2_view.unit == 'T2XgbClassifier':
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                if 'prob0' in t2_res.keys():
                    if t2_res['model']==self.tree_1v2:
                        is1 = t2_res['prob0']
                    elif t2_res['model']==self.tree_21v22:
                        is21 = t2_res['prob0']
                    elif t2_res['model']==self.tree_1113v12:
                        is1113 = t2_res['prob0']
            elif t2_view.unit == 'T2RunParsnip':
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                if 'classification' in t2_res.keys():
                    parsnip_class = {
                        parsnip_taxonomy[parsnip_type]: prob
                        for parsnip_type, prob in t2_res['classification'].items()
                    }

        # Did XGBClassifier run?
        if is1 is None or is21 is None or is1113 is None:
            self.logger.debug('No XGB result, file none report')
            # This is meant to say that we really do not know.
            class_report['classifications'].append(
                {
                    'classifierName': self.classifier_name,
                    'classifierParams': self.classifier_version,
                    'classId': 0,
                    'probability': 1.,
                }
            )
            class_report['t2_submit'] = self.submit(class_report)
            return class_report


        # Create first set of probabilities
        prob1 = is1
        prob21 = (1-is1) * is21
        prob22 = (1-is1) * (1-is21)

        # We can now finalize probabilities for Reucurring events (2)
        class_report['classifications'].append(
            {
                'classifierName': self.classifier_name,
                'classifierParams': self.classifier_version,
                'classId': 21,
                'probability': prob21,
            }
        )
        class_report['classifications'].append(
            {
                'classifierName': self.classifier_name,
                'classifierParams': self.classifier_version,
                'classId': 22,
                'probability': prob22,
            }
        )

        # Did Parsnip run?
        if parsnip_class is None:
            self.logger.debug('No Parsnip result, file simple report')
            class_report['classifications'].append(
                {
                    'classifierName': self.classifier_name,
                    'classifierParams': self.classifier_version,
                    'classId': 1,
                    'probability': prob1,
                }
            )
            return {
                "report": class_report,
                "t2_submit": self.submit(class_report),
            }


        # Time to integrate Parsnip results.
        # The operative question is whether we should make use of the
        # is1113 probability.
        # Something like
        # parsnip_prob12: float = # sum up prob in this subset
        # parsnip_prob1113: float = 1.-parsnip_prob12# sum up prob for the 11113 classes
        # prob111 = is1 * is1113 * parsnip[111] / parsnip_prob1113
        # ...
        # prob121 = is1 * (1-is1113) * parsnip[121] / parsnip_prob1113

        # However, lets first complete the base version
        for klass, parsnip_prob in parsnip_class.items():
            class_report['classifications'].append(
                {
                    'classifierName': self.classifier_name,
                    'classifierParams': self.classifier_version,
                    'classId': klass,
                    'probability': parsnip_prob*prob1,
                }
            )
        return {
            "report": class_report,
            "t2_submit": self.submit(class_report),
        }
