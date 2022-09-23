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

# T2XgbClassifier can also yield direct evaluations for some cases
# when the main run fails. We here map these to elasticc codes
direct_evaluations = {
    'AGN': 22,
    'uLens': 124,
    'EB': 21
}


# Prior section - can be applied as demo, but not matching the simulations.

# Elasticc redshift priors.
# Warning: Based on normalized distributions of training sample, no rate info taken into account.
# Bins based on z 0.1 bins starting at -0.09 and extending to 3.4
# Probabilities assumed to add up to 1000 for each redshift bin.
zmap = {'CART': {0: 9, 1: 124, 2: 135, 3: 143, 4: 143, 5: 127, 6: 97, 7: 62, 8: 33, 9: 16, 10: 8, 11: 6, 12: 6, 13: 5, 14: 5, 15: 6, 16: 7, 17: 6, 18: 3, 19: 1, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'uLens': {0: 304, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SLSN-I': {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 11, 6: 20, 7: 35, 8: 55, 9: 80, 10: 112, 11: 153, 12: 215, 13: 307, 14: 421, 15: 536, 16: 645, 17: 746, 18: 827, 19: 879, 20: 911, 21: 936, 22: 958, 23: 976, 24: 988, 25: 995, 26: 998, 27: 999, 28: 999, 29: 999, 30: 999, 31: 999, 32: 999, 33: 999, 34: 999, 35: 999}, 'dwarf-nova': {0: 304, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SNIa91bg': {0: 5, 1: 78, 2: 100, 3: 127, 4: 151, 5: 156, 6: 132, 7: 88, 8: 45, 9: 17, 10: 6, 11: 3, 12: 2, 13: 2, 14: 1, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SNIa': {0: 1, 1: 19, 2: 26, 3: 39, 4: 60, 5: 92, 6: 134, 7: 181, 8: 217, 9: 225, 10: 200, 11: 150, 12: 93, 13: 48, 14: 23, 15: 11, 16: 5, 17: 2, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'ILOT': {0: 24, 1: 257, 2: 223, 3: 173, 4: 114, 5: 59, 6: 23, 7: 6, 8: 1, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SNIax': {0: 4, 1: 64, 2: 77, 3: 95, 4: 116, 5: 136, 6: 145, 7: 137, 8: 111, 9: 77, 10: 46, 11: 23, 12: 9, 13: 3, 14: 1, 15: 1, 16: 1, 17: 2, 18: 1, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'TDE': {0: 3, 1: 46, 2: 56, 3: 69, 4: 84, 5: 100, 6: 113, 7: 123, 8: 133, 9: 141, 10: 145, 11: 142, 12: 130, 13: 108, 14: 83, 15: 63, 16: 50, 17: 43, 18: 38, 19: 34, 20: 30, 21: 25, 22: 20, 23: 14, 24: 9, 25: 4, 26: 1, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'KN': {0: 29, 1: 293, 2: 234, 3: 166, 4: 98, 5: 45, 6: 15, 7: 3, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SNII': {0: 3, 1: 47, 2: 59, 3: 75, 4: 94, 5: 112, 6: 126, 7: 131, 8: 127, 9: 117, 10: 105, 11: 95, 12: 86, 13: 77, 14: 70, 15: 66, 16: 63, 17: 59, 18: 48, 19: 33, 20: 21, 21: 11, 22: 5, 23: 2, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'SNibc': {0: 4, 1: 60, 2: 73, 3: 89, 4: 105, 5: 119, 6: 127, 7: 128, 8: 121, 9: 107, 10: 90, 11: 72, 12: 53, 13: 36, 14: 25, 15: 18, 16: 13, 17: 8, 18: 3, 19: 1, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'Mdwarf-flare': {0: 304, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}, 'PISN': {0: 0, 1: 7, 2: 10, 3: 16, 4: 24, 5: 38, 6: 62, 7: 100, 8: 152, 9: 215, 10: 284, 11: 352, 12: 401, 13: 409, 14: 368, 15: 296, 16: 211, 17: 131, 18: 76, 19: 49, 20: 36, 21: 26, 22: 15, 23: 6, 24: 1, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0, 30: 0, 31: 0, 32: 0, 33: 0, 34: 0, 35: 0}}
# Parsnip type priors
# Based on BTS _observed_ rate distributions, i.e. no knowledge of elasticc
# or LSST simulation properties.
# Probabilities assumed to add up to 1000
btsmap = {'SNIa': 663, 'uLens': 10, 'SLSN': 14, 'dwarf-nova': 10, 'SNIa91bg': 10, 'ILOT': 10, 'SNIax': 10, 'TDE': 10, 'KN': 10, 'SNII': 168, 'SNIbc': 58, 'Mdwarf-flare': 10, 'PISN': 10, 'CART': 10}



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

    # Use a redshift prior from elasticc alerts and obs rate prior from BTS.
    use_priors: bool = False


    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal["T2RunParsnip", "T2XgbClassifier"]]]

    def submit(self, report: dict) -> str:
        """
        Placeholder for actually doing a quick T2 submit.
        """

        # Possibly check total probability
        #psum = 0
        #for klass in report['classifications']:
        #    psum += klass['probability']
        #if not abs(psum-1)<0.01:
        #    self.log.info('Probability does not add.')

        return 'Not submitted.'

    def add_zprior(self, parsnip_prob: dict, z: float):
        """
        Adjust probabilities based on redshift prior.
        """
        # Out of range
        if  z>3.409 or z<0:
            return parsnip_prob

        zbin = int( (z + 0.09) / 0.1 )

        scaled_prob = {}
        for model_name, zdist in zmap.items():
            scaled_prob[parsnip_taxonomy[model_name]] = parsnip_prob[parsnip_taxonomy[model_name]] * float(zdist[zbin])/1000

        # Reweight probabilities
        p = sum( [v for v in scaled_prob.values()] )

        return { k:v/p for k,v in scaled_prob.items() }


    def add_rateprior(self, parsnip_prob: dict):
        """
        Modify fit probabilities based on observed rates.
        """

        scaled_prob = {}
        for model_name, model_prob in btsmap.items():
            scaled_prob[parsnip_taxonomy[model_name]] = parsnip_prob[parsnip_taxonomy[model_name]] * float(model_prob)/1000

        # Reweight probabilities
        p = sum( [v for v in scaled_prob.values()] )

        return { k:v/p for k,v in scaled_prob.items() }


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
        class_report: dict[str,UBson] = {
            'brokerName': self.broker_name,
            'brokerVersion': self.broker_version,
        }

        # use an alias variable to inform mypy that classifications is always a list
        class_report['classifications'] = classifications = []

        # Get alert info from T1Document if present
        for metarecord in compound['meta']:
            if isinstance(alert_id := metarecord.get('alert'), int):
                class_report['alertId'] = alert_id
                class_report['brokerIngestTimestamp'] = metarecord['ts']
            if isinstance(alert_ts := metarecord.get('alert_ts'), int):
                class_report['elasticcPublishTimestamp'] = alert_ts

        # Get diaSourceId
        # Get the last diaSource id in datapoints
        for dp in reversed(datapoints):
            if 'diaSourceId' in dp['body'].keys():
                class_report['diaSourceId'] = dp['body']['diaSourceId']
                break

        # 1. Obtain base classifications from XGBClassifier and Parsnip
        is1, is21, is1113, parsnip_class, parsnip_z, direct_eval = None, None, None, None, None, None

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
                else:
                    # Direct evaluation available even though XGB did not run
                    direct_eval = t2_res.get('direct_eval', None)

            elif t2_view.unit == 'T2RunParsnip':
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res

                if 'classification' in t2_res.keys():
                    parsnip_class = {
                        parsnip_taxonomy[parsnip_type]: prob
                        for parsnip_type, prob in t2_res['classification'].items()
                    }
                parsnip_z = t2_res.get('z_at_minchi', None)

        # Did XGBClassifier run?
        if is1 is None or is21 is None or is1113 is None:
            classId = 0   # Default for non-reply

            # Replace with  direct evaluation from XGB if present
            if direct_eval is not None:
                classId = direct_evaluations[direct_eval]

            self.logger.debug('No XGB result, file none report')
            # This is meant to say that we really do not know.
            classifications.append(
                {
                    'classifierName': self.classifier_name,
                    'classifierParams': self.classifier_version,
                    'classId': classId,
                    'probability': 1.,
                }
            )
            return {
                "report": class_report,
                "t2_submit": self.submit(class_report),
            }

        # Create first set of probabilities
        prob1 = is1
        prob21 = (1-is1) * is21
        prob22 = (1-is1) * (1-is21)

        # We can now finalize probabilities for Reucurring events (2)
        classifications.append(
            {
                'classifierName': self.classifier_name,
                'classifierParams': self.classifier_version,
                'classId': 21,
                'probability': prob21,
            }
        )
        classifications.append(
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
            classifications.append(
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


        # Scale parsnip with redshift prior if requested
        if self.use_priors:
            # Redshift
            if parsnip_z is not None:
                parsnip_class = self.add_zprior(parsnip_class, parsnip_z)
            # Rate
            parsnip_class = self.add_rateprior(parsnip_class)

        # However, lets first complete the base version
        for klass, parsnip_prob in parsnip_class.items():
            classifications.append(
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
