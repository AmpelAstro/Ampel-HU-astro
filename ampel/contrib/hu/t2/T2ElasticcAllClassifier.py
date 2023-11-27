#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2ElasticcAllClassification.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                22.08.2023
# Last Modified Date:  23.11.2023
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Any, Iterable, Literal, TypedDict
import numpy as np
import copy


from ampel.contrib.hu.t2.T2Elasticc2Classifier import BaseElasticc2Classifier
from ampel.contrib.hu.t2.T2Elasticc2Classifier import add_bts_rateprior
from ampel.contrib.hu.t2.T2Elasticc2Classifier import add_elasticc2_galcolprior



class T2ElasticcAllClassifier(BaseElasticc2Classifier):
    """

    Feature extraction, submission in Base.

    Classification structure:
    - 1. Extract rise/decline features. (already done)
    - 2. Try to run parsnip
    - 3. If that works, run allInclusive Xgb model.
    - 4. If that did not work, run all Xgb model.
    - 5. Create base results: 3 or 4.
    - 6. Add prior result variant.
    - 7. Add ping result variant.
    - 8. If (4) worked, create add this as specific result.
    (should there be ping or prior variants of these?)

    """

    # Setting for report to construct
    classifier_name: str = 'ElasticcMonster'
    classifier_version: str = 'v231123'


    # Ping limit (SNIa)
    pinglim: float = 0.5
    pingtrange: list[float] = [15., 60.]
    pingmaxmag: float = 22.

    def classify(self, base_class_dict, features, flux_table) -> tuple[list[dict], dict]:
        """

        """

        # 2. Try to parsnip
        if features['ndet']<1:
            parsnip_class: dict[str, Any] = {'Failed':'few det'}
        else:
            parsnip_class = self.get_parsnip_class('snlong', features, flux_table)

        # 3./4.  Run Xgb classifier
        # After validation, we could add a clause to not even run the xgb model
        if 'Failed' in parsnip_class.keys():
            xgb_prob = self.get_multixgb_class('all', features)
        else:
            # Add parsnip features to features
            features.update( {'parsnip_'+featname : featval for featname, featval in  parsnip_class['prediction'].items() } )
            features['parsnip_modelchidof'] = features['parsnip_model_chisq'] / features['parsnip_model_dof']
            xgb_prob = self.get_multixgb_class('all_incl', features)

        # 5.  Create core result
        multixgb_class =  copy.deepcopy(base_class_dict)
        if features['ndet']==0:
            multixgb_class['classifications'].append( {'classId': 2000, 'probability': 1.0} )
        else:
            multixgb_class['classifications'].extend( [
                    {'classId': fitclass, 'probability': float(fitprop)}
                        for fitclass, fitprop in xgb_prob.items()
                ]
            )
        multixgb_class['classifierName'] += 'AllIncl'
        output_classes = [multixgb_class]
        prob = {'allxgb': {str(k):v for k,v in xgb_prob.items()}  }

        # 6. Create post version of core
        all_post_class = copy.deepcopy(multixgb_class)
        all_post_class['classifierName'] += 'Post'
        # Z prior not quite logical, and z is used by parsnip/xgb - skip?
        # all_post_class = add_elasticc2_zprior(all_post_class, features['z'])
        all_post_class = add_bts_rateprior(all_post_class)
        all_post_class = add_elasticc2_galcolprior(all_post_class,
                    features.get("hostgal_mag_u-g"), features.get("mwebv") )
        output_classes.append( all_post_class )

        # 7. Create PING version.
        # List historic artifact when several models were used to gather SNIa probabilities
        snprob = [ classp['probability'] for classp in multixgb_class['classifications'] if classp['classId']==2222 ]
        # One should be defined
        if len(snprob)>0 and max(snprob)>self.pinglim:
            # Check other features
            if (self.pingtrange[0] <= features.get('t_lc',0) <= self.pingtrange[1]
                and self.pingmaxmag>features.get('mag_last',99)):
                ping_class = copy.deepcopy(base_class_dict)
                ping_class['classifierName'] += 'Ping'
                ping_class['classifications'].extend( [
                        {'classId': 2222, 'probability': 1.0},
                        ])
                output_classes.append(ping_class)


        # 8. Report version only included if parsnip ran
        if 'prediction' in parsnip_class.keys():
            parsnip_class =  copy.deepcopy(multixgb_class)
            parsnip_class['classifierName'] += 'Ltd'
            output_classes.append(parsnip_class)


        # Verification
        for report in output_classes:
            testsum = sum([d['probability'] for d in report['classifications']])
            if np.abs(testsum-1)>0.01:
                self.logger.info('prob ne 1',extra=report)


        return (output_classes, prob)
