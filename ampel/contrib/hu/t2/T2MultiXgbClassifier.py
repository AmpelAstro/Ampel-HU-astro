#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2MultiXgbClassifier.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                03.07.2022
# Last Modified Date:  12.12.2022
# Last Modified By:    jnordin@physik.hu-berlin.de

from typing import Literal, Union, Dict
from collections.abc import Sequence


import numpy as np
import xgboost as xgb
import os

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsTiedStateT2Unit import AbsTiedStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.content.T1Document import T1Document
from ampel.content.DataPoint import DataPoint
from ampel.view.T2DocView import T2DocView
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import T2TabulatorRiseDeclineBase

class T2MultiXgbClassifier(AbsTiedStateT2Unit, AbsTabulatedT2Unit, T2TabulatorRiseDeclineBase):
    """
    For a range of xgboost classifier models, find a classification.
    This unit differs from T2XgbClassifier in that this unit assumes
    that a classifier applies to any number of detections and that every
    unit in the input dict is run.

    Will test whether the first (null) model of each classifier is
    more likely, P(model 0)>0.5. What this means is determined by the training.

    E.g. For the elasticc1v2 model, the null model corresponds to Elasticc class 1
    (and model 1 to class 2.)

    """


    model_folder: str = '/home/jnordin/github/zweilasticc/ml_workshop/elasticc/iter100'
    # Dict with names and path (in model_folder) for each classifier
    model_dict: dict[str, str] = {}

    _classifiers: dict[str,xgb.XGBClassifier]

    # Columns (in order) used for training
    # Which columsn to use for training. But these seem not to align with
    # parquet file!!!!
    use_cols: list[str] = ['ndet', 'jd_det', 'jd_last', 't_predetect', 'mag_det',
        'band_det_id', 'mag_last', 'band_last_id', 't_lc', 'bool_pure',
        'jd_peak_lsstz', 'bool_rise', 'bool_fall', 'bool_peaked', 'bool_fastrise',
        'bool_fastfall', 'bool_hasgaps', 'det_bands', 'last_bands',
        'jd_peak_lsstr', 'rise_slope_lsstz', 'rise_slopesig_lsstz',
        'fall_slope_lsstz', 'fall_slopesig_lsstz', 'jd_peak', 't_rise', 't_fall',
        'peak_bands', 'jd_peak_lssty', 'lsstz-lssty_last', 'lsstz-lssty_peak',
        'fall_slope_lssty', 'fall_slopesig_lssty', 'jd_peak_lssti',
        'fall_slope_lssti', 'fall_slopesig_lssti', 'jd_peak_lsstu',
        'fall_slope_lsstr', 'fall_slopesig_lsstr', 'jd_peak_lsstg',
        'rise_slope_lsstg', 'rise_slopesig_lsstg', 'rise_slope_lsstu',
        'rise_slopesig_lsstu', 'rise_slope_lssti', 'rise_slopesig_lssti',
        'lsstu-lsstg_last', 'fall_slope_lsstu', 'fall_slopesig_lsstu',
        'lssti-lsstz_last', 'rise_slope_lssty', 'rise_slopesig_lssty',
        'fall_slope_lsstg', 'fall_slopesig_lsstg', 'lsstg-lsstr_last',
        'lsstr-lssti_last', 'lsstu-lsstg_peak', 'rise_slope_lsstr',
        'rise_slopesig_lsstr', 'lsstg-lsstr_peak', 'lsstz-lssty_det',
        'lssti-lsstz_peak', 'lsstr-lssti_peak', 'lssti-lsstz_det',
        'lsstg-lsstr_det', 'lsstr-lssti_det', 'lsstu-lsstg_det',
        'z', 'z_err', 'host_sep']

    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal["T2TabulatorRiseDecline", "T2ElasticcRedshiftSampler"]]]


    def post_init(self) -> None:
        """
        Load models
        """

        self._classifiers = {}

        for modelid, modelpath in self.model_dict.items():
            fname = os.path.join(self.model_folder, modelpath )
            model = xgb.XGBClassifier()
            model.load_model(fname=fname)
            self._classifiers[modelid] =  model



    def process(self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView]
    ) -> Union[UBson, UnitResult]:
        """

        Extract results from TabulatorRiseDecline and apply suitable model.

        Parameters
        -----------

        t2_records: List of T2Records from (here) TabulatorRiseDecline.

        Returns
        -------
        dict
        """


        t2data = {}
        # Parse t2views - should not be more than one.
        for t2_view in t2_views:
            self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
            # So far only knows how to parse phases from T2TabulatorRiseDecline
            if t2_view.unit == 'T2TabulatorRiseDecline':
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                t2data.update(t2_res)
            if t2_view.unit == 'T2ElasticcRedshiftSampler':
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                # For some reason we trained xgb using z, zerr and host_sep
                zdata = {'z':None, 'z_err':None, 'host_sep': None}
                if t2_res['z_source'] in ['HOSTGAL2_ZQUANT', 'HOSTGAL_ZQUANT', 'HOSTGAL_ZSPEC', 'default']:
                    if len(t2_res['z_samples'])==3:
                        # This was the sampling used for job_training
                        zdata['z'] = t2_res['z_samples'][1]
                        zdata['z_err'] = t2_res['z_samples'][1] - t2_res['z_samples'][0]
                        zdata['host_sep'] = t2_res['host_sep']
                    elif len(t2_res['z_samples'])==4 and t2_res['z_samples'][0]==0.01:
                        # ... unless the event was hostless, in which case
                        # xgboost use None values.
                        pass
                    else:
                        self.logger.info('Do not know how to handle z info', extra=t2_res)
                        return {'model': None, 'xgbsuccess': False}
                t2data.update(zdata)


        # If we did not find chained results from TabulatorRiseDecline,
        # then run the feature extraction
        if not 'ndet' in t2data.keys():
            flux_table = self.get_flux_table(datapoints)

            # Cut the flux table if requested
            if self.max_ndet>0 and len(flux_table)>self.max_ndet:
                flux_table = self.cut_flux_table(flux_table)

            # Calculate features
            features = self.compute_stats(flux_table)
            t2data.update( features )




        # No detections - none of the models will work
        if t2data['ndet']==0:
            # Explore whether negative detections can be used to parse type
            if 'nnegdet' in t2data.keys() and t2data['nnegdet']>0:
                if t2data['z'] is not None and t2data['z']>0.001:
                    # High redshifts are typically AGNs
                    return {'model': 'directEval', 'xgbsuccess': False, 'cause': 'No pos det', 'direct_eval':'AGN'}
                else:
                    # Low redshifts are (typically) either uLens or EB.
                    # The former have a large fraction of neg det (if any)
                    nfrac = float(t2data['nnegdet']) / t2data['alldet']
                    if nfrac > 0.2:
                        return {'model': 'directEval', 'xgbsuccess': False, 'cause': 'No pos det', 'direct_eval':'uLens'}
                    else:
                        return {'model': 'directEval', 'xgbsuccess': False, 'cause': 'No pos det', 'direct_eval':'EB'}
            # Otherwise, not much to go by
            return {'model': 'directEval', 'xgbsuccess': False, 'imodel':-1, 'cause': 'No sig. det'}
        if t2data['success']==False:
            return {'model': 'directEval', 'xgbsuccess': False, 'imodel':-1, 'cause': 'RiseDecline error.'}


        # Loop through and apply all models
        t2out: Dict[str, UBson] = {'model':'multiXgb', 'classifications':{} }
        for modelid, model in self._classifiers.items():
            # Run model
            # Can we find away round creating a numpy array?
            prob = model.predict_proba( np.array([ [t2data[col]
                                        if col in t2data.keys() else None
                                        for col in self.use_cols] ]) )
            t2out['classifications'][modelid] = {
                'prob0': float( prob[0][0] ),
                'is_0': (float(prob[0][0])>0.5)
            }


        return t2out
