#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2RunParsnip.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 24.09.2021
# Last Modified Date: 24.09.2021
# Last Modified By  : jnordin@physik.hu-berlin.de

from typing import Sequence, Literal, Any, Optional, Tuple, Union
import errno, os, re, backoff, copy, sys
import math
import gc
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsTiedLightCurveT2Unit import AbsTiedLightCurveT2Unit
from ampel.view.T2DocView import T2DocView
from ampel.view.LightCurve import LightCurve
from ampel.model.StateT2Dependency import StateT2Dependency
from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper


#import sncosmo
import parsnip
import lcdata



# All parsnip predictions that are not floats
dcast_pred = {
    'object_id' : str,
    'type' : str,
    'count' : int,
    'count_s2n_3' : int,
    'count_s2n_5' : int,
    'count_s2n_3_pre' : int,
    'count_s2n_3_rise' : int,
    'count_s2n_3_post' : int,
    'model_dof' : int,
}

dcast_class = {
    'object_id' : str,
}


class T2RunParsnip(AbsTiedLightCurveT2Unit):
    """
    Gathers information and runs the parsnip model and classifier.
    - Obtain model (read from file unless not in sncosmo registry)
    - Parse previous (chained) T2results for redshift or phase limits.
    - Creates suitable photometry, using the converter provided and  phase limits.
    - Defines model appropritately, including fit ranges and fixed parameters.
    - Run fit, potentially iterative in case of non-convergence.
    - Plot output if requested

    TODO:
    - Add option for redoing fits with disturbed initial conditions to avoid local minima
    - Add plot option
    - Add option for masking data?
    - Add MW correction
    """
    # Name (in case standard) or path to parsnip model to load
    parsnip_model: str
    # Path to classifier to apply to lightcurve fit. If not set, no classification will be done.
    parsnip_classifier: Optional[str]

    # Redshift usage options. Current options
    # T2MatchBTS : Use the redshift published by BTS and  synced by that T2.
    # T2DigestRedshifts : Use the best redshift as parsed by DigestRedshift.
    # None : run sncosmo template fit with redshift as free parameter OR use backup_z if set
    redshift_kind: Optional[str]
    # If loading redshifts from DigestRedshifts, provide the max ampel z group to make use of.
    # (note that filtering based on this can also be done for a potential t3)
    max_ampelz_group: int = 3
    # It is also possible to use fixed redshift whenever a dynamic redshift kind is not possible
    backup_z: Optional[float]
    # Finally, the provided lens redshift might be multiplied with a scale
    # Useful for lensing studies, or when trying multiple values
    scale_z: Optional[float]

    # Further fit parameters
    # Bounds - not yet implemented
    # sncosmo_bounds : Dict[ str, List[float] ] = {}
    # Remove MW dust absorption using SFD maps.
	# assumes that the position can be retrieved from the light_curve and
    # that the SFD_DIR env var is set to allow them to be found.
	# The default value of Rv will be used.
    # apply_mwcorrection : bool = False

    # Phase range usage. Current option:
    # T2PhaseLimit : use the jdmin jdmax provided in this unit output
    # None : use full datapoint range
    # (T2BayesianBlocks should be added)
    phaseselect_kind: Optional[str]

    # Save / plot parameters
    plot_suffix: Optional[str]
    plot_dir: Optional[str]


    # Which units should this be changed to
    t2_dependency: Sequence[StateT2Dependency[Literal[
                            "T2DigestRedshifts", "T2MatchBTS", "T2PhaseLimit"]]]







    def post_init(self)-> None:
        """
        Retrieve models.
        """

        # Load model and classifier
        self.model = parsnip.load_model(self.parsnip_model)
        self.classifier = None
        if self.parsnip_classifier:
            self.classifier = parsnip.Classifier.load(self.parsnip_classifier)



    def _get_redshift(self, t2_views) -> Tuple[Optional[float], Optional[str]]:
        """
        Can potentially also be replaced with some sort of T2DigestRershift tabulator?

        Assuming that only one instance of redshift sources exist
        """

        # Examine T2s for eventual information
        z: Optional[float] = None
        z_source: Optional[str] = None


        if self.redshift_kind in ['T2MatchBTS', 'T2DigestRedshifts']:
            for t2_view in t2_views:
                if not t2_view.unit == self.redshift_kind:
                    continue
                self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                # Parse this
                if self.redshift_kind == 'T2MatchBTS':
                    if 'bts_redshift' in t2_res.keys() and not t2_res['bts_redshift'] == '-':
                        z = float(t2_res['bts_redshift'])
                        z_source = "BTS"
                elif self.redshift_kind == 'T2DigestRedshifts':
                    if ('ampel_z' in t2_res.keys() and t2_res['ampel_z'] is not None
                            and t2_res['group_z_nbr'] <= self.max_ampelz_group):
                        z = float(t2_res['ampel_z'])
                        z_source = "AMPELz_group" + str(t2_res['group_z_nbr'])
        else:
            # Check if there is a fixed z set for this run, otherwise keep as free parameter
            if self.backup_z:
                z = self.backup_z
                z_source = "Fixed"
            else:
                z = None
                z_source = "Fitted"

        if z and self.scale_z:
            z *= self.scale_z
            z_source += " + scaled {}".format(self.scale_z)

        return z, z_source




    def _get_phaselimit(self, t2_views) -> Tuple[Optional[float], Optional[float]]:
        """
        Can potentially also be replaced with some sort of tabulator?

        """

        # Examine T2s for eventual information
        jdstart: Optional[float] = None
        jdend: Optional[float] = None

        if self.phaseselect_kind is None:
            jdstart = -np.inf
            jdend = np.inf
        else:

            for t2_view in t2_views:
                # So far only knows how to parse phases from T2PhaseLimit
                if not t2_view.unit == 'T2PhaseLimit':
                    continue
                self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                jdstart = t2_res['t_start']
                jdend = t2_res['t_end']

        return jdstart, jdend


    def _get_sncosmo_table(self, light_curve, jdstart, jdend) -> Table:
        """
        Generate sncosmo-like table from the provided lightcurve.

        This step is to be replaced by a tabulator.
        """

        # Datapoint selection filter
        dp_filter = [
            {'attribute': 'jd', 'operator': '>=', 'value': jdstart},
            {'attribute': 'jd', 'operator': '<=', 'value': jdend}
        ]
        phot = np.asarray(light_curve.get_ntuples(('jd', 'magpsf', 'sigmapsf',
		                                           'fid'), filters=dp_filter))
        phot_tab = Table(phot, names=('jd', 'magpsf', 'sigmapsf', 'fid'))
        phot_tab['band'] = 'ztfband'
        for fid, fname in zip([1, 2, 3], ['ztfg', 'ztfr', 'ztfi']):
            phot_tab['band'][phot_tab['fid'] == fid] = fname
        phot_tab['flux'] = 10 ** (-(phot_tab['magpsf'] - 25) / 2.5)
        phot_tab['fluxerr'] = np.abs(phot_tab['flux'] * (-phot_tab['sigmapsf'] / 2.5 * np.log(10)))
        phot_tab['zp'] = 25
        phot_tab['zpsys'] = 'ab'
        phot_tab.sort('jd')

        return phot_tab



    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(self,
                light_curve: LightCurve, t2_views: Sequence[T2DocView]
                ) -> Union[UBson, UnitResult]:
        """

        Fit the parameters of the initiated snocmo_model to the light_curve
        provided. Depending on the configuration, the provided T2DovViews
        are used to look for redshift information and any phase (time)
        limits for the fit.


        Parameters
        -----------
        light_curve: "ampel.view.LightCurve" instance.
        See the LightCurve docstring for more info.

        t2_records: List of T2Records from the following units (if available)
        T2DigestRedshifts (redshift parsed from catalogs)
        T2MatchBTS (redshifts synced from BTS page)
        T2PhaseLimit (fit time-limits as determined from lightcurve)

        Returns
        -------
        dict
        """

        # Initialize output dict
        t2_output: dict[str, UBson] = {"model" : self.parsnip_model,
		                              "classifier" : self.parsnip_classifier}

        # Obtain redshift
        z, z_source = self._get_redshift(t2_views)
        t2_output['z'] = z
        t2_output['z_source'] = z_source
        # A source class of None indicates that a redshift source was required, but not found.
        if t2_output['z_source'] is None:
            return t2_output

        # Check for phase limits
        (jdstart, jdend) = self._get_phaselimit(t2_views)
        t2_output['jdstart'] = jdstart
        t2_output['jdend'] = jdend
        if t2_output['jdstart'] is None:
            return t2_output

        # Obtain photometric table
        sncosmo_table = self._get_sncosmo_table(light_curve,
                                                t2_output['jdstart'], t2_output['jdend'])
        self.logger.debug('Sncosmo table {}'.format(sncosmo_table))

        # Fitting section

        # Deal with redshift and fit
        if z is None:
            sys.exit("Implement redshift fitting from parsnip_ztf_classification_no_redshift.ipynb")

        sncosmo_table.meta["redshift"] = z
        # Not sure why we need to define a dataset
        dataset = lcdata.from_light_curves([sncosmo_table])
        lc_predictions = self.model.predict_dataset(dataset)
        lc_classifications = self.classifier.classify(lc_predictions)
        foo = dict(lc_predictions[0][lc_predictions.colnames[1:]])
        t2_output["prediction"] = {k: dcast_pred[k](v)
		                           if k in dcast_pred and v is not None
								   else float(v) for k, v in foo.items()}
        foo = dict(lc_classifications[0][lc_classifications.colnames[1:]])
        t2_output["classification"] = {k: dcast_class[k](v)
		                               if k in dcast_class and v is not None
									   else float(v) for k, v in foo.items()}

        # Plot
        if self.plot_suffix and self.plot_dir:
            tname = ZTFIdMapper.to_ext_id(light_curve.stock_id)

            fig = plt.figure()
            ax = plt.gca()

            parsnip.plot_light_curve(dataset.light_curves[0], self.model, ax=ax)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 't2parsnip_%s.%s'%(tname, self.plot_suffix)))

            plt.close('fig')
            plt.close('all')
            del(fig)
            gc.collect()




        return t2_output
