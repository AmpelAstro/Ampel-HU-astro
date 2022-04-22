#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2RunParsnip.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 24.09.2021
# Last Modified Date: 06.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

from typing import Sequence, Literal, Any, Optional, Tuple, Union
import errno, os, re, backoff, copy, sys
import math
from scipy.stats import chi2
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


import parsnip
import lcdata

# The following three only used if correcting for MW dust
from sfdmap import SFDMap  # type: ignore[import]
import sncosmo             # type: ignore[import]
import extinction          # type: ignore[import]



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
    - Add option for masking data?
    """

    # Name (in case standard) or path to parsnip model to load
    parsnip_model: str
    # Path to classifier to apply to lightcurve fit. If not set, no classification will be done.
    parsnip_classifier: Optional[str]

    # Redshift usage options. Current options
    # T2MatchBTS : Use the redshift published by BTS and  synced by that T2.
    # T2DigestRedshifts : Use the best redshift as parsed by DigestRedshift.
    # T2ElasticcRedshiftSampler: Use a list of redshifts and weights from the sampler.
    # None : run sncosmo template fit with redshift as free parameter OR use backup_z if set
    redshift_kind: Optional[str]
    # If loading redshifts from DigestRedshifts, provide the max ampel z group to make use of.
    # (note that filtering based on this can also be done for a potential t3)
    max_ampelz_group: int = 3
    # It is also possible to use fixed redshift whenever a dynamic redshift kind is not possible
    # This could be either a single value or a list
    fixed_z: Optional[Union[float, Sequence[float]]]
    # Finally, the provided lens redshift might be multiplied with a scale
    # Useful for lensing studies, or when trying multiple values
    scale_z: Optional[float]

    # Remove MW dust absorption.
    # MWEBV is either derived from SFD maps using the position from light_curve
    # (assuming the SFD_DIR env var is set) OR retrieved from stock (ELASTICHOW?)
    # The default value of Rv will be used.
    apply_mwcorrection: bool = False


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
        "T2ElasticcRedshiftSampler",
        "T2DigestRedshifts",
        "T2MatchBTS",
        "T2PhaseLimit"]]]




    def post_init(self)-> None:
        """
        Retrieve models.
        """

        # Load model and classifier
        self.model = parsnip.load_model(self.parsnip_model)
        self.classifier = None
        if self.parsnip_classifier:
            self.classifier = parsnip.Classifier.load(self.parsnip_classifier)

        if self.apply_mwcorrection:
            self.dustmap = SFDMap()



    def _get_redshift(self, t2_views) -> Tuple[Optional[Sequence[float]], Optional[str]]:
        """
        Can potentially also be replaced with some sort of T2DigestRershift tabulator?

        Assuming that only one instance of redshift sources exist
        """

        # Examine T2s for eventual information
        z: Optional[Sequence[float]] = None
        z_source: Optional[str] = None
        z_weights: Optional[Sequence[float]] = None


        if self.redshift_kind in ['T2MatchBTS', 'T2DigestRedshifts', 'T2ElasticcRedshiftSampler']:
            for t2_view in t2_views:
                if not t2_view.unit == self.redshift_kind:
                    continue
                self.logger.debug('Parsing t2 results from {}'.format(t2_view.unit))
                t2_res = res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                # Parse this
                if self.redshift_kind == 'T2MatchBTS':
                    if 'bts_redshift' in t2_res.keys() and not t2_res['bts_redshift'] == '-':
                        z = [float(t2_res['bts_redshift'])]
                        z_source = "BTS"
                elif self.redshift_kind == 'T2DigestRedshifts':
                    if ('ampel_z' in t2_res.keys() and t2_res['ampel_z'] is not None
                            and t2_res['group_z_nbr'] <= self.max_ampelz_group):
                        z = [float(t2_res['ampel_z'])]
                        z_source = "AMPELz_group" + str(t2_res['group_z_nbr'])
                elif self.redshift_kind == 'T2ElasticcRedshiftSampler':
                    z = t2_res['z_samples']
                    z_source = t2_res['z_source']
                    z_weights = t2_res['z_weights']
        else:
            # Check if there is a fixed z set for this run, otherwise keep as free parameter
            if self.fixed_z:
                if isinstance(self.fixed_z, list):
                    z = self.fixed_z
                else:
                    z = [self.fixed_z]
                z_source = "Fixed"
            else:
                z = None
                z_source = "Fitted"

        if z and self.scale_z:
            z = [onez*self.scale_z for onez in z]
            z_source += " + scaled {}".format(self.scale_z)


        return z, z_source, z_weights




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


    def _deredden_mw_extinction(self, ebv, phot_tab, rv=3.1) -> Table:
        """
        For an input photometric table, try to correct for mw extinction.
        """

        # Find effective wavelength for all filters in phot_tab
        filterlist = set( phot_tab['band'] )
        eff_wave = [sncosmo.get_bandpass(f).wave_eff for f in filterlist]

        # Determine flux correction (dereddening) factors
        flux_corr = 10**(0.4* extinction.ccm89( np.array(eff_wave), ebv*rv, rv) )

        # Assign this appropritately to Table
        phot_tab['flux_original'] = phot_tab['flux']
        phot_tab['fluxerr_original'] = phot_tab['fluxerr']
        for k, band in enumerate(filterlist):
            phot_tab['flux'][(phot_tab['band']==band)] *= flux_corr[k]
            phot_tab['fluxerr'][(phot_tab['band']==band)] *= flux_corr[k]

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

        # Potentially correct for dust absorption
        if self.apply_mwcorrection:
            # Get ebv from coordiantes.
            # Here there should be some option to read it from journal/stock etc
            mwebv = self.dustmap.ebv(*light_curve.get_pos(ret="mean"))
            t2_output['mwebv'] = mwebv
            sncosmo_table = self._deredden_mw_extinction(mwebv, sncosmo_table)


        ## Obtain redshift(s) from catalog fit or a RedshiftSample
        z, z_source, z_weights = self._get_redshift(t2_views)
        t2_output['z'] = z
        t2_output['z_source'] = z_source
        t2_output['z_weights'] = z_weights
        # A source class of None indicates that a redshift source was required, but not found.
        if t2_output['z_source'] is None:
            return t2_output


        # Fitting section

        # If redshift should be fitted, we start with getting samples
        if z_source == 'Fitted':
            sys.exit('Parsnip redshift fit depracated??, use list of fixed_z.')
            #z, z_probabilities = self.model.predict_redshift_distribution(
            #                    sncosmo_table, max_redshift=self.max_fit_z)

        # Create a list of lightcurves, each at a discrete redshift
        lcs = []
        for redshift in z:
            use_lc = sncosmo_table.copy()
            use_lc.meta['object_id']  = f'parsnip_z{redshift:4f}'
            use_lc.meta['redshift'] = redshift
            lcs.append(use_lc)
        lc_dataset = lcdata.from_light_curves(lcs)

        # Peform actual model predictions and classifications
        lc_predictions = self.model.predict_dataset(lc_dataset)
        lc_classifications = self.classifier.classify(lc_predictions)


        # Cast result for storage and look at relative probabilities
        t2_output["predictions"] = {}
        t2_output["classifications"] = {}
        for i, redshift in enumerate(z):
            foo = dict(lc_predictions[i][lc_predictions.colnames[1:]])
            t2_output["predictions"][str(redshift)] = {k: dcast_pred[k](v)
                                        if k in dcast_pred and v is not None
								        else float(v) for k, v in foo.items()}
            # Not sure whether the dof could change? Normalizing now
            t2_output["predictions"][str(redshift)]['chi2pdf'] = chi2.pdf(
                                        foo['model_chisq']/foo['model_dof'], 1)
            foo = dict(lc_classifications[i][lc_classifications.colnames[1:]])
            t2_output["classifications"][str(redshift)] = {k: dcast_class[k](v)
		                                  if k in dcast_class and v is not None
									      else float(v) for k, v in foo.items()}

        # Marginalize over the redshift
        # p(c|y) = Integral[p(c|z,y) p(z|y) dz]
        types = lc_classifications.colnames[1:]
        dtype = lc_classifications[types[0]].dtype
        probabilities = lc_classifications[types].as_array().view((dtype, len(types)))
        # Now we could normalize these z prob and normalize types over redshifts
        z_probabilities = np.array( [lcfit['chi2pdf']
                    for redshift, lcfit in t2_output["predictions"].items()] )
        # Take redshift probabilities into account, if available
        if z_weights is not None:
            z_probabilities *= z_weights
        integrated_probabilities = z_probabilities.dot(probabilities)
        integrated_probabilities /= np.sum(integrated_probabilities)
        t2_output["marginal_lc_classifications"] = dict(zip(types, integrated_probabilities))
        # Find the best redshifts
        t2_output["z_at_minchi"] = z[np.argmax(z_probabilities)]
        # Map these to singular value predictions/lc_classifications
        # (wastes DB space, but possible to filter based on)
        t2_output["prediction"] = t2_output["predictions"][str(t2_output["z_at_minchi"])]
        t2_output["classification"] = t2_output["classifications"][str(t2_output["z_at_minchi"])]

        # Plot
        if self.plot_suffix and self.plot_dir:
            tname = ZTFIdMapper.to_ext_id(light_curve.stock_id)

            fig = plt.figure()
            ax = plt.gca()

            # Set redshift to best value and plot this fit
            lc_dataset.light_curves[0].meta['redshift'] = t2_output["z_at_minchi"]

            parsnip.plot_light_curve(lc_dataset.light_curves[0], self.model, ax=ax)
            plt.tight_layout()
            plt.savefig(os.path.join(self.plot_dir, 't2parsnip_%s.%s'%(tname, self.plot_suffix)))

            plt.close('fig')
            plt.close('all')
            del(fig)
            gc.collect()




        return t2_output
