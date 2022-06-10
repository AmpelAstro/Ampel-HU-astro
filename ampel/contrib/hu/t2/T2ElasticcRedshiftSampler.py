#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2ElasticcRedshiftSampler.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 12.04.2022
# Last Modified Date: 12.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

from typing import Sequence, Union
#import errno, os, re, backoff, copy, sys
#import math
#from scipy.stats import chi2
#import gc
#import numpy as np
#from astropy.table import Table
#import matplotlib.pyplot as plt


from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.content.DataPoint import DataPoint
#from ampel.content.StockDocument import StockDocument
from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
#from ampel.abstract.AbsStockT2Unit import AbsStockT2Unit

#from ampel.view.T2DocView import T2DocView
#from ampel.view.LightCurve import LightCurve
#from ampel.model.StateT2Dependency import StateT2Dependency
#from ampel.ztf.util.ZTFIdMapper import ZTFIdMapper


class T2ElasticcRedshiftSampler(AbsPointT2Unit):
    """
    Parse the elasticc diaSource host information and
    returns a list of redshifts and weights.

    Note 1:
    Elasticc alerts contain both quantiles and mean value + error.
    According to information late May (Kessler), uncertainties will likely
    be normal, meaning that there will be no additional information in the
    quantiles. But they _might_ change this.

    Note 2:
    Elasticc training lightcurves do not seem to consistently have any photo-z
    sort_values, but rather an encoding as "REDSHIFT_FINAL(_ERR)".
    This might change with the test stream, but there could also be hostless sne.

    Note 3:
    Values for two potential host galaxies are provided. Entries such as
    e.g. 'HOSTGAL2_SNSEP' allows to calculate the probability. Again,
    this information not regularly provided. We can calculate this here,
    and still provide a combined set of redshifts with probabilities weighted
    for the two possible goals.

    Host galaxy properties (mass, sfr) also provided. We could use this to also
    provide constraints on the potential types. This can be saved into the T2Doc
    and considered in the T3.


    Note 4:
    The default setting for normal error distribution will be to use
    [-1.5 sig, 0, 1.5 sig] with relative weights [0.19684199, 0.60631602, 0.19684199]

    """

    # How many redshifts to return (up to some maximum)
    # (prop not all values allowed)
    max_zs: int = 3

    # Default redshifts and weights - in particular used for hostless events
    # Values set semi-randomly now.
    default_zs: Sequence[float] = [0.2, 0.4, 0.6, 1]
    default_weights: Sequence[float] = [0.2, 0.3, 0.3, 0.2]

    # Check for final/spectroscopic redshifts, and use if present
    # Usually would be the case, but for testing we might not want to
    use_final_z: bool = True


    def get_hostprob(self, hostinfo: dict) -> (float, float):
        """
        Use the information from the alert info to guess relative prob
        for the two different host galaxies.
        """

        if hostinfo['HOSTGAL_PHOTOZ']>0 and hostinfo['HOSTGAL2_PHOTOZ']<0:
            # The easy case, only have first galaxy
            return (1.0, 0.0)
        elif hostinfo['HOSTGAL_PHOTOZ']<0 and hostinfo['HOSTGAL2_PHOTOZ']<0:
            # "Hostless" - if these exist
            self.logger.info('Hostless')
            return (0.0, 0.0)
        elif hostinfo['HOSTGAL_PHOTOZ']>0 and hostinfo['HOSTGAL2_PHOTOZ']>0:
            # Here we actually need to evaluate the relative probabilities
            self.logger.info('Two potential hosts')
            print(hostinfo)
            sys.exit('does this exist in test data? this is where we fig out what to do')
            return (0.0, 0.0)

        raise DataError('This combination of host data should not exist: {}'.format(hostinfo))



    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    def process(self, datapoint: DataPoint) -> Union[UBson, UnitResult]:
        """

        Parses the provided datapoint for information regarding
        which redshifts should be sampled (and weighted).


        Parameters
        -----------
        datapoint: "DataPoint" instance.

        Returns
        -------
        dict
        """

        dp = datapoint['body']

        # Eventually we assume that we can use the directive ingest setting
        # to make sure that we here get the diaSource object.

        # Guess at relative host galaxy probabilities
        (probH1, probH2) = self.get_hostprob(dp)

        # Depending on the above relative probabilities, choose which redshifts
        z, dz, zsource = 0.0, 0.0, None

        # Final (simulated) data available and used?
        if self.use_final_z:
            if dp.get('REDSHIFT_HELIO',-1)>0:
                z, dz, zsource = dp.get('REDSHIFT_HELIO'), dp.get('REDSHIFT_HELIO_ERR'), 'REDSHIFT_HELIO'
            elif dp.get('REDSHIFT_FINAL',-1)>0:
                z, dz, zsource = dp.get('REDSHIFT_FINAL'), dp.get('REDSHIFT_FINAL_ERR'), 'REDSHIFT_FINAL'

        # Both sources need to be weighted together
        if zsource is None and probH1>0 and probH2>0:
            raise NotImplementedError

        if zsource is None and probH1>0:
            # Use data for first host, spec if available
            if dp.get('HOSTGAL_SPECZ',-1)>0:
                z, dz, zsource = dp.get('HOSTGAL_SPECZ'), dp.get('HOSTGAL_SPECZ_ERR'), 'HOSTGAL_SPECZ'
            else:
                z, dz, zsource = dp.get('HOSTGAL_PHOTOZ'), dp.get('HOSTGAL_PHOTOZ_ERR'), 'HOSTGAL_PHOTOZ'

        # TODO: Retrieve the weighted host galaxy properties and add to t2output

        if zsource is not None:
            # Finish up gaussian case
            t2_output: dict[str, UBson] = {"z_source": zsource,
		                                  "z_samples": [z-1.5*dz, z, z+1.5*dz],
                                          "z_weights": [0.19684199, 0.60631602, 0.19684199],
                                          }
            return t2_output

        # Final cases should be the hostless (default)
        if probH1>0 or probH2>0:
            raise DataError('Unexpected host prob from  {}'.format(dp))


        # Initialize output dict for defult/hostless SNe
        t2_output: dict[str, UBson] = {"z_source": "default",
		                              "z_samples": self.default_zs,
                                      "z_weights": self.default_weights,
                                      }

        return t2_output
