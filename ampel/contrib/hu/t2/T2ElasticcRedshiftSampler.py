#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2ElasticcRedshiftSampler.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 12.04.2022
# Last Modified Date: 12.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

from typing import Sequence, Union, TypedDict, cast
import numpy as np

from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.content.DataPoint import DataPoint
from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit

class QuantileEntry(TypedDict):
    w: float
    q: list[str]

class RedshiftSamples(TypedDict):
    z_source: str
    host_sep: float
    galaxy_color: dict[str,float] | None
    z_samples: list[float]
    z_weights: list[float]

# Quantile subselections and weights
# Dictionaries of redhift weights together with the accompaning quantile entries.
QUANT_WEIGHTS: dict[int, list[QuantileEntry]] = {
    3: [
        {'w':0.2, 'q': ['000', '010', '020', '030', '040']},
        {'w':0.6, 'q': ['030', '040', '050', '060', '070']},  # Will always be 050?
        {'w':0.2, 'q': ['060', '070', '080', '090', '100']},
        ],
    5: [
        {'w':0.2, 'q': ['000', '010', '020']},
        {'w':0.2, 'q': ['020', '030', '040']},
        {'w':0.2, 'q': ['060', '070', '080']},
        {'w':0.2, 'q': ['040', '050', '060']},
        {'w':0.2, 'q': ['080', '090', '100']},
        ]
    }

# For the gaussian case, which sigma deviations should be used for
# each number of samples?
NORM_SIGMAS: dict[int,list[float]] = {
    1: [0],
    3: [-1.5,0,1.5],
    5: [-2.,-1.,0,1.,2.],
}


def get_elasticc_hostprob(hostinfo: dict) -> tuple[float, float, float]:
    """
    Use the information from the alert info to guess relative prob
    for the two different host galaxies.

    To calculate delta/DLR one would need, except SN and galaxy position,
    some estimate on size, orientation and ellipticity (e.g. moments).
    As these are not contained in the latest test alerts we will for know
    choose the host with the smallest separation / radius.
    """

    # Does sqradius mean what I believe?
    if hostinfo.get('hostgal_sqradius',-99)>=0:
        hostgal_sigsep = hostinfo['hostgal_snsep'] / np.sqrt(hostinfo['hostgal_sqradius'])
    else:
        hostgal_sigsep = -99.
    if hostinfo.get('hostgal2_sqradius',-99)>=0:
        hostgal2_sigsep = hostinfo['hostgal2_snsep'] / np.sqrt(hostinfo['hostgal2_sqradius'])
    else:
        hostgal2_sigsep = -99.

    # SNsep can be exactly zero.
    if hostinfo.get('hostgal_snsep',-99)>=0 and hostinfo.get('hostgal2_snsep',-99)<0:
        # The easy case, only have first galaxy
        return (1.0, 0.0, hostinfo['hostgal_snsep'])
    elif hostinfo.get('hostgal_snsep',-99)<0 and hostinfo.get('hostgal2_snsep',-99)<0:
        # "Hostless" - if these exist
        return (0.0, 0.0, -99.)


    # In principle one could calculate relative weights based on this sigma
    # distances, but skipping for now.
    if hostgal_sigsep<hostgal2_sigsep:
        return (1.0, 0., hostinfo['hostgal_snsep'])
    else:
        return (0.0, 1.0, hostinfo['hostgal2_snsep'])


def get_elasticc_redshift_samples(dp: dict, nbr_samples: int, use_galcol: list[str],
                                    use_final_z: bool=False,
                                    default_zs: list[float] = [0.01, 0.2, 0.5, 0.8],
                                    default_weights: list[float] = [0.4, 0.2, 0.2, 0.2])-> RedshiftSamples:
    """
    Extract a set of samples from the photometric redshift distribution
    of the most likely host galaxy.
    """

    # Guess at relative host galaxy probabilities
    (probH1, probH2, host_sep) = get_elasticc_hostprob(dp)

    # Depending on the above relative probabilities, choose which redshifts
    z, dz, zsource, galcol = 0.0, 0.0, None, None

    # Final (simulated) data available and used?
    if use_final_z:
        # TODO: implement galaxy color determination - but not using final z anyway?
        if dp.get('z_final',-1)>0:
            z, dz, zsource = dp.get('z_final',-1), dp.get('z_final_err', 0), 'Z_FINAL'
        elif dp.get('redshift_helio',-1)>0:
            z, dz, zsource = dp.get('redshift_helio',-1), dp.get('redshift_helio_err', 0), 'REDSHIFT_HELIO'

    # Both sources need to be weighted together
    if zsource is None and probH1>0 and probH2>0:
        raise NotImplementedError

    # Get spectroscopiy host redshift if available
    if zsource is None and probH1>0 and dp.get('hostgal_zspec',-1)>0:
        z, dz, zsource = dp.get('hostgal_zspec',-1), dp.get('hostgal_zspec_err',0), 'HOSTGAL_ZSPEC'
    if zsource is None and probH2>0 and dp.get('hostgal2_zspec',-1)>0:
        z, dz, zsource = dp.get('hostgal2_zspec',-1), dp.get('hostgal2_zspec_err',0), 'HOSTGAL2_ZSPEC'

    # Get galaxy color for most likely host
    # Hostless events will get prefix hostgal2, but if statement should fail
    prefix = "hostgal_mag_" if probH1 > 0 else "hostgal2_mag_"
    galcol = {}
    for thiscol in use_galcol:
        if ( (bluemag:= dp.get(prefix+thiscol[0], 99)) <90 and
            (redmag:= dp.get(prefix+thiscol[1], 99)) <90  ):
            galcol["{}-{}".format(thiscol[0],thiscol[1])] = bluemag-redmag

    # Finish up gaussian case
    if zsource is not None:
        # Calculate weights
        pulls = NORM_SIGMAS[nbr_samples]
        weights = np.exp(-0.5 * np.array(pulls)**2) / np.sqrt(2)
        weights /= np.sum(weights)
        t2_output: RedshiftSamples = {
            "z_source": zsource,
            "host_sep": host_sep,
            "galaxy_color": galcol,
            "z_weights": list(weights),
            "z_samples": [z + p*dz for p in pulls]
        }

        return cast(UBson, t2_output)

    # Extract values from photo-z quantiles
    # (assuming these always exist)
    # Final cases should be the hostless (default)
    if probH1>0 or probH2>0:
        t2_output = {
            "z_source": "HOSTGAL_ZQUANT" if probH1 > 0 else "HOSTGAL2_ZQUANT",
            "host_sep": host_sep,
            "galaxy_color": galcol,
            "z_samples": [],
            "z_weights": []
        }
        prefix = "hostgal_zphot_q" if t2_output["z_source"] == "HOSTGAL_ZQUANT" else "hostgal2_zphot_q"
        # Get some redshifts
        if nbr_samples==1:
            t2_output["z_samples"] = [dp.get(f"{prefix}050", -1)]
            t2_output["z_weights"] = [1.0]
        if nbr_samples==3 or nbr_samples==5:
            t2_output["z_samples"], t2_output["z_weights"] = [], []
            for weight_quantiles in QUANT_WEIGHTS[nbr_samples]:
                weight = weight_quantiles['w']
                quantiles = weight_quantiles['q']
                t2_output["z_samples"].append( np.mean( [qv
                            for q in quantiles if (qv := dp.get(prefix+q,-9)) is not None and qv >-9] ) ) # -9 seems to be null
                t2_output["z_weights"].append(weight)


        return cast(UBson, t2_output)

    # Only left with output dict for defult/hostless SNe
    return cast(UBson, {
            "z_source": "default",
            "host_sep": np.nan,
            "galaxy_color": None,
		    "z_samples": default_zs,
            "z_weights": default_weights,
    })


class T2ElasticcRedshiftSampler(AbsPointT2Unit):
    """
    Parse the elasticc diaSource host information and
    returns a list of redshifts and weights.

    Note 1:
    Elasticc alerts contain both quantiles and mean value + error.
    According to information late May (Kessler), uncertainties will likely
    be normal, meaning that there will be no additional information in the
    quantiles. But they _might_ change this.
    Update from Rob June 10th: They will _likely_ change this.

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

    Update June 14:
    From ghost https://iopscience.iop.org/article/10.3847/1538-4357/abd02b/pdf
    we take it that the first order solution is to take the galaxy with the smallest
    DDLR.

    Note 4:
    The default setting for normal error distribution will be to use
    [-1.5 sig, 0, 1.5 sig] with relative weights [0.19684199, 0.60631602, 0.19684199]


    Note 5:
    Rewrite based on the (lack of) fields in the test alert.

    """

    # How many redshifts samples to return ?
    # Implemented are 1,3 or 5. Would be traightforward to extend to 11
    # based on elasticc quantiles, but that seems excessive.
    nbr_samples: int = 3

    # Default redshifts and weights - in particular used for hostless events
    # Values set semi-randomly now.
    default_zs: list[float] = [0.01, 0.2, 0.5, 0.8]
    default_weights: list[float] = [0.4, 0.2, 0.2, 0.2]

    # Check for final/spectroscopic redshifts, and use if present
    # Usually would be the case, but for testing we might not want to
    use_final_z: bool = False

    # Which bands should be used to estimate the galaxy color?
    # Should these be corrected for MW absorption?
    # Quick trial did not find a strong correlation.
    use_galcol: list[list[str]] = [['u','i']]

    # Can potentially have different metrics for how to choose host galaxy.
    # (could even involve spanning the joint redshifts according to location prob)
    # minDDLR: Select the host with smallest DDLR
    #host_selection: Literal['minDDLR'] = 'minDDLR'






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

        return get_elasticc_redshift_samples(dp, self.nbr_samples, self.use_galcol,
                                            self.use_final_z, self.default_zs, self.default_weights)
