#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2RunPossis.py
# License:             BSD-3-Clause
# Author:              jnordin@physik.hu-berlin.de
# Date:                11.12.2021
# Last Modified Date:  11.12.2021
# Last Modified By:    jnordin@physik.hu-berlin.de


import numpy as np
import sncosmo # type: ignore[import]
from sfdmap import SFDMap  # type: ignore[import]
from typing import List, Dict, Any, Optional, Tuple, Union, Sequence, Literal
import errno, os, backoff, copy

#from ampel.types import UBson
from ampel.struct.UnitResult import UnitResult
from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo


class T2RunPossis(T2RunSncosmo):
    """

    Load one of the POSSIS models and create an sncosmo_model
    model for fit by T2RunSncosmo.

    """

    # Parameters determining which POSSIS model will be read

    # Currently references to sample model included in Ampel-HU-astro
    possis_dir: str = 'data/kilonova_models'
    model_gen: str = 'bns_m3_3comp'
    mej_dyn: float = 0.01
    mej_wind: float = 0.09
    phi: int = 45
    cos_theta: float = 0.3    # Typically 0., 0.1, ...1.0

    sncosmo_model_name: str = '_'.join(map(str, [model_gen, mej_dyn,
        mej_wind, phi, cos_theta]) )

    # Fix time to specific explosion timestamp
    explosion_time_jd: Optional[float]


    def post_init(self)-> None:
        """
        Retrieve POSSIS model.

        Note that this could be done once at instance init.
        """

        # Find file
        fname = os.path.join(self.possis_dir, self.model_gen,
                     "nph1.0e+06_mejdyn{:05.3f}_mejwind{:05.3f}_phi{}.txt".format(self.mej_dyn, self.mej_wind, self.phi))
        if not os.path.exists(fname):
            self.logger.info('Model not found', extra={'fname':fname})
            return UnitResult(code=DocumentCode.T2_MISSING_INFO)


        # Read model Parameters from first three lines
        with open(fname, 'r') as fh:
            lines = fh.readlines()
            nobs = int(lines[0])
            nwave = int(lines[1])
            line3 = lines[2].split(' ')
            ntime = int(line3[0])
            t_i = float(line3[1])
            t_f = float(line3[2])
            model_cos_theta = np.linspace(0, 1, nobs)  # 11 viewing angles
            phase = np.linspace(t_i, t_f, ntime)  # epochs

        # Limit to one angle
        # Note: U. Feindt developed model where angle was fit, left out for know
        theta_mask = np.isclose(self.cos_theta, model_cos_theta)
        if not sum(theta_mask)==1:
            self.logger.info('Angle not identified',
                extra={'model_cos_theta':model_cos_theta})
            return {}
#            return UnitResult(code=DocumentCode.T2_MISSING_INFO)

        # Read model data
        mdata = np.genfromtxt(fname, skip_header=3)
        wave = mdata[0 : int(nwave), 0] # noqa
        flux = []
        for i in range(int(nobs)):
            flux.append(mdata[i * int(nwave) : i * int(nwave) + int(nwave), 1:]) # noqa
        flux = np.array(flux).T

        # Reduce to one angle
        flux_1angle = flux[:,:,theta_mask].squeeze()
        # Create model
        source = sncosmo.TimeSeriesSource(phase, wave, flux_1angle, name=self.sncosmo_model_name)

        # Setup model, with or without MW correction
        if self.apply_mwcorrection:
            dust = sncosmo.models.CCM89Dust()
            self.sncosmo_model = sncosmo.Model(
                source=source, effects=[dust], effect_names=["mw"], effect_frames=["obs"]
            )
            self.dustmap = SFDMap()
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)
            self.fit_params.remove("mwebv")
        else:
            self.sncosmo_model = sncosmo.Model(source=source)
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)


        # If redshift _should_ be provided we remove this from fit parameters
        if self.redshift_kind is not None or self.backup_z is not None:
            self.fit_params.remove("z")

        # If explosion time should be fixed, do so
        if self.explosion_time_jd is not None:
            self.sncosmo_model.set(t0=self.explosion_time_jd)
            self.fit_params.remove("t0")

        self.default_param_vals = self.sncosmo_model.parameters

        # retry on with exponential backoff on "too many open files"
        self.process = backoff.on_exception(
            backoff.expo,
            OSError,
            giveup=lambda exc: exc.errno != errno.EMFILE,
            logger=self.logger,
            max_time=300,
        )(self.process)
