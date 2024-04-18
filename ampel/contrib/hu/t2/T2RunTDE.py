#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel/contrib/hu/t2/T2RunTDE.py
# License:             BSD-3-Clause
# Author:              simeon.reusch@desy.de
# Date:                17.05.2022
# Last Modified Date:  19.04.2023
# Last Modified By:    simeon.reusch@desy.de


import copy
import errno

# from astropy.time import Time
import backoff
import numpy as np
import sncosmo  # type: ignore[import]
from astropy import constants as c
from astropy import units as u
from astropy.units import Quantity
from sfdmap2.sfdmap import SFDMap  # type: ignore[import]

from ampel.contrib.hu.t2.T2RunSncosmo import T2RunSncosmo

# from ampel.enum.DocumentCode import DocumentCode

# from ampel.model.StateT2Dependency import StateT2Dependency

# from urllib.request import urlopen
# from urllib.parse import urljoin


class TDESource(sncosmo.Source):
    _param_names: list = ["risetime", "decaytime", "temperature", "amplitude"]
    param_names_latex: list = [
        "Rise Time~[log10 day]",
        "Decay Time~[log10 day]",
        "Temperature~[log10~K]",
        "Amplitude",
    ]

    def __init__(
        self,
        phase: np.ndarray,
        wave: np.ndarray,
        name: str = "TDE",
        version: str = "1.0",
    ) -> None:
        self.name: str = name  # mandatory for SNCosmo
        self.version: str = version  # mandatory for
        self._phase: np.ndarray = phase
        self._wave: np.ndarray = wave

        # Fit parameters
        # defaults: peaktime = 0, rise=1.584 / decay=2.278 / T=4.0 / peakflux = 1e-25
        self._parameters: np.ndarray = np.array([1.584, 2.278, 4.0, 1e-25])

    @staticmethod
    def _planck_lam(self, wave: np.ndarray, T: np.ndarray) -> np.float64 | np.ndarray:
        """
        Calculate the spectral radiance of a blackbody
        :wave: np.ndarray, array containing wavelength in AA
        :T: np.ndarray, array containing temperatures in K
        """
        wave_aa: Quantity = wave * u.AA
        wave_m: Quantity = wave_aa.to(u.m)

        prefactor = 2 * c.h * c.c**2 / wave_m**5

        prefactor = np.tile(prefactor, (len(T), 1)).transpose()

        exponential_term = (
            c.h.value * c.c.value * 1 / np.outer(wave_m.value, T) / c.k_B.value
        )

        # returns spectral radiance: J s-1 sr-1 m-3
        return prefactor * 1 / (np.exp(exponential_term) - 1) / u.sr

    @staticmethod
    def _cc_bol_lam(self, wave: float | np.ndarray, T: np.ndarray):
        return self._planck_lam(self, wave, T) * u.sr

    @staticmethod
    def _gauss(x: float, sigma: float | np.float64) -> float | np.float64:
        """
        Calculate a Gaussian
        """
        return np.exp(-0.5 * x**2 / (sigma**2))

    @staticmethod
    def _gauss_exp(self, phases: np.ndarray) -> np.ndarray:
        risetime = self._parameters[0]
        decaytime = self._parameters[1]
        peakflux = self._parameters[3]

        # Gaussian rise
        a1 = peakflux
        b1 = 10**risetime
        b2 = 10**decaytime
        a2 = a1 * self._gauss(0, b1)

        phases_rise = phases[(phases <= 0)]
        phases_decay = phases[(phases > 0)]

        vals_rise = a1 * self._gauss(phases_rise, b1)

        phases_decay = phases[(phases > 0)]

        # exponential decay
        vals_decay = a2 * np.exp(-(phases_decay) / b2)

        return np.concatenate((vals_rise, vals_decay))

    def _temp_evolution(self, phase: np.ndarray) -> np.ndarray:
        """
        Create an array with a linear temperature evolution
        """
        temp = self._parameters[2]

        return (10**temp) + phase * 0

    def _flux(self, phase: np.ndarray, wave: np.ndarray) -> np.ndarray:
        """
        Calculate the model flux for given wavelengths and phases
        """
        t_evo = self._temp_evolution(phase=phase)

        phase_iter = np.asarray(phase) if np.ndim(phase) == 0 else phase

        bb_lam = self._cc_bol_lam(self, T=t_evo, wave=wave)

        rise_decay = self._gauss_exp(self, phases=phase_iter)

        model_flux = (rise_decay * bb_lam).transpose()

        return model_flux.to(u.erg / u.s / u.cm**2 / u.AA)


class T2RunTDE(T2RunSncosmo):
    """
    Create a TDE model and fit to a LightCurve object as process is called.

    Create a TDE sncosmo_model for fit by T2RunSncosmo.
    """

    sncosmo_model_name: str = "tde"

    def post_init(self) -> None:
        """
        Create TDE model.
        """

        phase = np.linspace(-50, 100, 20)
        wave = np.linspace(1000, 10000, 5)

        # initialize the TDE source
        tde_source = TDESource(phase, wave, name=self.sncosmo_model_name)

        # Setup model, with or without MW correction
        if self.apply_mwcorrection:
            dust = sncosmo.models.CCM89Dust()
            self.sncosmo_model = sncosmo.Model(
                source=tde_source,
                effects=[dust],
                effect_names=["mw"],
                effect_frames=["obs"],
            )
            self.dustmap = SFDMap()
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)
            self.fit_params.remove("mwebv")
        else:
            self.sncosmo_model = sncosmo.Model(source=tde_source)
            self.fit_params = copy.deepcopy(self.sncosmo_model.param_names)

        # If redshift _should_ be provided we remove this from fit parameters
        if self.redshift_kind is not None or self.backup_z is not None:
            self.fit_params.remove("z")

        self.default_param_vals = self.sncosmo_model.parameters

        # retry on with exponential backoff on "too many open files"
        self.process = backoff.on_exception(  # type: ignore[assignment]
            backoff.expo,
            OSError,
            giveup=lambda exc: not isinstance(exc, OSError)
            or exc.errno != errno.EMFILE,
            logger=self.logger,  # type: ignore[arg-type]
            max_time=300,
        )(self.process)
