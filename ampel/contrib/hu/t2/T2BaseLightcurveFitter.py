#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/contrib/hu/t2/T2BaseLightcurveFitter.py
# License           : BSD-3-Clause
# Author            : jnordin@physik.hu-berlin.de
# Date              : 24.09.2021
# Last Modified Date: 22.04.2022
# Last Modified By  : jnordin@physik.hu-berlin.de

from collections.abc import Sequence

# The following three only used if correcting for MW dust
import extinction  # type: ignore[import]
import numpy as np
import sncosmo  # type: ignore[import]
from astropy.table import Table
from sfdmap2.sfdmap import SFDMap  # type: ignore[import]

from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.base.decorator import abstractmethod
from ampel.content.DataPoint import DataPoint
from ampel.content.T1Document import T1Document
from ampel.contrib.hu.t2.T2DigestRedshifts import T2DigestRedshifts
from ampel.struct.UnitResult import UnitResult
from ampel.types import UBson
from ampel.view.T2DocView import T2DocView


class T2BaseLightcurveFitter(T2DigestRedshifts, AbsTabulatedT2Unit, abstract=True):
    """

    Base class for constructing lightcurve fitters.
    Includes step common to most lightcurve fitters:
    - Obtain a table of flux values
    - Get a redshift (fixed or from catalogs through T2DigestRedshifts.
    - Correct flux for MW reddening
    - Restricting table to phase range as determined from other units


    """

    # Adding default redshift selection values, corresponding to usage of "good" redshifts from catalogs
    redshift_kind = "T2DigestRedshifts"
    max_redshift_category: int = 3

    # Remove MW dust absorption.
    # MWEBV is either derived from SFD maps using the position from light_curve
    # (assuming the SFD_DIR env var is set) OR retrieved from stock (ELASTICHOW?)
    # The default value of Rv will be used.
    # Using this requires extinction, sfdmap and SNCOSMO to be installed. The latter is used to determine effective wavelengths
    apply_mwcorrection: bool = False

    # Phase range usage. Current option:
    # T2PhaseLimit : use the jdmin jdmax provided in this unit output
    # None : use full datapoint range
    phaseselect_kind: None | str

    def post_init(self) -> None:
        """
        Retrieve models and potentially dustmaps.
        """

        if self.apply_mwcorrection:
            self.dustmap = SFDMap()
        # Load e.g. model files as needed
        # e.g. self.model = parsnip.load_model(self.parsnip_model, threads=1)

    def _get_phaselimit(self, t2_views) -> tuple[None | float, None | float]:
        """
        Can potentially also be replaced with some sort of tabulator?

        """

        # Examine T2s for eventual information
        jdstart: None | float = None
        jdend: None | float = None

        if self.phaseselect_kind is None:
            jdstart = -np.inf
            jdend = np.inf
        else:
            for t2_view in t2_views:
                # So far only knows how to parse phases from T2PhaseLimit
                if t2_view.unit != "T2PhaseLimit":
                    continue
                self.logger.debug(f"Parsing t2 results from {t2_view.unit}")
                t2_res = (
                    res[-1] if isinstance(res := t2_view.get_payload(), list) else res
                )
                jdstart = t2_res["t_start"]
                jdend = t2_res["t_end"]

        return jdstart, jdend

    def _deredden_mw_extinction(self, ebv, phot_tab, rv=3.1) -> Table:
        """
        For an input photometric table, try to correct for mw extinction.
        Resuires extinction & sncosmo to be loaded, and that sncosmo knows the band wavelength.
        """

        # Find effective wavelength for all filters in phot_tab
        filterlist = set(phot_tab["band"])
        eff_wave = [sncosmo.get_bandpass(f).wave_eff for f in filterlist]

        # Determine flux correction (dereddening) factors
        flux_corr = 10 ** (0.4 * extinction.ccm89(np.array(eff_wave), ebv * rv, rv))

        # Assign this appropritately to Table
        phot_tab["flux_original"] = phot_tab["flux"]
        phot_tab["fluxerr_original"] = phot_tab["fluxerr"]
        for k, band in enumerate(filterlist):
            phot_tab["flux"][(phot_tab["band"] == band)] *= flux_corr[k]
            phot_tab["fluxerr"][(phot_tab["band"] == band)] *= flux_corr[k]

        return phot_tab

    def get_fitdata(
        self,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> tuple[Table | None, dict]:
        """

        Obtain data necessary for fit:
        - Flux table, possiby constrained in time.
          (corrected for MW extinction if set)
        - Redshift. Possibly multiple values, possibe None

        Returns
        -------
        dict
        """

        # Fit data info
        fitdatainfo: dict[str, UBson] = {}

        # Check for phase limits
        (jdstart, jdend) = self._get_phaselimit(t2_views)
        fitdatainfo["jdstart"] = jdstart
        fitdatainfo["jdend"] = jdend
        if fitdatainfo["jdstart"] is None:
            return (None, fitdatainfo)

        # Obtain photometric table
        sncosmo_table = self.get_flux_table(datapoints)
        sncosmo_table = sncosmo_table[
            (sncosmo_table["time"] >= jdstart) & (sncosmo_table["time"] <= jdend)
        ]

        # Potentially correct for dust absorption
        # Requires filters to be known by sncosmo (and the latter installed)
        if self.apply_mwcorrection:
            # Get ebv from coordiantes.
            # Here there should be some option to read it from journal/stock etc
            mwebv = self.dustmap.ebv(*self.get_pos(datapoints, which="mean"))
            fitdatainfo["mwebv"] = mwebv
            sncosmo_table = self._deredden_mw_extinction(mwebv, sncosmo_table)

        ## Obtain redshift(s) from T2DigestRedshifts
        zlist, z_source, z_weights = self.get_redshift(t2_views)
        fitdatainfo["z"] = zlist
        fitdatainfo["z_source"] = z_source
        fitdatainfo["z_weights"] = z_weights
        # A source class of None indicates that a redshift source was required, but not found.
        if not isinstance(zlist, list) or z_source is None:
            return (None, fitdatainfo)

        return (sncosmo_table, fitdatainfo)

    # ==================== #
    # AMPEL T2 MANDATORY   #
    # ==================== #
    @abstractmethod
    def process(
        self,
        compound: T1Document,
        datapoints: Sequence[DataPoint],
        t2_views: Sequence[T2DocView],
    ) -> UBson | UnitResult:
        """

        Fit a model to the lightcurve of this transient.
        See T2DemoLightcurveFitter

        Returns
        -------
        dict
        """

        raise NotImplementedError
        return None
