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

    How should we do with potentially wrong host association?
    However we hope this can be done, we assume the work in selecting
    is done in this unit.

    It actually turns out that this is now stored as a datapoint, so we will
    need to go through these. Since the ordering is not known we will have to
    parse all!

    """

    # How many redshifts to return (up to some maximum)
    # (prop not all values allowed)
    max_zs: int = 3
    # Hostless redshifts: which redshifts and weights to use if no host found.
    # Make it depend on host brightness? Set based on training sample distribution.
    default_zs: Sequence[float] = [0.2, 0.4, 0.6, 1]
    default_weights: Sequence[float] = [0.1, 0.2, 0.5, 0.2]




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

        # Eventually we assume that we can use the directive ingest setting
        # to make sure that we here get the diaSource object.
        print(datapoint)
        print(datapoint['body'])

        # Initialize output dict
        t2_output: dict[str, UBson] = {"z_source": "default",
		                              "z_samples": self.default_zs,
                                      "z_weights": self.default_weights,
                                      }

        return t2_output
