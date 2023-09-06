#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File:                ampel-hu-astro/ampel/contrib/hu/t2/T2HealpixProb.py
# License:             BSD-3-Clause
# Author:              jn <jnordin@physik.hu-berlin.de>
# Date:                28.03.2018
# Last Modified Date:  28.03.2021
# Last Modified By:    jn <jnordin@physik.hu-berlin.de>

from typing import Any, ClassVar, Literal, Iterable, Union
from functools import cached_property
from ampel.types import UBson

from extcats import CatalogQuery
from pymongo import MongoClient
from pymongo.errors import AutoReconnect
import backoff


from astropy.coordinates import SkyCoord
from astropy.table import Table
from extcats.catquery_utils import get_closest, get_distances
from numpy import asarray, degrees, mean
import numpy as np

from ampel.abstract.AbsPointT2Unit import AbsPointT2Unit
from ampel.content.DataPoint import DataPoint
from ampel.enum.DocumentCode import DocumentCode
from ampel.struct.UnitResult import UnitResult
from ampel.struct.Resource import Resource
from ampel.model.DPSelection import DPSelection
from ampel.content.T1Document import T1Document
from ampel.abstract.AbsStateT2Unit import AbsStateT2Unit
from ampel.abstract.AbsTabulatedT2Unit import AbsTabulatedT2Unit
from ampel.contrib.hu.util.AmpelHealpix import AmpelHealpix, deres


class T2HealpixProb(AbsStateT2Unit, AbsTabulatedT2Unit):
    """
    Retrive cumulative probability for object to be associated to Healpix map.
    Will be read using AmpelHealpix. 

    Probability is the _min prob contour needed to include position_. 

    Assumes that (at least) two kinds of datapoints exists:
    - One describing details regarding the Healpix map we wish to use.
    - At least one describing an optical alert, with Ra+DEC
    
    Assumes Healpix map retrievable through AmpelHealpix.
    {map_name, map_dir, map_url, map_hash} given directly.
    
    Problem: as two different datapoints are needed this can no longer be an AbsPointT2Unit. In practice this might not change much
    since we do not expect a lot of datapoints per object, but it could lead to problems. Keep in mind
    
    """

    # Name (signifier)
    map_name: str       # Assumed to agree with map name 
    healpix_map: None | AmpelHealpix = None
    map_hash: None | str = None

    pvalue_limit: float = 0.9

    # URL for healpix retrieval
#    map_url: None | str = None
#    map_dir: None | str = None   # Will first look here

    def load_map(self, map_info):

        # Load and process map
        self.healpix_map = AmpelHealpix(map_name=map_info.get('map_name'), map_url=map_info.get('map_url'), save_dir = map_info.get('map_dir'))
        self.map_hash = self.healpix_map.process_map()
        if not self.map_hash==map_info['hash']:
            raise ValueError("Healpix hash changed - modified map?")
        

    def process(self,
		compound: T1Document,
		datapoints: Iterable[DataPoint],
	) -> Union[UBson, UnitResult]:
        """
        :returns: cumulative probability in Healpix:

        {
            'cumprob': 0.666,
            'map_name': 'HealpixID14',
            'map_hash': 'x',
            "map_dist": 123.3,
            "map_dist_unc": 40.3,
        }
        """
        
        # Loading the map not during init, but as first alert is processed, as the map id is not initially known.
        # Should be in a tabulator?
        # 1. Look for the healpix dp
        # Alternatively, one could look for this in tags
        healdps = [dp for dp in datapoints if self.map_name == dp['body'].get('map_name') ]
        if not len(healdps)==1:
            raise ValueError("Zero or multiple Healpix maps connected to state.")         
        # 2. Load map if not already loaded, if so check that hash is still the same.
        if self.map_hash:
            if not self.map_hash == healdps[0]['body'].get('hash'):
                raise ValueError("Hash changed w.r.t. loaded map")         
        else:
            self.load_map(healdps[0]['body'])


        #print("HEALPIXPROB DDDDDDDDDDDDDDDDDDDDDDDDDDDD::", healdps[0].items())
	    
	# 3. Otherwise, find the position of the max lum dp.
        pos = self.get_positions(datapoints)   # (jd, ra, dec)
	
	# 4. Use this to return the prob. 
        out_dict: dict[str, Any] = {'map_name': self.map_name, 'map_hash': self.map_hash, 'trigger_time': self.healpix_map.trigger_time}
        out_dict['cumprob'] = self.healpix_map.get_cumprob(mean([dp[1] for dp in pos]), mean([dp[2] for dp in pos]))

    # 5. Propagate used area, total alerts (unfiltered)
        # Combine pixels when possible
        # pixels = self.healpix_map.get_pixelmask(self.pvalue_limit)
        # deresdict = deres(self.healpix_map.nside, pixels)
        # healpix_regions = [
        #     {"nside": nside, "pixels": members} for nside, members in deresdict.items()
        # ]

        # hp_area = 0
        # for region in healpix_regions:
        #     npix_from_nside = 12 * region["nside"]**2
        #     hp_area += len(region["pixels"]) / npix_from_nside
        # hp_area *= 360**2 / np.pi
        # print("HEALPIX AREA:", hp_area)
        out_dict["map_area"] = healdps[0]['body']["map_area"]
        out_dict["unfiltered_alerts"] = healdps[0]['body']["alert_count_nofilter"]
        out_dict["queried_alerts"] = healdps[0]['body']["alert_count_query"]
    
    # 6. return map distance & uncertainty
        map_dist, map_dist_unc = self.healpix_map.get_mapdist()
        out_dict["map_dist"] = map_dist
        out_dict["map_dist_unc"] = map_dist_unc


        return out_dict
