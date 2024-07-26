#!/usr/bin/env python
# coding: utf-8
# Run a T2 (T2TabulatorRiseDecline) on dataset provided by lcdata
# Split lightcurves into alerts (i.e. multiple per transient)


import sys
import sncosmo
import lcdata
import numpy as np
from hashlib import blake2b
from ampel.view.ReadOnlyDict import ReadOnlyDict
from bson import encode
from ampel.types import Tag
from ampel.alert.AmpelAlert import AmpelAlert

from ampel.content.T1Document import T1Document
from ampel.view.LightCurve import LightCurve
from ampel.contrib.hu.t2.T2TabulatorRiseDecline import T2TabulatorRiseDecline
from ampel.log.AmpelLogger import AmpelLogger
from ampel.ztf.view.ZTFT2Tabulator import ZTFT2Tabulator

import pandas as pd
from ampel.content.DataPoint import DataPoint




# Settings
fname = "/home/jnordin/data/noiztf/ztf_train_bts_noisified.h5"
tags: list[Tag] = ["NOIZTF","ZTF"]
minsig: float = 3.0
max_transients: int = 1000



# Setting up basic units
ZTF_FILTER_MAP = {
                "ZTF_g": 1, "ZTF_r": 2, "ZTF_i": 3, 
                "ztfg": 1, "ztfr": 2, "ztfi": 3,
                 }
NONFLOAT_DTYPES = {
    "fid": int,
    "band": str,
    "rcid": int,
    "zpsys": str,
    "programid": int,
}
logger = AmpelLogger.get_logger()
t2 = T2TabulatorRiseDecline( logger=logger, significant_bands=['ztfg','ztfr','ztfi'], sigma_det=minsig, t_cadence=3.0 )
t2.post_init()
t2._tab_engines.append( ZTFT2Tabulator() )






# Load the lcdata lightcurves
bts_lc = lcdata.read_hdf5(fname)




def sncosmo2ztpdps(tab, metadict, filtercolumn='band', timeformat='jd', 
                  nonfloat_columns={}):
    """
    Convert an sncosmo table to a list of dps
    as expected for ZTF ingestion.

    Metadict assumed to contain object specific information:
    ra, dec

    Required variable fields:
    - magpsf
    - sigmapsf
    - fid
    - jd
    - programid ("1")
    - rcid ("1")
    - ra
    - dec
    - candid
    """

    if timeformat=='jd':
        tab['jd'] = tab['time']
    else:
        sys.exit('time format not found')
        
    tab['ra'] = metadict['ra']
    tab['dec'] = metadict['dec']
    tab['programid'] = 1
    tab['rcid'] = 1
    # Assuming filter id encoded as in ZTF_FILTER_MAP
    tab['fid'] = [ZTF_FILTER_MAP[b] for b in tab['band']]

    # Convert into list of datapoints
    pps = []
    for row in tab.iterrows():
        pp = {
                k: nonfloat_columns[k](v) if (k in nonfloat_columns and v is not None) else float(v)
                for k, v in zip(tab.colnames, row)
            }

        # It seems magpsf is sometimes not there - fix with updated fluxes etc
        if not pp['magpsf'] >0:
            continue

        
        pp_hash = blake2b(encode(pp), digest_size=7).digest()
        pp["candid"] = int.from_bytes(pp_hash, byteorder=sys.byteorder)
        if pp['candid']>0:
            pp['id'] = pp['candid']
        else:
            sys.exit('hash int negative')

        pps.append(
            DataPoint(id=pp["candid"], tag=tags, 
            body=ReadOnlyDict(pp) )
        )

    # Ensure dps are ordered (should be)
    pps = sorted( pps, key=lambda d: d["body"]['jd'] )


    
    return pps



# Start stepping through transients
results = []
for k in range(len(bts_lc.meta)):
    print(k)

    lc = bts_lc.get_sncosmo_lc(k)
    meta = bts_lc.meta[k]

    # Just for compliance
    t1d = T1Document(stock=meta['object_id'], link=0)
    
    # Convert to (sorted) datapoints
    pps = sncosmo2ztpdps(bts_lc.get_lc(k), bts_lc.meta[k], nonfloat_columns=NONFLOAT_DTYPES)

    # Step through datapoints 
    for l in range(len(pps)):
        alert_pps = pps[0:l+1]
        if (alert_pps[-1]['body']['flux'] / alert_pps[-1]['body']['fluxerr'] )<minsig:
            continue

        # Do t2 processing
        t2out = t2.process(t1d, alert_pps)
        results.append( {'object_id':meta['object_id'], **t2.process(t1d, alert_pps) } )

    if k>max_transients:
        break


# Collect and export
df = pd.DataFrame.from_dict(results)
df.to_csv('foo.csv')



