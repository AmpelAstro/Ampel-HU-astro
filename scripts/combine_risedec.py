#!/usr/bin/env python
# coding: utf-8
# Combine different run_risedecline output files

import pandas as pd
import glob
import lcdata

# Base parsnip file to star from
#fname_parsnip = "/home/jnordin/data/noiztf/v240809/ztf_train_bts_train.h5"
fname_parsnip = "/home/jnordin/data/noiztf/v240809/ztf_train_bts_test.h5"
bname = fname_parsnip.split('/')[-1].split('.')[0]
print(bname)


# Create map with classes
print('creating class dict')
typedict = {}
namedict = {}
for fname in [
                  fname_parsnip              
             ]:
                  bts_lc = lcdata.read_hdf5(fname)
                  typedict.update( { row[0]:row[-1] for row in bts_lc.meta.iterrows() } )
                  namedict.update( { row[0]:row[5] for row in bts_lc.meta.iterrows() } )
                  print(len(typedict))

print(set(typedict.values()))

# Create a map for taxonomy (inheriting elasticc)
taxdict = {
    'slsn': 2242,
    'tde': 2243, 
    'sn_ii': 2224, 
    'sn_iin': 2227, 
    'sn_ia': 2222, 
    'sn_ibc': 2223,
    'sn_91bg': 2226,
}

# Combine to large batch of lightcurves
print('Reading and joining')
filestub = 'risedecline_'+bname+'*batch*csv'
print(filestub)
files = glob.glob(filestub)

print(files)
dft = pd.concat([pd.read_csv(fname) for fname in files])

# Change format for bool entries
for boolkind in ['rise','fall','peaked','fastrise','fastfall','hasgaps']:
    dft['bool_'+boolkind] = dft['bool_'+boolkind].astype( bool )


print('adding class column')
dft['class'] = dft['object_id'].apply(lambda x: typedict[x])
print(set(dft['class']))
dft['ztfid'] = dft['object_id'].apply(lambda x: namedict[x])
dft['taxid'] = dft['class'].apply(lambda x: taxdict[x])

# Saving output
dft.to_parquet('risedecline_'+bname+'_combined.parquet')
