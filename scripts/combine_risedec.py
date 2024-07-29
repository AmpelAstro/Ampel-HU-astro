#!/usr/bin/env python
# coding: utf-8
# Combine different run_risedecline output files

import pandas as pd
import glob
import lcdata

# Create map with classes
print('creating class dict')
typedict = {}
namedict = {}
for fname in [
                  "/home/jnordin/data/noiztf/ztf_train_bts_noisified.h5",
                  "/home/jnordin/data/noiztf/ztf_train_bts_test.h5"
              ]:
                  bts_lc = lcdata.read_hdf5(fname)
                  typedict.update( { row[0]:row[-1] for row in bts_lc.meta.iterrows() } )
                  namedict.update( { row[0]:row[6] for row in bts_lc.meta.iterrows() } )
                  print(len(typedict))


# Combine to large batch of lightcurves
print('Reading and joining')
filestub = 'risedec*csv'
files = glob.glob(filestub)

print(files)
dft = pd.concat([pd.read_csv(fname) for fname in files])

print('adding class column')
dft['class'] = dft['object_id'].apply(lambda x: typedict[x])
dft['ztfid'] = dft['object_id'].apply(lambda x: namedict[x])

# Saving output
dft.to_parquet('risedecline_jul28.parquet')
