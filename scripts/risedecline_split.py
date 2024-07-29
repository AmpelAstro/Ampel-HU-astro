#!/usr/bin/env python
# coding: utf-8
# Parse a full risedecline feature setup (as from combine_risedec) and:
# - Select subets with specific ndet, time or class ranges.
# - Select max samples per original transient.
# [- Add parsnip information.]
# - Consistently divide into training / validation.
# [- Map class (labels)] Not implemented.


import pandas as pd
import numpy as np

ndetrange = [0,200]
timerange = [0,50]
maxalert = 20


valfrac = 0.2
random_state = 42


df = pd.read_parquet( 'risedecline_jul28.par')
print('startsize', df.shape)

im = (df['ndet']>=ndetrange[0]) & (df['ndet']<=ndetrange[1]) & (df['t_lc']>=timerange[0]) & (df['t_lc']<=timerange[1])
df = df[im]

print('after cuts', df.shape)

stocks = df['ztfid'].unique()
print('original transients', len(stocks))

# Randomize transients for validation sample
valstock = np.random.choice(stocks, size=int(valfrac*len(stocks)), replace=False)

# Go through each transient and subselect (hmm...)
train, validate = [], []
dfgroup = df.groupby('ztfid')
k = 0
for groupstock, group in dfgroup:
  if k % 1000==0:
    print(k, groupstock, group.shape)
    
  if groupstock in valstock:
    savelist = validate
  else: 
    savelist = train
    
  if group.shape[0]>maxalert:
    print('... yep, many alerts from event', group.shape[0])
    savelist.append( group.sample(n=maxalert, random_state=random_state ) )
  else:
    savelist.append( group )
    
  k += 1
  
df_train = pd.concat(train)
df_val = pd.concat(validate)
  
df_train.to_parquet('test_train.parquet')
df_val.to_parquet('test_val.parquet')

print(df_train.shape)
print(df_val.shape)


