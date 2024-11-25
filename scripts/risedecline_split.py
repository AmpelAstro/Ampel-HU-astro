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
from astropy.io.misc.hdf5 import read_table_hdf5


# Filenames
fname_risedecline = 'risedecline_aug13.parquet'
fname_parsnip = '/home/jnordin/data/noiztf/v240724/btsnoisyztf_predictions_scale3_subsample90_iin_slsn_val_v5_split.h5'

# Parameters for selecting which RiseDecline runs to use
ndetrange = [5,20]
timerange = [0,50]
# Maximum number of rows / transient
maxalert = 20
# Cut TDEs?
cut_tde = True
# Fraction to set aside for validation
valfrac = 0.2

# Columns to transfer from parsnip results
parsnip_cols = ['object_id','color', 'color_error', 'amplitude', 'amplitude_error', 's1', 's1_error', 's2', 's2_error', 's3', 's3_error', 'total_s2n', 'count', 'count_s2n_3', 'count_s2n_5', 'count_s2n_3_pre', 'count_s2n_3_rise', 'count_s2n_3_fall', 'count_s2n_3_post', 'model_chisq', 'model_dof', 'luminosity', 'luminosity_error', 'chisqdof']
# Column to cut from final feature sets
cut_cols = ['Unnamed: 0', 'object_id', 'jd_det', 'jd_last', 'jd_min', 'jd_peak_ztfr', 'jd_peak_ztfg', 'jd_peak_ztfi', 'success', 'jd_peak', 'cause', 'class', 'ztfid','alldet','nnegdet','frac_pos']
# Random seed when dividing into training/validation
random_state = 42



# Collect parsnip information
parsnip_tab = read_table_hdf5( fname_parsnip )

# Add specific chi/dof column
# Only use results with at least one degree of freedom...
parsnip_tab = parsnip_tab[ (parsnip_tab['model_dof']>0) ]
parsnip_tab['chisqdof'] = parsnip_tab['model_chisq'] / parsnip_tab['model_dof']
parsnip_tab = parsnip_tab[parsnip_cols]

df_parsnip = parsnip_tab.to_pandas()
df_parsnip['object_id'] = df_parsnip['object_id'].str.decode('utf-8')



df = pd.read_parquet( fname_risedecline )
print('startsize', df.shape)

# TDEs might or might not be included in the feature sample
if cut_tde:
  df = df[ ~(df['taxid']==2243) ]
  print('after cutting tdes', df.shape)

im = (df['ndet']>=ndetrange[0]) & (df['ndet']<=ndetrange[1]) & (df['t_lc']>=timerange[0]) & (df['t_lc']<=timerange[1])
df = df[im]

print('after cuts', df.shape)

# Join with parsnip information
df = pd.merge(df, df_parsnip, on="object_id",how="left")


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

df_train = df_train.drop(columns=cut_cols)
df_val = df_val.drop(columns=cut_cols)

# Final step - remove columns
print(list(df_val.columns))

df_train.to_parquet('features_noiztf_v0_train.parquet')
df_val.to_parquet('features_noiztf_v0_validate.parquet')





