#!/usr/bin/env python
# coding: utf-8
# Parse a full risedecline feature setup (as from combine_risedec) and:
# - Select subets with specific ndet, time or class ranges.
# - Select max samples per original transient.
# [- Add parsnip information.]
# [- Map class (labels)] Not implemented.
#
# This version will require parsnip information from file, and will use all of those, For training and validation, you will need to run twice.


import pandas as pd
import numpy as np
from astropy.io.misc.hdf5 import read_table_hdf5


# Filenames
fname_risedecline = 'risedecline_aug13.parquet'
#fname_parsnip = '/home/jnordin/data/noiztf/v240809/noiztfpred_train_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde.h5'
fname_parsnip = '/home/jnordin/data/noiztf/v240809/noiztfpred_val_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde.h5'

bname = fname_parsnip.split('/')[-1].split('.')[0]

print(bname)


# Parameters for selecting which RiseDecline runs to use
ndetrange = [5,20]
timerange = [0,50]
# Maximum number of rows / transient
maxalert = 1
# Cut TDEs?
cut_tde = True

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


print('Parsnip input file shape', df_parsnip.shape)

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





# Go through each transient and subselect (hmm...)
selected = []
dfgroup = df.groupby('ztfid')
k = 0
for groupstock, group in dfgroup:
  if k % 1000==0:
    print(k, groupstock, group.shape)
        
  if group.shape[0]>maxalert:
#    print('... yep, many alerts from event', group.shape[0])
    selected.append( group.sample(n=maxalert, random_state=random_state ) )
  else:
    selected.append( group )
    
  k += 1
  
df_out = pd.concat(selected)

# For this setup, we require parsnip information
pmask = (df_out['model_dof']>0)
print('parsnipmax', len(pmask), sum(pmask))

# Final step - remove columns
df_out = df_out.drop(columns=cut_cols)
print(list(df_out.columns))

df_out.to_parquet('features_noiztf_'+bname+'.parquet')





