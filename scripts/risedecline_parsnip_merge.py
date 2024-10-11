#!/usr/bin/env python
# coding: utf-8
# Parse a full risedecline feature setup (as from combine_risedec) and:
# - Select subets with specific ndet, time or class ranges.
# - Select max samples per original transient.
# - Add parsnip information.
# - Consistently divide into training / validation if so wished.
# [- Map class (labels)] Not implemented.


import pandas as pd
import numpy as np
from astropy.io.misc.hdf5 import read_table_hdf5
import matplotlib.pyplot as plt


training = True	
mode = 5 	

# Filenames
if training:
  # Fraction to set aside for validation
  # train
  fname_risedecline = 'risedecline_aug13.parquet'
  fname_parsnip = '/home/jnordin/data/noiztf/v240809/noiztfpred_train_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde.h5'
  fname_out = 'noiztfpred_train_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde'
  valfrac = 0.0
else:
  # test
  # Fraction to set aside for validation
  fname_risedecline = 'risedecline_ztf_train_bts_test_combined.parquet'
  fname_parsnip = '/home/jnordin/data/noiztf/v240809/noiztfpred_val_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde.h5'
  fname_out = 'noiztfpred_val_sncut_scale5_z3_cadscale0.5_subamp0.9_seed0_notde'
  valfrac = 1.0


# Parameters, usually overwritten below
ndetrange = [0,50]
timerange = [0,999]
maxalert = 10
# Only use final alert for each transient (to match that parsnip only lives here)
only_last = False
# Cut down to one classifier (None, parsnip or risedecline)
classifier = None


# Kinds of test samples created
# Varying [alertsize]x[classifier]x[archive]

if mode==1:
  ## Only early data, risedecline
  ndetrange = [0,4]
  timerange = [0,10]
  maxalert = 10
  only_last = False
  classifier='risedecline'
  fname_out = fname_out+'_early'
elif mode==2:
  # archive dataset both classifiers, 
  ndetrange = [5,999]
  timerange = [0,999]
  only_last = True
  maxalert = 10 # should not matter
  classifier = None # both of them
elif mode==3:
  # archive dataset, parsnip
  ndetrange = [5,999]
  timerange = [0,999]
  only_last = True
  maxalert = 10 # should not matter
  classifier = 'parsnip' # both of them
elif mode==4:
  # archive dataset, risedecline
  ndetrange = [5,999]
  timerange = [0,999]
  only_last = True
  maxalert = 10 # should not matter
  classifier = 'risedecline' # both of them
elif mode==5:
  # mixed dataset, risedecline
  ndetrange = [5,50]
  timerange = [0,999]
  only_last = False
  maxalert = 10 
  classifier = 'risedecline' # both of them
  fname_out = fname_out+'_late'


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



df = pd.read_parquet( fname_risedecline )
print('startsize', df.shape)

# TDEs might or might not be included in the feature sample
if cut_tde:
  df = df[ ~(df['taxid']==2243) ]
  print('after cutting tdes', df.shape)

im = (df['ndet']>=ndetrange[0]) & (df['ndet']<=ndetrange[1]) & (df['t_lc']>=timerange[0]) & (df['t_lc']<=timerange[1])
df = df[im]

if only_last:
  #dfg = df.groupby('object_id')
  #idx = dfg['ndet'].idxmax()	
  df = df.sort_values(['object_id', 'ndet']).groupby('object_id').tail(1)
  fname_out = fname_out + '_lastalert'


print('after cuts', df.shape)

# Join with parsnip information
df = pd.merge(df, df_parsnip, on="object_id",how="left")




# Cut sample if only including one classifier
if classifier=='risedecline':
  print('cutting parsnip')
  df = df.drop(columns=[p for p in parsnip_cols if not p=='object_id'])
  fname_out = fname_out + '_risedecline'
elif classifier=='parsnip':
  print('cutting risedecline')
  print(df.shape)
  df = df[ [p for p in parsnip_cols if not p=='object_id'] + ['taxid'] + cut_cols ]
  print(df.shape)
  fname_out = fname_out + '_parsnip'
else:
  fname_out = fname_out + '_risdecpar'

stocks = df['ztfid'].unique()
print('original transients', len(stocks))

simid = df['object_id'].unique()
print('simulated transients', len(simid))




# Go through each transient and subselect (hmm...)
train, validate = [], []


# Randomize transients for validation sample	
#valstock = np.random.choice(stocks, size=int(valfrac*len(stocks)), replace=False)
# This groups according to original transient
#dfgroup = df.groupby('ztfid')

# But we rather wish to sample from the ztf transients?
valstock = np.random.choice(simid, size=int(valfrac*len(simid)), replace=False)
dfgroup = df.groupby('object_id')

print('chosen for validation', len(valstock))



k = 0
for groupstock, group in dfgroup:
  if k % 1000==0:
    print(k, groupstock, group.shape)
    
  if groupstock in valstock:
    savelist = validate
  else: 
    savelist = train
    
  if group.shape[0]>maxalert:
    #print('... yep, many alerts from event', group.shape[0])
    savelist.append( group.sample(n=maxalert, random_state=random_state ) )
  else:
    savelist.append( group )
    
  k += 1


# We now have a mixed setup (noo...) where we also have a flag for complete
# Standard procedure is to only have training or validation samples...
if training:
  train = train+validate
  validate = []
else:
  validate = train+validate
  train = []


if len(train)>0:
  print('storing training data for events:', len(train))  
  df_train = pd.concat(train)
  df_train = df_train.drop(columns=cut_cols)
  df_train.to_parquet('features_'+fname_out+'_train.parquet')
else:
  print('No training set created')

if len(validate)>0:
  print('storing training data for events:', len(validate))  
  df_val = pd.concat(validate)
  df_val = df_val.drop(columns=cut_cols)
  df_val.to_parquet('features_'+fname_out+'_validate.parquet')
else:
  print('No val set created')







