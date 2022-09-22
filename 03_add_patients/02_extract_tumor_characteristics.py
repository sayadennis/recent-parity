import os
import sys
import re
import numpy as np
import pandas as pd

din='/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/02_interim'
dout='/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/02_interim'

pathol = pd.read_csv(f'{din}/02_interim/isolated_sections_pathology.csv')
for i in range(pathol.shape[0]):
    for j in range(pathol.shape[1]):
        if pathol.iloc[i,j]=='_':
            pathol.iloc[i,j] = np.nan
        if pathol.iloc[i,j]=='Pending':
            pathol.iloc[i,j] = np.nan

dup = pathol.iloc[[sum(pathol['ir_id']==x)>1 for x in pathol['ir_id']],:]
dup_ids = dup['ir_id'].unique()

def consolidate(df_slice): # slice of data frame. there are 2-8 rows (11/18/2021)
    consolidated = {}
    # get column names
    colnames = df_slice.columns
    # iterate through columns
    isnan = pd.isnull(df_slice)
    for colname in colnames:
        nonnan_slice = df_slice.iloc[~isnan[colname].values,:]
        if np.all(isnan.loc[:,colname]): # all NaN
            continue
        elif np.sum(~isnan.loc[:,colname])==1: # only one is not NaN
            consolidated[colname] = list(df_slice.loc[:,colname].values)[np.where(~isnan.loc[:,colname])[0][0]]
        elif np.sum(~isnan.loc[:,colname])==1: # only one is not NaN
            consolidated[colname] = list(df_slice.loc[:,colname].values)[np.where(~isnan.loc[:,colname])[0][0]]
        else: # multiple non-NaN values
            if np.all([nonnan_slice.loc[i,colname]==nonnan_slice.loc[nonnan_slice.index[0],colname] for i in nonnan_slice.index]): # if all non-NaN values match
                consolidated[colname] = df_slice.loc[df_slice.index[0],colname]
            elif colname in ['Grade (cumulative)', 'Nuclear Pleomorphism/Atypia', 'Mitotic Rate']: # take largest value
                consolidated[colname] = np.max(df_slice[colname].values)
            elif colname=='diagnosis_year': # take earliest year
                consolidated[colname] = np.min(df_slice[colname].values)
    return consolidated

singles = pathol.iloc[[sum(pathol['ir_id']==x)==1 for x in pathol['ir_id']],:]

for dup_id in dup_ids:
    singles = singles.append(consolidate(dup.iloc[dup['ir_id'].values==dup_id,:]), ignore_index=True)

singles.to_csv(f'{dout}/formatted_pathology_results.csv', index=False, header=True)

