import os
import sys
import re
import numpy as np
import pandas as pd

sys.path.append('recent_parity/process_notes')

import ParseNotesParity

dn = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/02_interim'

data = pd.read_csv(f'{dn}/isolated_gynecological_history_section_initial_counsel.csv', sep='!')

#######################
#### Clean up text ####
#######################

for i in range(data.shape[0]):
    data.loc[i,'gyn_history'] = ParseNotesParity.clean_string(data.loc[i,'gyn_history'])

data.drop_duplicates(inplace=True, ignore_index=True)

#################################
#### Read and record results ####
#################################

readfunc_dict = {
    'Initial Consultation' : ParseNotesParity.read_initcons,
    'Consultation 1' : ParseNotesParity.read_con_newpat,
    'Consultation 2' : ParseNotesParity.read_con_newpat,
    'New Patient History and Physical' : ParseNotesParity.read_con_newpat
}

def get_info(ir_id, title_type, notestring):
    readfunc = readfunc_dict[title_type]
    return readfunc(ir_id, title_type, notestring)


gyn_history = pd.DataFrame(
    None, 
    columns=[
        'ir_id', 'title_type', 'gravida', 'para', 'menarche', # 'diagnosis_year', 
         'menopause', 'lmp', 'agefirstpreg', 'ageoflastpreg', 'lastnursed',
         'duration'
    ]
)

for i in range(data.shape[0]):
    ir_id = data.iloc[i,0]
    title_type = data.iloc[i,1]
    year = data.iloc[i,2]
    notestring = data.iloc[i,3]
    resultdict = get_info(ir_id, title_type, notestring)
    resultdict['diagnosis_year'] = year
    gyn_history = gyn_history.append(resultdict, ignore_index=True)

# only select columns that are necessary
gyn_history = gyn_history[[
    'ir_id', 'gravida', 'para', 'menarche', 'diagnosis_year', # 'title_type', 
    'agefirstpreg', 'ageoflastpreg', 'lastnursed', 'duration'
]]

# drop rows where diagnosis year is less than 2010 or over 2020
gyn_history = gyn_history.iloc[[(gyn_history.loc[i,'diagnosis_year'] >= 2010) & (gyn_history.loc[i,'diagnosis_year'] <= 2020) for i in gyn_history.index],:]
# gyn_history = gyn_history.groupby(by='ir_id').min().reset_index()

gyn_history.drop_duplicates(ignore_index=True, inplace=True)
gyn_history = gyn_history.iloc[~np.all(pd.isnull(gyn_history.iloc[:,2:]), axis=1).values,:]
# gyn_history.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/formatted_parity_results.csv', index=False)

##############################################################################################
#### Consolidate rows to create one-to-one relationship between test results and patients ####
##############################################################################################

dup = gyn_history.iloc[[sum(gyn_history['ir_id']==x)>1 for x in gyn_history['ir_id']],:]
dup_ids = dup['ir_id'].unique()

def consolidate(df_slice): # slice of data frame. there are either 2, 3, or 4 rows (6/10/2021)
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
        else: # multiple non-NaN values
            if np.all([nonnan_slice.loc[i,colname]==nonnan_slice.loc[nonnan_slice.index[0],colname] for i in nonnan_slice.index]): # if all non-NaN values match
                consolidated[colname] = df_slice.loc[df_slice.index[0],colname]
            elif colname in ['menarche', 'agefirstpreg', 'ageoflastpreg']: # average the ages
                valuelist = df_slice[colname].values
                try:
                    consolidated[colname] = np.mean(df_slice[colname].iloc[df_slice[colname].values!='unknown'])
                except:
                    print('Could not take mean for menarche for patient %s: %s' % (df_slice['ir_id'].values[0], list(df_slice[colname])))
            elif colname in ['gravida', 'para']: # take largest value
                consolidated[colname] = np.max(df_slice[colname].values)
            elif colname=='diagnosis_year': # take earliest year
                consolidated[colname] = np.min(df_slice[colname].values)
            else:
                if colname!='title_type':
                    # print('Values do not match for patient: %s and feature: %s' % (df_slice['ir_id'].values[0], colname))
                    consolidated[colname] = str(list(df_slice.loc[:,colname].unique()))
                else:
                    consolidated[colname] = None
                # take the positive value. i manually checked which patients came up in this case (there were two, 6/10/2021)
                # consolidated[colname] = 'positive'
    return consolidated

singles = gyn_history.iloc[[sum(gyn_history['ir_id']==x)==1 for x in gyn_history['ir_id']],:]

for dup_id in dup_ids:
    singles = singles.append(consolidate(dup.iloc[dup['ir_id'].values==dup_id,:]), ignore_index=True)

singles.to_csv(f'{dn}/formatted_parity_results.csv', index=False, header=True)


## utility function to see a slice of the dataframe based on ir_id 
# def see_id(table, ir_id):
#     return table.iloc[table['ir_id'].values==ir_id,:]
