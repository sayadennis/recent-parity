import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

#####################################
#### Read in demographic columns ####
#####################################

with open(os.path.join(dn, 'column_names_demographic.txt'), 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())

########################################
#### Read and lightly clean up data ####
########################################

demo = pd.read_csv(os.path.join(dn, 'sql_query_results_patient_demographic.csv'), header=None, names=colnames)
gen = pd.read_csv(os.path.join(dn, 'formatted_genetic_testing_results.csv'))
par = pd.read_csv(os.path.join(dn, 'formatted_parity_results.csv'))

gen_region = pd.read_csv(os.path.join(dn, 'formatted_genetic_testing_results_regions.csv'))
# gen_region.drop(['palb', 'cdh', 'stk'], axis=1, inplace=True)
coldict = {}
for colname in gen_region.columns[1:]:
    coldict[colname] = colname + '_mutregion'

gen_region.rename(coldict, axis=1, inplace=True)

# Demo: change birthday column so that it's only showing the day, not time -- e.g. "1970-12-29 00:00:00.00" -> "1970-12-29"
for i in demo.index:
    demo.loc[i,'birth_date'] = demo.loc[i,'birth_date'].split()[0]

# Demo: take diagnosis year to be the smallest value (earliest year) for duplicate rows 
demo.drop('diagnosis_year', axis=1, inplace=True)
# for ir_id in demo['ir_id'].unique():
#     indices = demo.iloc[demo['ir_id'].values==ir_id,:].index
#     demo.loc[indices,'diagnosis_year'] = np.min(demo.loc[indices,'diagnosis_year'])

demo.drop_duplicates(ignore_index=True, inplace=True)
# demo = demo.iloc[:,[x!='diagnosis_year' for x in demo.columns]].drop_duplicates(ignore_index=True)

# Genetic: remove duplicate gene names e.g. 'palb' and 'palb2'
# gen.drop(['palb', 'cdh', 'stk'], axis=1, inplace=True)

# Parity History: remove diagnosis year column because the demographic table has this information already
# par.drop('diagnosis_year', axis=1, inplace=True)

##############################
#### Join all four tables ####
##############################

joined = demo.set_index('ir_id').join(par.set_index('ir_id'),how='inner').join(gen.set_index('ir_id'), how='left').join(gen_region.set_index('ir_id'), how='left').reset_index()

##############################
#### Add age of diagnosis ####
##############################

joined['age_at_diagnosis'] = None
for i in joined.index:
    birth_year = int(joined.loc[i,'birth_date'][:4])
    diag_year = int(joined.loc[i,'diagnosis_year'])
    joined.loc[i,'age_at_diagnosis'] = diag_year - birth_year

joined = joined[[ # just to reorder so that "age of diagnosis" is right next to "diagnosis year"
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name',
       'birth_date', 'gravida', 'para', 'menarche', 'diagnosis_year', 'age_at_diagnosis',
       'agefirstpreg', 'ageoflastpreg', 'lastnursed', 'duration', 'general',
       'brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11',
       'chek2', 'atm', 'general_mutregion', 'brca_mutregion',
       'brca1_mutregion', 'brca2_mutregion', 'palb2_mutregion',
       'tp53_mutregion', 'pten_mutregion', 'cdh1_mutregion', 'stk11_mutregion',
       'chek2_mutregion', 'atm_mutregion'
]]

print('Total number of patients: %s' % joined.shape[0])
print('Number of patients whose age at diagnosis is 51+: %s' % joined.iloc[joined['age_at_diagnosis'].values>50,:].shape[0])

joined = joined.iloc[joined['age_at_diagnosis'].values<=50,:]

joined.to_csv(os.path.join(dn, 'recent_parity_combined_data.csv'), index=False, header=True)

