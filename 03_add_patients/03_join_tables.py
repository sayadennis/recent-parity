import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

din = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing'
dout = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/03_processed'

#####################################
#### Read in demographic columns ####
#####################################

with open(f'{din}/02_interim/column_names_demographic.txt', 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())

########################################
#### Read and lightly clean up data ####
########################################

demo = pd.read_csv(f'{din}/01_raw/demographic_data_patients_wo_genetic_testing.csv', header=None, names=colnames)
pathol = pd.read_csv(f'{din}/02_interim/formatted_pathology_results.csv')
parity = pd.read_csv(f'{din}/02_interim/formatted_parity_results.csv')

# gen_region = pd.read_csv(os.path.join(dn, 'formatted_genetic_testing_results_regions.csv'))
# # gen_region.drop(['palb', 'cdh', 'stk'], axis=1, inplace=True)
# coldict = {}
# for colname in gen_region.columns[1:]:
#     coldict[colname] = colname + '_mutregion'

# gen_region.rename(coldict, axis=1, inplace=True)

# Demo: change birthday column so that it's only showing the day, not time -- e.g. "1970-12-29 00:00:00.00" -> "1970-12-29"
for i in demo.index:
    demo.loc[i,'birth_date'] = demo.loc[i,'birth_date'].split()[0]

# Parity History: remove diagnosis year column because the demographic table has this information already
parity.drop('diagnosis_year', axis=1, inplace=True)
pathol = pathol.drop(['first_name', 'last_name', 'birth_date'],axis=1)
for i in range(demo.shape[0]):
    demo.loc[i,'earliest_diagnosis'] = int(demo.loc[i,'earliest_diagnosis'].split()[0].split('-')[0])

##############################
#### Join all four tables ####
##############################

joined = demo.set_index('ir_id').join(parity.set_index('ir_id'),how='left').join(pathol.set_index('ir_id'), how='left').reset_index()

##############################
#### Add age of diagnosis ####
##############################

joined['age_at_diagnosis'] = None
for i in joined.index:
    birth_year = int(joined.loc[i,'birth_date'][:4])
    diag_year = int(joined.loc[i,'earliest_diagnosis'])
    joined.loc[i,'age_at_diagnosis'] = diag_year - birth_year

#### indicate which list which patient belongs to ####
list1 = list(pd.read_csv(f'{din}/list1_bcdiag.csv')['ir_id'].values)
list2 = list(pd.read_csv(f'{din}/list2_bcdiag_provenc.csv')['ir_id'].values)
list3 = list(pd.read_csv(f'{din}/list3_bcdiag_provenc_gennote.csv')['ir_id'].values)
joined['category'] = None
for i in joined.index:
    if joined.loc[i,'ir_id'] in list1:
        joined.loc[i,'category'] = 'bcdiag'
    elif joined.loc[i,'ir_id'] in list2:
        joined.loc[i,'category'] = 'bcdiag_provenc'
    elif joined.loc[i,'ir_id'] in list3:
        joined.loc[i,'category'] = 'bcdiag_provenc_gennote'


joined = joined[[ # just to reorder so that "age of diagnosis" is right next to "diagnosis year"
    'ir_id', 'category', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name',
       'birth_date', 'gravida', 'para', 'menarche', 'earliest_diagnosis', 'age_at_diagnosis',
       'agefirstpreg', 'ageoflastpreg', 'lastnursed', 'duration', 
       'Specimen Submitted', 'Specimen Dimensions', 'Tumor Size',
       'Histologic Type', 'Grade (cumulative)', 'Tubule Formation',
       'Nuclear Pleomorphism/Atypia', 'Mitotic Rate',
       'Lymphovascular Invasion', 'DCIS as Extensive Intraductal Component',
       'DCIS Measurement/Proportion', 'LCIS', 'Calcifications',
       'Locations of Calcifications', 'Margins of Excision', 'Invasive Cancer',
       'Distance to Margin', 'Margin Widely Free (more than 0.5 cm)', 'DCIS',
       '# Axillary Lymph Nodes Examined', 'Number of Positive Versus Total',
       'Size of Largest Metastasis', 'Extranodal Extension',
       'Breast Tumor Markers', 'ER', 'PR', 'HER2', 'Ki-67', 'p53', 'TNM'
    #    'general', 'brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm', 
    #    'general_mutregion', 'brca_mutregion', 'brca1_mutregion', 'brca2_mutregion', 'palb2_mutregion',
    #    'tp53_mutregion', 'pten_mutregion', 'cdh1_mutregion', 'stk11_mutregion', 'chek2_mutregion', 'atm_mutregion'
]]

print('Total number of patients: %s' % joined.shape[0])
print('Number of patients whose age at diagnosis is 51+: %s' % joined.iloc[joined['age_at_diagnosis'].values>50,:].shape[0])

joined = joined.iloc[joined['age_at_diagnosis'].values<=50,:]

joined.to_csv(f'{dout}/recent_parity_combined_data.csv', index=False, header=True)
