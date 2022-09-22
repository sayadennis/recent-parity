import os
import numpy as np
import pandas as pd
import re
from datetime import datetime

din='/share/fsmresfiles/breast_cancer_pregnancy/explore_cohort'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

colnames=['ir_id', 'Epic_MRN', 'Powerchart_MRN', 'first_name', 'last_name', 'dob', 'earliest_diagnosis']

bc = pd.read_csv(f'{din}/cohort_by_breastcanceronly.csv', header=None, names=colnames)
prov = pd.read_csv(f'{din}/cohort_by_prov_name.csv', header=None, names=colnames)
genn = pd.read_csv(f'{din}/cohort_by_genetic_notes.csv', header=None, names=colnames)

## Clean data
for df in [bc, prov, genn]:
    # add age of diagnosis
    df['age_of_dx'] = [(datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f') - datetime.strptime(y, '%Y-%m-%d %H:%M:%S.%f')).days // 365 for (x,y) in zip(df['earliest_diagnosis'].values, df['dob'].values)]
    # add year of diagnosis
    df['earliest_diagnosis'] = [datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f').date() for x in df['earliest_diagnosis'].values]
    # reformat dob to exclude time (only year, month, and day)
    df['dob'] = [datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f').date() for x in df['dob'].values]
    # eliminate unnecessary columns and remove duplicates
    # df.drop(['Epic_MRN', 'Powerchart_MRN', 'earliest_diagnosis'], axis=1, inplace=True)
    df.drop_duplicates(ignore_index=True, inplace=True)


#########################################################################
#### Get lists of patients who should have genetic testing but don't ####
#########################################################################

## List 1: breast cancer diagnosis at age <= 35 but neither provider encounter nor genetic note 
# get patients who are in BC but not in PROV or GENN
bconly = []
for ir_id in bc['ir_id'].unique():
    if np.all((ir_id not in prov['ir_id'].values) & (ir_id not in genn['ir_id'].values)):
        bconly.append(ir_id)

bconly = bc.iloc[[ir_id in bconly for ir_id in bc['ir_id'].values],:]
bconly = bconly.iloc[bconly['age_of_dx'].values<=35,:] # only select ones whose earliest diagnosis is at age 35 or younger 
bconly = bconly.iloc[bconly['year_of_dx'].values<=2019,:] # only select ones who were diagnosed before 2019 (to ensure enough time for genetic testing) 

bconly.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/list1_bcdiag.csv')

## List 2: breast cancer diagnosis at age <= 35 and has provider encounter but no genetic note

bcprov = []
for ir_id in prov['ir_id'].unique():
    if ir_id not in genn['ir_id'].values:
        bcprov.append(ir_id)

bcprov = prov.iloc[[ir_id in bcprov for ir_id in prov['ir_id'].values],:]
bcprov = bcprov.iloc[bcprov['age_of_dx'].values<=35,:] # only select ones whose earliest diagnosis is at age 35 or younger 

bcprov.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/list2_bcdiag_provenc.csv')

## List 3: breast cancer diagnosis at age <= 35 and has provider encounter and genetic note but no test results

struc = pd.read_csv(
    f'{din}/exploratory_notes_genetic.csv', header=None,
    names=['ir_id', 'Epic_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'dob', 'sid1', 'sid2', 'sid3', 'source_table', 'created', 'updated', 'note_text']
)

separated_sections = pd.DataFrame(None, columns=['ir_id', 'test_results']) # , 'interpretation' 'created_year'

## Loop through patient ir_id's and isolate sections 
for i in struc['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(struc['note_text'].iloc[[x==i for x in struc['ir_id']]].values)
    # subsetyears = list(struc['created_year'].iloc[[x==i for x in struc['ir_id']]].values)
    for j in range(len(subsetlines)):
        # find first occurrence of 'Test results:'
        testsearch = re.split('Test results:', subsetlines[j], 1)
        if len(testsearch) == 1: # if there is no 'Test results:'
            continue
        else:
            testresults_text = testsearch[1]
            intesearch = re.split('Interpretation:', testresults_text, 1)
            if len(intesearch) == 1: # if there is no 'Interpretation:', ignore this instance of 'Test results:'
                continue
            else:
                # save the sliced out 'Test results:' section 
                testresults = intesearch[0].strip()
                separated_sections = separated_sections.append({
                    'ir_id' : i, 
                    # 'created_year' : int(subsetyears[j]),
                    'test_results' : testresults
                }, ignore_index=True)

separated_sections.drop_duplicates(inplace=True, ignore_index=True)

genfound = list(separated_sections['ir_id'].unique()) # patient IDs for those whose genetic testing results were identified 

bcprovgen = []
for ir_id in genn['ir_id'].unique():
    if ir_id in prov['ir_id'].values:
        if ir_id not in genfound:
            bcprovgen.append(ir_id)

bcprovgen = genn.iloc[[ir_id in bcprovgen for ir_id in genn['ir_id'].values],:]
bcprovgen = bcprovgen.iloc[bcprovgen['age_of_dx'].values<=35,:] # only select ones whose earliest diagnosis is at age 35 or younger 

bcprovgen.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/list3_bcdiag_provenc_gennote.csv')
