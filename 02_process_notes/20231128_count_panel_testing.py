import os
import sys
import re
import pickle
from datetime import datetime

import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/01_ssms_raw'

with open(os.path.join(dn, 'column_names_note_type_detail.txt'), 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())

##################################
#### Genetic counseling notes #### 
##################################

## Read in SQL query results
struc = pd.read_csv(
    os.path.join(dn, 'genetic/sql_query_results_notes_genetic_counsel.csv'), 
    header=None, names=colnames
)

## Store isolated sections in this dataframe 
separated_sections = pd.DataFrame(None, columns=['ir_id', 'test_results']) # , 'interpretation' 'created_year'

patterns = [
    "Test results:",
    "Interpretation:",
    "Risk Assessment"
]

texts = pd.DataFrame('', index=struc['ir_id'].unique(), columns=patterns[:-1])

## Loop through patient ir_id's and isolate sections 
for i in struc['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(struc['note_text'].iloc[[x==i for x in struc['ir_id']]].values)
    patient_notelist = []
    # subsetyears = list(struc['created_year'].iloc[[x==i for x in struc['ir_id']]].values)
    for j in range(len(subsetlines)):
        # find first occurrence of 'Test results:'
        testsearch = re.split('Test results:', subsetlines[j], 1)
        if len(testsearch) == 1: # if there is no 'Test results:'
            continue
        else:
            for ix, pattern in enumerate(patterns[:-1]):
                sliced_string = subsetlines[j].split(patterns[ix])[-1]
                sliced_string = sliced_string.split(patterns[ix+1])[0]
                if len(sliced_string) > len(texts.loc[i,pattern]):
                    texts.loc[i,pattern] = sliced_string


genelist = ['BRCA', 'BRCA1', 'BRCA2', 'PALB2', 'TP53', 'PTEN', 'CDH1', 'STK11', 'CHEK2', 'ATM']

mention = pd.DataFrame(0, index=struc['ir_id'].unique(), columns=genelist)

for i in texts.index:
    for genename in genelist:
        if genename in texts.loc[i,'Test results:']:
            mention.loc[i,genename] = 1

##################################
#### Get REDCap data to aling ####
##################################

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'

redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

# remove exluded patients
redcap = redcap.iloc[(redcap.exclude_demo.values!=1) & (redcap.exclude_tum.values!=1),:]

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# mark family history as none for patients who are missing it 
redcap.fam_hx = redcap.fam_hx.fillna(0.)

redcap = redcap.iloc[(redcap.exclude_demo.values!=1) & (redcap.exclude_tum.values!=1),:]

include_ir_ids = []

struc.birth_date = [x.split()[0] for x in struc.birth_date]

for i in redcap.index:
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    subsetlines = struc.iloc[(struc.EPIC_mrn.values==epic_mrn) & (struc.birth_date.values==dob)]
    if subsetlines.shape[0]>0:
        include_ir_ids.append(subsetlines['ir_id'].iloc[0])

mention = mention.loc[include_ir_ids,:]

mention.iloc[mention.BRCA.values==1,:].sum(axis=0)/606

# Estimating what proportion of patients had panel testing: of the patients that have genes mentioned in their genetic testing results section, what proportion has mention of any other gene?  
np.any(mention.iloc[:,3:], axis=1).sum()/606


