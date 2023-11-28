import os
import sys
import re
import numpy as np
import pandas as pd
from datetime import datetime

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
separated_sections.to_csv(os.path.join(dn, 'genetic_counsel/isolated_sections_genetic_testing.csv'), sep='!', header=True, index=False)


###############################
#### Initial counsel notes #### 
###############################

## Read in SQL query results
struc = pd.read_csv(
    os.path.join(dn, 'initialcounsel/sql_query_results_notes_initial_counsel.csv'),
    header=None, names=colnames
)

## Store isolated sections in this dataframe 
separated_sections = pd.DataFrame(None, columns=['ir_id', 'title_type', 'created_year', 'gyn_history'])

## Loop through patient ir_id's and isolate sections 
for i in struc['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(struc['note_text'].iloc[[x==i for x in struc['ir_id']]].values)
    subsetyears = list(struc['created_year'].iloc[[x==i for x in struc['ir_id']]].values)
    sectionnotfound = []
    for j in range(len(subsetlines)):
        note = subsetlines[j]
        if 'Initial Consultation' in note:
            if 'GYNECOLOGICAL HISTORY:' in note:
                begin = note.split('GYNECOLOGICAL HISTORY:')[1]
                gyn_his = begin.split('REVIEW OF SYSTEMS:')[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Initial Consultation', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        elif 'Consultation' in note:
            # find the section of interest
            if 'Breast cancer risk factors:' in note:
                begin = re.split('Breast cancer risk factors:', note, 1)[1]
                gyn_his = re.split('FAMILY HISTORY|Family history:|Family History of Breast or Ovarian Cancer:|SOCIAL HISTORY:|Family Hx', begin)[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Consultation 1', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            elif 'Breast Cancer Risk Factors/Gynecologic History:' in note:
                begin = re.split('Breast Cancer Risk Factors/Gynecologic History:', note, 1)[1]
                gyn_his = begin.split('Family Cancer History:')[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Consultation 2', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        elif 'New Patient History and Physical' in note:
            if 'BREAST RISK DATA' in note:
                begin = note.split('BREAST RISK DATA')[1]
                gyn_his = re.split('Other Risk exposure:|Other R isk exposure|SOCIAL HISTORY:|Social History|Social Hx:', begin)[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'New Patient History and Physical', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        else:
            continue
    if np.all(sectionnotfound==0):
        print('Section of interest not found in patient %s' % i)


separated_sections.drop_duplicates(inplace=True, ignore_index=True)
separated_sections.to_csv(os.path.join(dn, 'initial_counsel/isolated_gynecological_history_section_initial_counsel.csv'), sep='!', header=True, index=False)

# init = separated_sections.iloc[[x=='Initial Consultation' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# con1 = separated_sections.iloc[[x=='Consultation 1' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# con2 = separated_sections.iloc[[x=='Consultation 2' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# newp = separated_sections.iloc[[x=='New Patient History and Physical' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
