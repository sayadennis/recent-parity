import os
import numpy as np
import pandas as pd
import re

din='/share/fsmresfiles/breast_cancer_pregnancy/explore_cohort'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

fn_colname='/share/fsmresfiles/breast_cancer_pregnancy/data/demo_and_dx/column_names_demographic.txt'

colnames=[]
with open(fn_colname, 'r') as f:
    lines = f.readlines()

for line in lines:
    colnames.append(line.strip())

cohort_dict = {
    'BC diagnosis' : pd.read_csv(f'{din}/cohort_by_breastcanceronly.csv', header=None, names=colnames),
    'Provider names' : pd.read_csv(f'{din}/cohort_by_prov_name.csv', header=None, names=colnames),
    'Billing' : pd.read_csv(f'{din}/cohort_by_billing_codes.csv', header=None, names=colnames),
    'Department' : pd.read_csv(f'{din}/cohort_by_encounter_department.csv', header=None, names=colnames),
    'Notes' : pd.read_csv(f'{din}/cohort_by_genetic_notes.csv', header=None, names=colnames)
}


cohort_dict = {
    'BC diagnosis' : list(pd.read_csv(f'{din}/cohort_by_breastcanceronly.csv', header=None, names=colnames)['ir_id'].unique()),
    'Provider names' : list(pd.read_csv(f'{din}/cohort_by_prov_name.csv', header=None, names=colnames)['ir_id'].unique()),
    'Billing' : list(pd.read_csv(f'{din}/cohort_by_billing_codes.csv', header=None, names=colnames)['ir_id'].unique()),
    'Department' : list(pd.read_csv(f'{din}/cohort_by_encounter_department.csv', header=None, names=colnames)['ir_id'].unique()),
    'Notes' : list(pd.read_csv(f'{din}/cohort_by_genetic_notes.csv', header=None, names=colnames)['ir_id'].unique())
}

co = pd.DataFrame(index=None, columns=['BC diagnosis', 'Provider names', 'Billing', 'Department', 'Notes'])

for key in cohort_dict.keys():
    for ir_id in cohort_dict[key]:
        if ir_id in co.index:
            co.loc[ir_id, key] = 1
        else:
            co.loc[ir_id,:] = 0
            co.loc[ir_id, key] = 1

co.to_csv(f'{dout}/cooccurrence_recent_parity_cohorts.csv')

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

set1 = set(cohort_dict['Billing'])
# set2 = set(cohort_dict['Department'])
set2 = set(cohort_dict['Provider names'])
set3 = set(cohort_dict['Notes'])

venn3([set1, set2, set3], ('Billing', 'Provider names', 'Genetic Counseling Notes'))
plt.savefig(f'{dout}/venn_recent_parity_cohorts.png')
plt.close()

#####################################################
#### Looking at overlaps including note contents ####
#####################################################

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/notes'

with open(os.path.join(dn, 'column_names_note_type_detail.txt'), 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())

gennote = pd.read_csv(f'{din}/exploratory_notes_genetic.csv', header=None, names=colnames)
prov = pd.read_csv(f'{din}/exploratory_notes_prov_name.csv', header=None, names=colnames)

ol = 0
for ir_id in gennote['ir_id'].unique():
    if ir_id in prov['ir_id']:
        ol += 1

## Plot Venn diagram showing overlap between the two ##
from matplotlib_venn import venn2
venn2((len(gennote['ir_id'].unique), len(prov['ir_id'].unique()), ol), set_labels=('Selected by genetic notes', 'Selected by provider encounters'))
plt.savefig(f'{dout}/venn_gennotes_prov_overlap.png')
plt.close()

## Plot how many of those have genetic results ## 

## Fist for genetic notes dataset 
separated_sections = pd.DataFrame(None, columns=['ir_id', 'test_results']) # , 'interpretation' 'created_year'

## Loop through patient ir_id's and isolate sections 
for i in gennote['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(gennote['note_text'].iloc[[x==i for x in gennote['ir_id']]].values)
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


## Next for provider-selected dataset 

separated_sections = pd.DataFrame(None, columns=['ir_id', 'test_results']) # , 'interpretation' 'created_year'

## Loop through patient ir_id's and isolate sections 
for i in prov['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(prov['note_text'].iloc[[x==i for x in prov['ir_id']]].values)
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


