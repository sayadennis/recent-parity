import re
import pickle
import numpy as np
import pandas as pd
from collections import defaultdict

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

#################################################
#### Identify patients who were missing HER2 ####
#################################################

datadir = 'data/05_import_to_redcap'
fn = '20220616_redcap_import_data.csv'
import_to_redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

correct_her2 = list(import_to_redcap.iloc[pd.isnull(import_to_redcap.her2_status).values,:].index)
# will split this list between patients with plain text note vs. not at the end 

############################################
#### Identify who is missing tumor char ####
############################################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-02-06_1107.csv'
redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# #### Getting text pattern to plug into SSMS query ####
# with open('/home/srd6051/2023-02-06_mrn_dob_sql_query_pattern_all.txt', 'w') as f:
#     for mrn_string, dob_string in zip(redcap.iloc[~pd.isnull(redcap.epic_mrn).values,:].epic_mrn.astype(int).astype(str), redcap.iloc[~pd.isnull(redcap.epic_mrn).values,:].dob):
#         f.write(f"OR (west_mrn LIKE '%{mrn_string[-6:]}') AND birth_date='{dob_string} 00:00:00.00'\n")

tum_char = ['histology', 'tumor_size', 'histologic_grade', 
            'tumor_staging_category', 'node_staging_category', 
            'er_status', 'pr_status', 'her2_status', 'ki67'] # , 'num_ln_positive'

missing = pd.DataFrame(index=redcap.index, columns=tum_char)
for item in tum_char:
    missing[item] = pd.isnull(redcap[item]).astype(int)

missing['any'] = [np.any(missing.loc[i]).astype(int) for i in missing.index]
missing.to_csv(f'{dn}/summary_tables/missing_tumor_char.csv')

########################################
#### Identify their pathology notes ####
########################################

pathol = pd.read_csv(
    f'{dn}/data/01_ssms_raw/pathology/20221213_cohort_pathol_reports.csv',
    # f'{dn}/data/01_ssms_raw/pathology/notes_pathology.csv',
    header=None
)

with open(f'{dn}/data/01_ssms_raw/pathology/colnames_20221213_cohort_pathol_reports.txt', 'r') as f:
    lines = f.readlines()

colnames = lines[0].split()
pathol.columns = colnames

# # rename some columns so that it's easier for Andrea to see 
# pathol_colnames = [
#     0,
#     'Epic_MRN',
#     'Powerchart_MRN',
#     'first_name',
#     'last_name',
#     'dob',
#     6,7,8,9,10,
#     'year',
#     12,13,
#     'note',
#     15,16,
#     'note_text'
# ]

# pathol.columns = pathol_colnames

pathol['dob'] = pd.to_datetime(pathol['birth_date'])
redcap['dob'] = pd.to_datetime(redcap['dob'])

patterns = [
	'surgical pathology final report',
	'surg path final report',
	'surgical pathology addendum report',
	'surg path addendum report',
	'tumor staging report',
	'breast cancer tumor markers',
	'breast cancer staging summary',
]

found_path = []

for pattern in patterns:
    for i in pathol.index:
        if re.search(pattern, pathol.loc[i,'pathol_note_text'], re.IGNORECASE):
            if i not in found_path:
                found_path.append(i)

# pathol = pathol.loc[found_path,:]

redcap_id_missing_but_no_pathol = []

correct_her2_dict = defaultdict(list)
for i in correct_her2:
    if not missing.loc[i,'any']:
        correct_her2_dict['all'].append(i)

for i in redcap.index:
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    first_name = redcap.loc[i,'first_name']
    last_name = redcap.loc[i,'last_name']
    # if there is a record in pathol that matches these... 
    matches = pathol.iloc[(
        (pathol['first_name'].values==first_name)
        & (pathol['last_name'].values==last_name) 
        & (pathol['dob'].values==dob)
    ),:]
    if matches.shape[0]>0:
        matches.to_csv(f'{dn}/chart_review_resources/pathology_notes_separated_by_patient_new/redcap_record_id_{i}.csv', index=False)
        if (i in correct_her2_dict['all']) & (i not in correct_her2_dict['has_pathol']):
            correct_her2_dict['has_pathol'].append(i)
    elif missing.loc[i,'any']:
        redcap_id_missing_but_no_pathol.append(i)
        if (i in correct_her2_dict['all']) & (i not in correct_her2_dict['no_pathol']):
            correct_her2_dict['no_pathol'].append(i)
    elif (i in correct_her2_dict['all']) & (i not in correct_her2_dict['no_pathol']):
            correct_her2_dict['no_pathol'].append(i)


with open(f'{dn}/chart_review_resources/redcap_id_missing_but_no_pathol_notes_new.txt', 'w') as f:
    for item in redcap_id_missing_but_no_pathol:
        f.write(f'{item}\n')

for key in correct_her2_dict.keys():
    with open(f'{dn}/chart_review_resources/her2_correction_redcap_id_{key}.txt', 'w') as f:
        for item in correct_her2_dict[key]:
            f.write(f'{item}\n')


###################################################
#### Create assistive files for family history ####
###################################################

initialcounsel = pd.read_csv(f'{dn}/data/01_ssms_raw/initialcounsel/notes_initial_counsel.csv', header=None)

with open(f'{dn}/data/01_ssms_raw/column_names_note_type_detail.txt', 'r') as f:
    lines = f.readlines()

colnames = [line.strip() for line in lines]

initialcounsel.columns = colnames
initialcounsel['dob'] = pd.to_datetime(initialcounsel['birth_date'])

for i in redcap.index:
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    first_name = redcap.loc[i,'first_name']
    last_name = redcap.loc[i,'last_name']
    # if there is a record in pathol that matches these... 
    matches = initialcounsel.iloc[(
        (initialcounsel['first_name'].values==first_name)
        & (initialcounsel['last_name'].values==last_name) 
        & (initialcounsel['dob'].values==dob)
    ),:]
    if matches.shape[0]>0:
        matches.to_csv(f'{dn}/chart_review_resources/initialcounsel_notes_separated_by_patient_new/redcap_record_id_{i}.csv', index=False)
