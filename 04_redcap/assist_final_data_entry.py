import re
import pickle
import numpy as np
import pandas as pd

############################################
#### Identify who is missing tumor char ####
############################################

dn = '/share/fsmresfiles/breast_cancer_pregnancy'
datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-20_1504.csv'
redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

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
    f'{dn}/data/01_ssms_raw/pathology/notes_pathology.csv',
    header=None
)

# rename some columns so that it's easier for Andrea to see 
pathol_colnames = [
    0,
    'Epic_MRN',
    'Powerchart_MRN',
    'first_name',
    'last_name',
    'dob',
    6,7,8,9,10,
    'year',
    12,13,
    'note',
    15,16,
    'note_text'
]

pathol.columns = pathol_colnames

pathol['dob'] = pd.to_datetime(pathol['dob'])
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
        if re.search(pattern, pathol.loc[i,'note_text'], re.IGNORECASE):
            if i not in found_path:
                found_path.append(i)

pathol = pathol.loc[found_path,:]

redcap_id_missing_but_no_pathol = []

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
        matches.to_csv(f'{dn}/chart_review_resources/pathology_notes_separated_by_patient/redcap_record_id_{i}.csv', index=False)
    elif missing.loc[i,'any']:
            redcap_id_missing_but_no_pathol.append(i)

with open(f'{dn}/chart_review_resources/redcap_id_missing_but_no_pathol_notes.txt', 'w') as f:
    for item in redcap_id_missing_but_no_pathol:
        f.write(f'{item}\n')
