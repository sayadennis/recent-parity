import re
import pickle
import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Load data ####
###################

## Load notes ##
notesdir = 'data/01_ssms_raw/pathology'

# Read data for promising notes 
notes = pd.read_csv(f'{dn}/{notesdir}/20221109_cohort_pathol_reports.csv', header=None)

# Add column names 
with open(f'{dn}/{notesdir}/colnames_20221109_cohort_pathol_reports.txt', 'r') as f:
    lines = f.readlines()

colnames = [x.strip() for x in lines]
notes.columns = colnames

## Load REDCap data ##
redcapdir = 'data/06_exported_from_redcap'
redcap = pd.read_csv(f'{dn}/{redcapdir}/FrequencyAndResultsO_DATA_2022-10-20_1504.csv')

with open(f'{dn}/{redcapdir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

#######################
#### Add REDCap ID ####
#######################

notes['birth_date'] = [x.split()[0] for x in notes.birth_date.values]
notes['birth_date'] = pd.to_datetime(notes['birth_date'])
redcap['dob'] = pd.to_datetime(redcap['dob'])

##########################################
#### Assess presence of text patterns ####
##########################################

# Define staging summary text patterns 
patterns = [
	'surgical pathology final report',
	'surg path final report',
	'surgical pathology addendum report',
	'surg path addendum report',
	'tumor staging report',
	'breast cancer tumor markers',
	'breast cancer staging summary',
	'surgpath',
	'surgical pathology report',
	'synoptic report',
	'in-situ breast carcinoma checklist',
	'surgical pathology report',
]

# DataFrame to record number of occurrences 
presence = pd.DataFrame(
    None,
    columns=patterns,
    index=redcap.record_id,
)

missing = pd.read_csv(f'{dn}/summary_tables/missing_tumor_char.csv', index_col=0)

f_no_notes = open(f'{dn}/chart_review_resources/redcap_id_missing_but_no_pathol_notes.txt', 'w')
f_has_notes = open(f'{dn}/chart_review_resources/redcap_id_missing_but_has_pathol_notes.txt', 'w')

# record occurrences 
for i in redcap.index:
    if missing.loc[i,'any']==1:
        pattern_exists = []
        first_name = redcap.loc[i,'first_name']
        last_name = redcap.loc[i,'last_name']
        dob = redcap.loc[i,'dob']
        notes_matches = notes.iloc[(
            (notes['first_name'].values==first_name)
            & (notes['last_name'].values==last_name) 
            & (notes['birth_date'].values==dob)
        ),:]
        if notes_matches.shape[0]==0:
            f_no_notes.write(f'{i}\n')
            continue
        #
        f_has_notes.write(f'{i}\n')
        for pattern in patterns:
            matches = [len(re.findall(pattern, note_text, re.IGNORECASE))>0 for note_text in notes_matches.pathol_note_text.values]
            pattern_exists.append(int(np.any(matches)))
        presence.loc[i,:] = pattern_exists
        notes_matches.to_csv(f'{dn}/chart_review_resources/pathology_notes_separated_by_patient/redcap_record_id_{i}.csv')

print(((~pd.isnull(presence)).sum(axis=1)==0).sum(), 'patients have no notes.')
print(
    ((pd.isnull(presence).sum(axis=1)==0) & (presence.sum(axis=1)==0)).sum(), 
    'patients had some pathol notes but no known table title patterns.'
)

f_no_notes.close()
f_has_notes.close()

## How many of the patients with missing tumor characteristics have these notes available? 

# for ir_id in notes['ir_id'].unique():
#     pattern_exists = []
#     sub = notes.iloc[notes.ir_id.values==ir_id,:]

#     dob = sub['dob'].values[0]
#     first_name = sub['first_name'].values[0]
#     last_name = sub['last_name'].values[0]
#     # if there is a record in pathol that matches these... 
#     redcap_matches = redcap.iloc[(
#         (redcap['first_name'].values==first_name)
#         & (redcap['last_name'].values==last_name) 
#         & (redcap['dob'].values==dob)
#     ),:]
#     if redcap_matches.shape[0]>0:
#         redcap_id = redcap_matches.record_id.values[0]
#     else:
#         print('No match in REDCap for ir_id:', ir_id)
#         continue

#     for pattern in patterns:
#         matches = [len(re.findall(pattern, note_text, re.IGNORECASE))>0 for note_text in sub.note_text.values]
#         pattern_exists.append(np.any(matches))
#         sub_matches = sub.iloc[matches,:]
#         dep_type_pairs = sub_matches[['department_specialty','note_type']].drop_duplicates(ignore_index=True).itertuples(index=False)
#         for dep_spec, note_type in dep_type_pairs:
#             if f'{dep_spec}-{note_type}' in presence.columns:
#                 presence.loc[pattern,f'{dep_spec}-{note_type}'] += 1
#             else:
#                 presence = pd.concat((
#                     presence, 
#                     pd.DataFrame(
#                         0,
#                         index=presence.index,
#                         columns=[f'{dep_spec}-{note_type}'],
#                     )
#                 ), axis=1)
#                 presence.loc[pattern,f'{dep_spec}-{note_type}'] += 1
#     if not np.any(pattern_exists):
#         sub.to_csv(f'{dn}/chart_review_resources/any_notes_separated_by_patient/redcap_record_id_{redcap_id}.csv')

# presence.to_csv('text_pattern_presence.csv')
