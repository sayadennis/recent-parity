import pickle
import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

##########################
#### Load REDCap data ####
##########################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-02-13_1219.csv'
redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

redcap.dob = pd.to_datetime(redcap.dob)

# #### Getting text pattern to plug into SSMS query ####
# with open('/home/srd6051/2023-02-10_mrn_dob_sql_query_pattern_all.txt', 'w') as f:
#     for mrn_string, dob_string in zip(redcap.iloc[~pd.isnull(redcap.epic_mrn).values,:].epic_mrn.astype(int).astype(str), redcap.iloc[~pd.isnull(redcap.epic_mrn).values,:].dob):
#         f.write(f"OR ((west_mrn LIKE '%{mrn_string[-6:]}') AND birth_date='{dob_string} 00:00:00.00')\n")

#################################
#### Load SSMS Query results ####
#################################

datadir = 'data/exploratory'
fn = '2023-02-06_mrn_dob_pathol_query_results.txt'
ssms = pd.read_csv(f'{dn}/{datadir}/{fn}', sep='\t')
ssms.birth_date = pd.to_datetime(ssms.birth_date)

#######################################################################################
#### Get list of REDCap IDs of patients that exist in REDCap but not in SSMS table ####
#######################################################################################

haspath_redcap_id = []
nopath_redcap_id = []
for i in redcap.index:
    redcap_id = redcap.loc[i,'record_id']
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    if ssms.iloc[(ssms.west_mrn.values==epic_mrn) & (ssms.birth_date.values==dob)].shape[0]==0:
        nopath_redcap_id.append(redcap_id)
    else:
        haspath_redcap_id.append(redcap_id)

# np.random.seed(11)
# random_subset = np.random.choice(nopath_redcap_id, size=20)
# random_subset = np.sort(random_subset)
with open('/home/srd6051/recent_parity/2023-02-08_redcap_id_without_pathology_reports_in_edw.txt', 'w') as f:
    for redcap_id in nopath_redcap_id:
        if redcap.iloc[redcap.record_id.values==redcap_id,:]['tumor_characteristics_complete'].values[0]!=2:
            f.write(f'{redcap_id}\n')

#######################################################################################
#### How many more patients are missing tumor char and also missing path reports?? ####
#######################################################################################

tumchar_incomplete = redcap.iloc[
    redcap.tumor_characteristics_complete.map({0: 'Incomplete', 1: 'Unverified', 2: 'Complete'}).values!='Complete',:
]

upset_data = pd.DataFrame(
    index=redcap.record_id.values, 
    columns=[
        'Tumor char incomplete',
        'NAT incomplete',
        'Fam Hx incomplete',
        'Missing MRN',
        'Cannot find pathol on EDW',
    ]
)

upset_data['Tumor char incomplete'] = (redcap.tumor_characteristics_complete==0).astype(int).values
upset_data['NAT incomplete'] = (redcap.treatment_complete!=2).astype(int).values
upset_data['Fam Hx incomplete'] = (redcap.geneticsfam_hx_complete!=2).astype(int).values
upset_data['Missing MRN'] = pd.isnull(redcap.epic_mrn).astype(int).values

upset_data.loc[nopath_redcap_id,'Cannot find pathol on EDW'] = 1
upset_data['Cannot find pathol on EDW'] = upset_data['Cannot find pathol on EDW'].fillna(0)

###############################################
#### Prepare data for Takahiro's RCB entry ####
###############################################

pathol = ssms

with open(f'{dn}/data/01_ssms_raw/pathology/colnames_20221213_cohort_pathol_reports.txt', 'r') as f:
    lines = f.readlines()

colnames = lines[0].split()
pathol.columns = colnames

for i in redcap.index:
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    first_name = redcap.loc[i,'first_name']
    last_name = redcap.loc[i,'last_name']
    # if there is a record in pathol that matches these... 
    matches = pathol.iloc[(
        (pathol['first_name'].values==first_name)
        & (pathol['last_name'].values==last_name) 
        & (pathol['birth_date'].values==dob)
    ),:]
    if matches.shape[0]>0:
        matches.to_csv(f'{dn}/chart_review_resources/pathology_notes_separated_by_patient_new/redcap_record_id_{i}.csv', index=False)

# Get a list of patients who has pathol records but are missing RCB data 
missing_rcb = redcap.iloc[
    (redcap.nat.values==1) & (pd.isnull(redcap.rcb).values  | pd.isnull(redcap.rcb_category).values),:
]

takahiro = missing_rcb.iloc[[x in haspath_redcap_id for x in missing_rcb.record_id],:].record_id.values

with open('/home/srd6051/recent_parity/2023-02-08_redcap_id_needs_rcb_from_fsm_pathol.txt', 'w') as f:
    for item in takahiro:
        f.write(f'{item}\n')


###################################################
#### Create assistive files for family history ####
###################################################

initialcounsel = pd.read_csv(f'{dn}/data/01_ssms_raw/initialcounsel/2023-02-10-cohort_initial_counsel_notes.csv', header=None)

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
