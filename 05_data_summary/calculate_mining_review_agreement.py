import numpy as np
import pandas as pd

reviewed = pd.read_excel(f'/share/fsmresfiles/breast_cancer_pregnancy/data/03_reviewed_erica_takahiro/takahiro_erica_reviewed_200.xlsx')
mined = pd.read_csv(f'/share/fsmresfiles/breast_cancer_pregnancy/data/05_import_to_redcap/20220616_redcap_import_data.csv')

reviewed['DOB'] = pd.to_datetime(reviewed['DOB'])
mined['dob'] = pd.to_datetime(mined['dob'])

def match_patient(reviewed_sample, mined_df):
    # get reviewed data 
    reviewed_dob = reviewed_sample['DOB']
    reviewed_mrn = reviewed_sample['MRN']
    reviewed_last_name = reviewed_sample['Patient Name'].split()[-1]
    # match to mined 
    matched = mined_df.iloc[
        (reviewed_last_name==np.array([x.split()[-1] for x in mined_df['last_name'].values]))
        & (reviewed_mrn==mined_df['epic_mrn'].values)
        & (reviewed_dob==mined_df['dob']).values,
        :
    ]
    return matched

demo = {
    'Age at Breast CA Dx' : 'age_at_diagnosis',
}

gyn = {
    'Number Pregnancies (G)' : 'number_pregnancies', 
    '# Live Births (P)' : 'number_births', 
    'Age at first pregnancy' : 'age_at_first_pregnancy',
    'Age at most recent pregnancy' : 'age_at_most_recent_pregnancy',
}

tumor = {
    'Histology (IDC, ILC, DCIS, other)' : 'histology',
    'Tumor size (mm)' : 'tumor_size', 
    'Histologic grade (1-3)' : 'histologic_grade',
    'ER status (pos/neg)' : 'er_status',
    'PR status (pos/neg)' : 'pr_status', 
    'HER2Neu status (pos/neg/not performed)' : 'her2_status',
    'Ki67 level (low/int/high/not performed)' : 'ki67',
}

# genetic = {
#     'Genetic Testing Results' : '', # UNFORMATTED
# }

loss = pd.DataFrame(
    index=reviewed.index,
    columns=list(demo.values()) + list(gyn.values()) + list(tumor.values())
)

missing = []

for i in loss.index:
    reviewed_sample = reviewed.loc[i,:]
    matched_mined = match_patient(reviewed_sample, mined)
    if matched_mined.shape[0]>0:
        # demo 
        key = 'Age at Breast CA Dx'
        if ',' in str(reviewed_sample.loc[key]):
            loss.loc[i,demo[key]] = np.min([int(x) for x in reviewed_sample.loc[key].split(', ')]) - matched_mined.loc[:,demo[key]].values[0]
        elif pd.isnull(reviewed_sample.loc[key]):
            continue
        else:
            loss.loc[i,demo[key]] = reviewed_sample.loc[key]-matched_mined.loc[:,demo[key]].values[0]
        # gynecological 
        for key in gyn.keys():
            loss.loc[i,gyn[key]] = reviewed_sample.loc[key]-matched_mined.loc[:,gyn[key]].values[0]
        # tumor size
        key = 'Tumor size (mm)'
        if str(reviewed_sample.loc[key]).startswith(', '):
            reviewed_sample.loc[key] = int(reviewed_sample.loc[key][2:])
        elif ',' in str(reviewed_sample.loc[key]):
            reviewed_sample.loc[key] = np.max([int(x) for x in reviewed_sample.loc[key].split(',')])
            loss.loc[i,tumor[key]] = reviewed_sample.loc[key]/10 - matched_mined.loc[:,tumor[key]].values[0]
        elif '/' in str(reviewed_sample.loc[key]):
            reviewed_sample.loc[key] = np.max([int(x) for x in reviewed_sample.loc[key].split('/')])
            loss.loc[i,tumor[key]] = reviewed_sample.loc[key]/10 - matched_mined.loc[:,tumor[key]].values[0]
        elif type(reviewed_sample.loc[key])==str:
            continue
        else:
            loss.loc[i,tumor[key]] = reviewed_sample.loc[key]/10 - matched_mined.loc[:,tumor[key]].values[0]
        # histology 
        key = 'Histology (IDC, ILC, DCIS, other)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_histology = {1:'IDC', 2:'ILC', 3:'DCIS', 4:'Other'}[matched_mined.loc[:,tumor[key]].values[0]]
            if reviewed_sample.loc[key]=='Mixed':
                loss.loc[i,tumor[key]] = 0
            else:
                loss.loc[i,tumor[key]] = 0 if matched_histology.lower() in reviewed_sample.loc[key].lower() else 1
        # grade 
        key = 'Histologic grade (1-3)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_grade = {1:'No histologic grade (DCIS)', 2:'1', 3:'2', 4:'3'}[matched_mined.loc[:,tumor[key]].values[0]]
            if 'DCIS' in reviewed_sample.loc['Histology (IDC, ILC, DCIS, other)']:
                loss.loc[i,tumor[key]] = 0
            else:
                loss.loc[i,tumor[key]] = 0 if matched_grade.lower() in str(reviewed_sample.loc[key]) else 1
        # er
        key = 'ER status (pos/neg)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_er = {1:'Positive', 2:'Negative'}[matched_mined.loc[:,tumor[key]].values[0]].lower()
            loss.loc[i,tumor[key]] = 0 if str(reviewed_sample.loc[key]).lower() in matched_er else 1
        # pr
        key = 'PR status (pos/neg)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_pr = {1:'Positive', 2:'Negative'}[matched_mined.loc[:,tumor[key]].values[0]].lower()
            loss.loc[i,tumor[key]] = 0 if str(reviewed_sample.loc[key]).lower() in matched_pr else 1
        # her2 
        key = 'HER2Neu status (pos/neg/not performed)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_her2 = {1:'Positive', 2:'Negative', 3:'Not performed'}[matched_mined.loc[:,tumor[key]].values[0]].lower()
            loss.loc[i,tumor[key]] = 0 if str(reviewed_sample.loc[key]).lower()[-3:] in matched_her2 else 1
        # ki67
        key = 'Ki67 level (low/int/high/not performed)'
        if pd.isnull(matched_mined.loc[:,tumor[key]].values[0]):
            continue
        else:
            matched_ki67 = {1:'Low', 2:'Intermediate', 3:'High', 4:'Not performed'}[matched_mined.loc[:,tumor[key]].values[0]].lower()
            if ',' in str(reviewed_sample.loc[key]):
                loss.loc[i,tumor[key]] = 0 if np.any([x.strip() in matched_ki67 for x in str(reviewed_sample.loc[key]).lower().split(',')]) else 1
            elif '/' in str(reviewed_sample.loc[key]):
                loss.loc[i,tumor[key]] = 0 if np.any([x.strip() in matched_ki67 for x in str(reviewed_sample.loc[key]).lower().split('/')]) else 1
            else:
                loss.loc[i,tumor[key]] = 0 if str(reviewed_sample.loc[key]).lower() in matched_ki67 else 1
    else:
        missing.append(i)

loss.to_csv('/home/srd6051/20230401_patient_agreement_test.csv')

np.mean(loss**2, axis=0)
