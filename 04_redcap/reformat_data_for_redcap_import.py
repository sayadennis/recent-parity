import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

proj='/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Load data ####
###################

## REDCap features template 
redcap_features=list(pd.read_csv(f'{proj}/data/06_exported_from_redcap/FrequencyAndResultsOfGeneticTe_ImportTemplate_2022-05-31.csv').columns)
redcap_features.remove('Unnamed: 60')

## Import reviewed data (missing tumor characteristics)
rev_comb=pd.read_csv(f'{proj}/data/03_reviewed_erica_takahiro/tkhr_pat_wo_gentest_added/recent_parity_combined_data_reviewed_tkhr_added.csv')

## Import tumor characteristics 
singles=pd.read_csv(f'{proj}/data/exploratory/singles.csv')
notable=pd.read_csv(f'{proj}/data/exploratory/no_table.csv')

#########################################################
#### Fill in tumor characteristics in reviewed table ####
#########################################################

rev_comb['pTNM_T']=None; rev_comb['pTNM_N']=None
rev_comb['cTNM_T']=None; rev_comb['cTNM_N']=None

tumchar_map={
    'Histology' : 'cancer histology',
    'PathTumorSize' : 'tumor size',
    'PathGrade' : 'cancer pathologic sbr grade',
    'pTNM_T' : 'pathological tumor staging category',
    'pTNM_N' : 'pathological nodes staging category',
    'ER' : 'estrogen receptor status', 
    'PR' : 'progesterone receptor status', 
    'HER2' : 'her2 status', 
    'Ki-67' : 'ki67',
    'LymphInv' : 'lymphovascular invasion', 
    'NumLymphExamined' : 'number lymph nodes examined', 
    'NumPosLymph' : 'number lymph nodes positive tumor',
}

for i in rev_comb.index:
    if pd.isnull(rev_comb.loc[i,'EPIC_mrn']):
        first_name=rev_comb.loc[i,'FirstName']
        last_name=rev_comb.loc[i,'LastName']
        if ((singles['first_name']==first_name) & (singles['last_name']==last_name)).sum()==1:
            case=singles.iloc[((singles['first_name']==first_name) & (singles['last_name']==last_name)).values,:]
            for key in tumchar_map.keys():
                if pd.isnull(rev_comb.loc[i,key]):
                    rev_comb.loc[i,key]=case.loc[:,tumchar_map[key]].values[0]
        elif ((notable['first_name']==first_name) & (notable['last_name']==last_name)).sum()>0:
            cases=singles.iloc[((notable['first_name']==first_name) & (notable['last_name']==last_name)).values,:]
            for key in tumchar_map.keys():
                if pd.isnull(rev_comb.loc[i,key]):
                    if (~pd.isnull(cases[tumchar_map[key]])).sum()==1: # if one non-NaN value 
                        rev_comb.loc[i,key]=cases.iloc[~pd.isnull(cases[tumchar_map[key]]).values,:][tumchar_map[key]].values[0]
    else:
        mrn=rev_comb.loc[i,'EPIC_mrn']
        # look in singles 
        if mrn in singles.EPIC_mrn.values:
            case=singles.iloc[singles.EPIC_mrn.values==mrn,:]
            for key in tumchar_map.keys():
                if pd.isnull(rev_comb.loc[i,key]):
                    rev_comb.loc[i,key]=case.loc[:,tumchar_map[key]].values[0]
        elif mrn in notable.EPIC_mrn.values:
            cases=notable.iloc[notable.EPIC_mrn.values==mrn,:]
            for key in tumchar_map.keys():
                if pd.isnull(rev_comb.loc[i,key]):
                    if (~pd.isnull(cases[tumchar_map[key]])).sum()==1: # if one non-NaN value 
                        rev_comb.loc[i,key]=cases.iloc[~pd.isnull(cases[tumchar_map[key]]).values,:][tumchar_map[key]].values[0]

#####################################################
#### Split up staging categories to tumor + node ####
#####################################################

for i in rev_comb.index:
    if not pd.isnull(rev_comb.loc[i,'cTNM']): # if not none
        # tumor staging 
        if 'T0' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='T0'
        elif 'Tis' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='Tis'
        elif 'T1' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='T1'
        elif 'T2' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='T2'
        elif 'T3' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='T3'
        elif 'T4' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_T']='T4'
        # node staging 
        if (('Nx' in rev_comb.loc[i,'cTNM']) | ('NX' in rev_comb.loc[i,'cTNM'])):
            rev_comb.loc[i,'cTNM_N']='NX'
        elif 'N0' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_N']='N0'
        elif 'N1' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_N']='N1'
        elif 'N2' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_N']='N2'
        elif 'N3' in rev_comb.loc[i,'cTNM']:
            rev_comb.loc[i,'cTNM_N']='N3'


for i in rev_comb.index:
    if not pd.isnull(rev_comb.loc[i,'pTNM']): # if not none
        # tumor staging 
        if 'T0' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pT0'
        elif 'Tis' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pTis'
        elif 'T1' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pT1'
        elif 'T2' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pT2'
        elif 'T3' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pT3'
        elif 'T4' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_T']='pT4'
        # node staging 
        if (('Nx' in rev_comb.loc[i,'pTNM']) | ('NX' in rev_comb.loc[i,'pTNM'])):
            rev_comb.loc[i,'pTNM_N']='pNX'
        elif 'N0' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_N']='pN0'
        elif 'N1' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_N']='pN1'
        elif 'N2' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_N']='pN2'
        elif 'N3' in rev_comb.loc[i,'pTNM']:
            rev_comb.loc[i,'pTNM_N']='pN3'

# # Way to examine missing rate 
# print((pd.isnull(rev_comb).sum()/rev_comb.shape[0]).to_string())

########################
#### Rename columns ####
########################

feature_map={
    #### Demographics etc. ####
    'EPIC_mrn' : 'epic_mrn', 
    'FirstName' : 'first_name',
    'LastName' : 'last_name',
    'DOB' : 'dob',
    'AgeBCDx' : 'age_at_diagnosis', 
    'TimeBCDx' : 'year_of_diagnosis',
    #### Gyn history ####
    'Menarche' : 'menarche',
    'Gravida' : 'number_pregnancies', 
    'Para' : 'number_births', 
    'AgeFirstPreg' : 'age_at_first_pregnancy',
    'AgeLastPreg' : 'age_at_most_recent_pregnancy', 
    'DurationBF' : 'breastfeeding_duration',
    'LastNursed' : 'recency_of_lactation', 
    #### Tumor characteristics ####
    'Surgery' : 'surgery_type', 
    'Histology' : 'histology',
    'PathTumorSize' : 'tumor_size',
    'PathGrade' : 'histologic_grade',
    'NumLymphExamined' : 'num_lymph_nodes_taken', # DTYPE MIGHT NEED FIXING 
    'NumPosLymph' : 'num_ln_positive', 
    'ER' : 'er_status', 
    'PR' : 'pr_status', 
    'HER2' : 'her2_status',
    'Ki-67' : 'ki67', 
    'pTNM_T' : 'tumor_staging_category',
    'pTNM_N' : 'node_staging_category',
    #### Treatment ####
    'NACT' : 'nat', 
    'NACTReg' : 'nat_reg', 
    'HER2Therapy' : 'her2_therapy',
    'HER2TherapyType' : 'her2_therapy_type',
    'ClinTumorSize' : 'clin_tumor_size', 
    'ClinTumorSizeMeasure' : 'measure_clin_size',
    'CNBGrade' : 'cnb_grade', 
    'cTNM_T' : 'clin_tumor_stag_cat', 
    'cTNM_N' : 'clin_node_stag_cat', 
    'Response' : 'response_nat', 
    #### Genetic testing ####
    'HasPathoMut' : 'any_patho_mutation',
    'HasVUS' : 'any_vus_mutation', 
    # 'BRCA1_patho' : 'brca1_patho', 
    # 'BRCA2_patho' : 'brca2_patho', 
    # 'PALB2_patho' : 'palb2_patho', 
    # 'TP53_patho' : 'tp53_patho', 
    # 'PTEN_patho' : 'pten_patho',
    # 'CDH1_patho' : 'cdh1_patho', 
    # 'STK11_patho' : 'stk11_patho', 
    # 'CHEK2_patho' : 'chek2_patho', 
    # 'ATM_patho' : 'atm_patho', 
    # 'BRCA1_VUS' : 'brca1_vus', 
    # 'BRCA2_VUS' : 'brca2_vus', 
    # 'PALB2_VUS' : 'palb2_vus', 
    # 'TP53_VUS' : 'tp53_vus', 
    # 'PTEN_VUS' : 'pten_vus',
    # 'CDH1_VUS' : 'cdh1_vus', 
    # 'STK11_VUS' : 'stk11_vus', 
    # 'CHEK2_VUS' : 'chek2_vus', 
    # 'ATM_VUS' : 'atm_vus', 
}

data = rev_comb.rename(feature_map, axis=1)

###############################################################
#### Reformat data to match dtype with REDCap restrictions ####
###############################################################

## Go through features and reformat 
# 1 - Demographics 
# 1a - MRN
data.epic_mrn = data.epic_mrn.astype('Int64')

# 1b - DOB
def mdy_to_ymd(mdy_string):
    mn, dy, yr = mdy_string.split('/')
    if int(yr)>20:
        yr='19'+yr
    else:
        yr='20'+yr
    #
    if int(mn)<10:
        mn='0'+mn
    #
    if int(dy)<10:
        dy='0'+dy
    return yr+'-'+mn+'-'+dy

for i in data.index:
    if '/' in data.loc[i,'dob']:
        data.loc[i,'dob']=mdy_to_ymd(data.loc[i,'dob'])

# 2 - Gyn history 
# 2a - menarche 
for i in data.index:
    try:
        data.loc[i,'menarche']=int(data.loc[i,'menarche'])
    except:
        try:
            data.loc[i,'menarche']=int(float(data.loc[i,'menarche']))
        except:
            numfound=re.findall(str(data.loc[i,'menarche']), r'[0-9]+')
            if len(numfound)>0:
                data.loc[i,'menarche']=int(numfound[0])
            else:
                data.loc[i,'menarche']=None

# 2b - gravida, para, ages 
data.number_births=data.number_births.astype('Int64')
data.number_pregnancies=data.number_pregnancies.astype('Int64')
data.age_at_first_pregnancy=data.age_at_first_pregnancy.astype('Int64')
data.age_at_most_recent_pregnancy=data.age_at_most_recent_pregnancy.astype('Int64')

data['parous']=(data.number_births.astype('Int64')>0).astype('Int64')
data['years_since_pregnancy']=(data.age_at_diagnosis.astype(float)-data.age_at_most_recent_pregnancy.astype(float))
for i in data.index:
    if data.loc[i,'years_since_pregnancy']<0:
        data.loc[i,'years_since_pregnancy']=None

data.years_since_pregnancy=data.years_since_pregnancy.astype('Int64')

# 2c - lactation related
for i in data.index:
    dur=re.findall(r'[0-9]+', str(data.loc[i,'breastfeeding_duration']))
    if len(dur)>0:
        data.loc[i,'breastfeeding_duration']=int(dur[0])
    else:
        data.loc[i,'breastfeeding_duration']=None
    #
    rec=re.findall(r'[0-9]+', str(data.loc[i,'recency_of_lactation']))
    if len(rec)>0:
        data.loc[i,'recency_of_lactation']=int(rec[0])
    else:
        data.loc[i,'recency_of_lactation']=None

data['breastfed']=None
for i in data.index:
    if not pd.isnull(data.loc[i,'breastfeeding_duration']):
        data.loc[i,'breastfed']='Yes'

# 3 - Genetic testing data 
for i in data.index:
    if np.any(data.loc[i,['brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']]=='Pathogenic'):
        data.loc[i,'any_patho_mutation']='Present'
    else:
        data.loc[i,'any_patho_mutation']='Absent'
    if np.any(data.loc[i,['brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']]=='VUS'):
        data.loc[i,'any_vus_mutation']='Present'
    else:
        data.loc[i,'any_vus_mutation']='Absent'

# ## Fill empty fields with "Not performed"
# for i in data.index:
#     for gene_col in ['brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']:
#         if pd.isnull(data.loc[i,gene_col]):
#             data.loc[i,gene_col]='Not performed'


# 4 - Tumor characteristics 
# data.histologic_grade=data.histologic_grade.astype(str) # converting to string since we are adding 'No histologic grade (DCIS)' as an option 

for i in data.index:
    # 4a - histology 
    if pd.isnull(data.loc[i,'histology']):
        data.loc[i,'histology']=None
    elif 'IDC' in str(data.loc[i,'histology']):
        data.loc[i,'histology']='IDC'
    elif 'ILC' in str(data.loc[i,'histology']):
        data.loc[i,'histology']='ILC'
    elif 'DCIS' in str(data.loc[i,'histology']):
        data.loc[i,'histology']='DCIS'
        data.loc[i,'histologic_grade']='No histologic grade (DCIS)'
    else:
        data.loc[i,'histology']='Other'
    # 4b - tumor size 
    if pd.isnull(data.loc[i,'tumor_size']):
        data.loc[i,'tumor_size']=None
    elif data.loc[i,'tumor_size']=='.':
        data.loc[i,'tumor_size']=None
    elif ((type(data.loc[i,'tumor_size'])==np.float64) | (type(data.loc[i,'tumor_size'])==float)):
        data.loc[i,'tumor_size']=data.loc[i,'tumor_size']
    elif type(data.loc[i,'tumor_size'])==int:
        data.loc[i,'tumor_size']=float(data.loc[i,'tumor_size'])
    elif type(data.loc[i,'tumor_size'])==str:
        num_ls=re.findall(r'[0-9\.]+', data.loc[i,'tumor_size'])
        data.loc[i,'tumor_size']=max([float(x) for x in num_ls])
    else:
        print(f'Unknown pattern at index {i}: ', data.loc[i,'tumor_size'])
    # 4c - number of lymph nodes taken 
    if len(re.findall(r'[0-9]+', str(data.loc[i,'num_lymph_nodes_taken'])))>0:
        data.loc[i,'num_lymph_nodes_taken']=np.max([int(x) for x in re.findall(r'[0-9]+', str(data.loc[i,'num_lymph_nodes_taken']))])
    else:
        data.loc[i,'num_lymph_nodes_taken']=None
    # 4d - number of lymph nodes positive 
    if len(re.findall(r'[0-9]+', str(data.loc[i,'num_ln_positive'])))>0:
        data.loc[i,'num_ln_positive']=np.max([int(x) for x in re.findall(r'[0-9]+', str(data.loc[i,'num_ln_positive']))])
    else:
        data.loc[i,'num_ln_positive']=None
    # 4e - ER and PR
    for colname in ['er_status', 'pr_status']:
        if ((pd.isnull(data.loc[i,colname])) | (data.loc[i,colname]=='.')):
            data.loc[i,colname]=None
        elif (('positive' in str(data.loc[i,colname]).lower()) | ('+' in str(data.loc[i,colname]).lower())):
            data.loc[i,colname]='Positive'
        elif (('negative' in str(data.loc[i,colname]).lower()) | ('-' in str(data.loc[i,colname]).lower())):
            data.loc[i,colname]='Negative'
        elif len(re.findall(r'[0-9]+', str(data.loc[i,colname])))>0:
            if np.max([int(x) for x in re.findall(r'[0-9]+', str(data.loc[i,colname]))])>0:
                data.loc[i,colname]='Positive'
            else:
                data.loc[i,colname]='Negative'
        else:
            print(f'Unknown pattern at index {i} for {colname}: ', str(data.loc[i,colname]))
    # 4f - HER2
    if ((pd.isnull(data.loc[i,'her2_status'])) | (data.loc[i,'her2_status']=='.')):
        data.loc[i,'her2_status']=None
    elif (('positive' in str(data.loc[i,'her2_status']).lower()) | ('+' in str(data.loc[i,'her2_status']).lower())):
        data.loc[i,'her2_status']='Positive'
    elif (('negative' in str(data.loc[i,'her2_status']).lower()) | ('-' in str(data.loc[i,'her2_status']).lower()) | ('0' in str(data.loc[i,'her2_status']).lower())):
        data.loc[i,'her2_status']='Negative'
    elif len(re.findall(r'[0-9]+', str(data.loc[i,'her2_status'])))>0:
        if np.max([int(x) for x in re.findall(r'[0-9]+', str(data.loc[i,'her2_status']))])>0:
            data.loc[i,'her2_status']='Positive'
        else:
            data.loc[i,'her2_status']='Negative'
    else:
        print(f'Unknown pattern at index {i} for her2_status: ', str(data.loc[i,'her2_status']))
    # 4g - Ki-67
    if pd.isnull(data.loc[i,'ki67']):
        data.loc[i,'ki67']=None
    elif type(data.loc[i,'ki67'])==str:
        if len(re.findall(r'[0-9]+', data.loc[i,'ki67']))>0:
            num=np.max([float(x) for x in re.findall(r'[0-9]+', data.loc[i,'ki67'])])
            if num>=20:
                data.loc[i,'ki67']='High'
            elif num>=10:
                data.loc[i,'ki67']='Intermediate'
            else:
                data.loc[i,'ki67']='Low'
        else:
            data.loc[i,'ki67']='Not performed'
    else: # if not string or NA
        num=data.loc[i,'ki67']
        if num>=20:
            data.loc[i,'ki67']='High'
        elif num>=10:
            data.loc[i,'ki67']='Intermediate'
        else:
            data.loc[i,'ki67']='Low'
    # 4h - tumor staging category 
    if not pd.isnull(data.loc[i,'tumor_staging_category']): # if not none
        # tumor staging 
        for stage in ['ypT4','ypT3','ypT2','ypT1mi','ypT1','ypTis','ypT0','ypTx','pT4','pT3','pT2','pT1mi','pT1','pTis','pT0','pTx']:
            if stage in str(data.loc[i,'tumor_staging_category']):
                data.loc[i,'tumor_staging_category']=stage
                break
    # 4i - node staging category 
        # node staging 
    if not pd.isnull(data.loc[i,'node_staging_category']):
        for stage in ['ypN3','ypN2','ypN1mi','ypN1','ypN0','ypNX','pN3','pN2','pN1mi','pN1','pN0','pNX']:
            if stage in str(data.loc[i,'node_staging_category']):
                data.loc[i,'node_staging_category']=stage
                break
    # 4j - tumor grade
    if pd.isnull(data.loc[i,'histologic_grade']):
        data.loc[i,'histologic_grade']=None
    elif type(data.loc[i,'histologic_grade'])==str:
        if len(re.findall(r'[0-9]', data.loc[i,'histologic_grade']))>0:
            num=np.max([float(x) for x in re.findall(r'[0-9]', data.loc[i,'histologic_grade'])])
            data.loc[i,'histologic_grade']=num
        elif data.loc[i,'histologic_grade']=='.':
            data.loc[i,'histologic_grade']==None

# 5 - Treatment info
for i in data.index:
    # 5a - nat 
    # first evaluate based on what's already entered in the NAT column 
    if ((data.loc[i,'nat']==1) | ('+' in str(data.loc[i,'nat']))):
        data.loc[i,'nat']=1
    elif ((data.loc[i,'nat']==0) | ('-' in str(data.loc[i,'nat']))):
        data.loc[i,'nat']=0
    # next evaluate based on the pTNM prefix 
    if str(data.loc[i,'tumor_staging_category']).startswith('y'):
        data.loc[i,'nat']=1
    elif str(data.loc[i,'tumor_staging_category']).startswith('p'):
        data.loc[i,'nat']=0
    # 5d - clinical tumor size 
    if type(data.loc[i,'clin_tumor_size'])==str:
        if ((data.loc[i,'clin_tumor_size']=='.') | (data.loc[i,'clin_tumor_size']=='..')):
            data.loc[i,'clin_tumor_size']=None
        else:
            nums=re.findall(r'[0-9\.]+', data.loc[i,'clin_tumor_size'])
            if len(nums)>0:
                if nums[0]!='.':
                    data.loc[i,'clin_tumor_size']=np.max([float(x) for x in nums])
                else:
                    data.loc[i,'clin_tumor_size']=None
            else:
                data.loc[i,'clin_tumor_size']=None
    # Note - her2_therapy, her2_therapy_type, cnb_grade are all empty 

data.surgery_type=data.surgery_type.replace({'Tm' : 'Mastectomy (Tm)', 'Bcs' : 'Lumpectomy (Bcs)', 'L;Tm, R;Tm' : 'Mastectomy (Tm)'})
data.response_nat=data.response_nat.replace(
    {'.':pd.NA, 
    '..':pd.NA, 
    'CR':'No residual disease', 
    'CR(PR?)':'No residual disease',
    'CR(PR)':'No residual disease', 
    'SD':'Residual disease', 
    'PR':'Residual disease', 
    'PD':'Residual disease', 
    'TCx6+TAM':pd.NA,
    'SD( I&D/excisional biopsy was done before NAC, margin positive)':'Residual disease'}
)
data.nat_reg=data.nat_reg.replace(
    {'.':pd.NA, 
    'T':'Taxane/DOC/PAC/TC', 
    'Both A/T':'Both', 
    'Other(Carboplatin+Eribulin)':'Other', 
    'T(TCHP)':'Taxane/DOC/PAC/TC', 
    'A(+Carbo/Eribulin); On NU10B7 trial with Carbo/Eribulin.  Completed therapy and had pCR. s/p AC4.':'Anthracycline/FEC/EC/AC', 
    'Both A/T(522 Pembro/Taxol/Carbo/AC)':'Both'}
)

## Identify columns that are in REDCap features but not in our current table & add them 
for item in redcap_features:
    if item not in data.columns:
        data[item]=None

## Reorder features to match REDCap template 
data=data[redcap_features]

## Reindex + fill in record number 
data.record_id=data.index
data.number=data.index

#######################################
#### Fill in the "complete" checks ####
#######################################

category_colnames={
    'demographics' : ['number','epic_mrn','dob','first_name','last_name','age_at_diagnosis','year_of_diagnosis'],
    'parity_information' : ['parous','number_pregnancies','number_births','age_at_first_pregnancy','age_at_most_recent_pregnancy','years_since_pregnancy'],
    'tumor_characteristics' : ['histology','tumor_size','histologic_grade','tumor_staging_category','node_staging_category','er_status','pr_status','her2_status','ki67'],
    'geneticsfam_hx' : ['any_patho_mutation', 'any_vus_mutation', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm', 'fam_hx'],
    'treatment' : ['nat', 'nat_reg', 'her2_therapy', 'her2_therapy_type', 'response_nat']
}

for i in data.index:
    for key in category_colnames.keys():
        if np.all(~pd.isnull(data.loc[i,category_colnames[key]])):
            data.loc[i,f'{key}_complete']=2 # 2 indicates "Complete" 
        elif key=='parity_information':
            if pd.isnull(data.loc[i,'parous']): # parous is null -> incomplete 
                data.loc[i,f'{key}_complete']=0
            elif data.loc[i,'number_births']==0: # number_births is zero -> ages unnecessary (nulliparous)
                if np.all(~pd.isnull(data.loc[i,['parous','number_pregnancies','number_births']])):
                    data.loc[i,f'{key}_complete']=2
                else:
                    data.loc[i,f'{key}_complete']=0
            else:
                data.loc[i,f'{key}_complete']=0 # parous is nonzero -> incomplete (all other info needs to be entered)
        elif key=='treatment':
            if pd.isnull(data.loc[i,'nat']):
                data.loc[i,f'{key}_complete']=0
            elif data.loc[i,'nat']==0:
                data.loc[i,f'{key}_complete']=2
            else:
                data.loc[i,f'{key}_complete']=0
        else:
            data.loc[i,f'{key}_complete']=0 # 0 indicates "Incomplete" (I think 1 is "Unverified")


############################################
#### Re-map data to REDCap coding (omg) ####
############################################

data.breastfed=data.breastfed.replace({'Yes': 1, 'No': 2, 'Unknown': 3}).astype('Int64')

data.presentation=data.presentation.replace({'Self-detected': 1, 'Physician-detected': 2, 'Screening imaging-detected': 3})
data.surgery_type=data.surgery_type.replace({'Mastectomy (Tm)':1, 'Lumpectomy (Bcs)':2, 'Biopsy-only':3, 'Other':4}).astype('Int64')
data.histology=data.histology.replace({'IDC':1, 'ILC':2, 'DCIS':3, 'Other':4}).astype('Int64')
data.histologic_grade=data.histologic_grade.replace({'No histologic grade (DCIS)':1, '1':2, '2':3, '3':4, '.':None}).astype('Int64')
data.tumor_staging_category=data.tumor_staging_category.replace({'pTx': 1, 'pTis': 2, 'pT1': 3, 'pT1mi': 4, 'pT2': 5, 'pT3': 6, 'pT4': 7, 'ypTx': 8, 'ypT0': 9, 'ypTis': 10, 'ypT1': 11, 'ypT1mi': 12, 'ypT2': 13, 'ypT3': 14, 'ypT4': 15}).astype('Int64')
data.node_staging_category=data.node_staging_category.replace({'pNX': 1, 'pN0': 2, 'pN1': 3, 'pN1mi': 4, 'pN2': 5, 'pN3': 6, 'ypNX': 7, 'ypN0': 8, 'ypN1': 9, 'ypN1mi': 10, 'ypN2': 11, 'ypN3': 12}).astype('Int64')
data.er_status=data.er_status.replace({'Positive':1, 'Negative':2}).astype('Int64')
data.pr_status=data.pr_status.replace({'Positive':1, 'Negative':2}).astype('Int64')
data.her2_status=data.her2_status.replace({'Positive':1, 'Negative':2, 'Not performed':3}).astype('Int64')
data.ki67=data.ki67.replace({'Low':1, 'Intermediate':2, 'High':3, 'Not performed':4}).astype('Int64')

for colname in ['any_patho_mutation', 'any_vus_mutation']:
    data[colname]=data[colname].replace({'Present':1, 'Absent':2})

for colname in ['brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']:
    data[colname]=data[colname].replace({'Pathogenic':1, 'VUS':2, 'Negative':3, 'Not performed':4})

data.nat_reg=data.nat_reg.replace({'Anthracycline/FEC/EC/AC':1, 'Taxane/DOC/PAC/TC':2, 'Both':3, 'Other':4})
data.her2_therapy=data.her2_therapy.replace({'Present':1, 'Absent':2})
data.her2_therapy_type=data.her2_therapy_type.replace({'Trastuzumab':1, 'Trastuzumab + more':2})
data.measure_clin_size=data.measure_clin_size.replace({'Ultrasound (US)':1, 'Magnetic resonance imaging (MRI)':2, 'Mammography (MM)':3})
data.clin_tumor_stag_cat=data.clin_tumor_stag_cat.replace({'TX': 1, 'Tis': 2, 'T1': 3, 'T1mi': 4, 'T2': 5, 'T3': 6, 'T4': 7}).astype('Int64')
data.clin_node_stag_cat=data.clin_node_stag_cat.replace({'NX':1, 'N0':2, 'N1':3, 'N2':4, 'N3':5}).astype('Int64')
data.response_nat=data.response_nat.replace({'No residual disease':1, 'Residual disease':2, 'Unknown':3}).astype('Int64')

###################
#### Save data ####
###################

data.to_csv(f'{proj}/data/20220616_redcap_import_data.csv', index=False)
