import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mike_dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/tumor_char'
saya_dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/pathology'
dout = '/share/fsmresfiles/breast_cancer_pregnancy/data_summary/exploratory'

m_tr = pd.read_csv(f'{mike_dn}/Tumor_Registry_Tumor_Characterstics.tsv', sep='\t') # mike - tumor registry 
m_pt = pd.read_csv(f'{mike_dn}/Pathology_Cases_with_Surgeries_and_Abstractions.tsv', sep='\t') # mike - pathology cases abstractions 
s_pt = pd.read_csv(f'{saya_dn}/isolated_sections_pathology.csv') # saya - pathology cases abstractions 

m_pt_sizes = pd.read_csv(f'{mike_dn}/Pathology_Cases_with_Surgeries_and_Abstractions_tumor_size_suggestions.tsv', sep='\t')
m_pt_sizes.rename({'procedure_occurrence_stable_identifier_value' : 'accession number'}, axis=1, inplace=True)

# Join the original pathology table and the size suggestion table 
m_pt = pd.merge(m_pt, m_pt_sizes, how='left', on='accession number')

# m_tr['west mrn']; s_pt['EPIC_mrn']
# m_pt['nmhc mrn']; s_pt['Powerchart_mrn']

## Convert Epic MRN and Powerchart MRN to integers in my data table 
s_pt['EPIC_mrn'] = s_pt['EPIC_mrn'].astype(int)

s_pt['Powerchart_mrn'].fillna(0, inplace=True)
s_pt['Powerchart_mrn'] = s_pt['Powerchart_mrn'].astype(int)

## what is the number of patients that had records in my search? 
saya_records_present_epic_MRNs = s_pt['EPIC_mrn'].unique() # 384

## what is the number of patients that had records in Mike's tumor registry search? 
mike_tumreg_hashistology_epic_MRNs = m_tr['west mrn'].iloc[~pd.isnull(m_tr['histology']).values].unique() # 1160 
mike_tumreg_hassizeclin_epic_MRNs = m_tr['west mrn'].iloc[~pd.isnull(m_tr['Tumor Size Clin']).values].unique() # 460 
len(m_tr['west mrn'].iloc[~pd.isnull(m_tr['Grade Pathological']).values].unique())

###############################
#### Tumor registry method ####
###############################

data = m_tr[[
    'west mrn', # 1160
    'histology', 
    'Tumor Size Clin', 
    'Tumor Size Path', 
    'Tumor Size Summ',
    'CS Tum Size', 
    'Grade Clinical', 
    'Grade Pathological',
    'Grade',
    'Lymph Vascular Invasion', 
    'ERSummary', 
    'estrogen receptor assay',
    'PRSummary', 
    'progesterone receptor assay', 
    'HER2Overall Summ',
    'her2 summary result', 
    'Ki67', 
    'Nodes Positive', 
    'Nodes Examined',
    'Clin T TNM', 
    'Clin N TNM', 
    'Clin M TNM',
    'Path T TNM', 
    'Path N TNM', 
    'Path M TNM'
]].drop_duplicates()

total = len(data['west mrn'].unique())

n_records = {}
for colname in data.columns[1:]: # cols except for 'nmhc mrn' 
    n = len(data['west mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()) # number of records/patients with this data available
    n_records[colname] = n

plt.bar(n_records.keys(), n_records.values())
plt.xticks(rotation=40, ha='right')
plt.xlim(-1, len(n_records.keys()))
plt.hlines(y=total, linestyles='dashed', xmin=-1, xmax=len(n_records.keys())+1, colors='red')
plt.title('Data Availability from Tumor Registry')
plt.ylabel('Number of patients with records')
plt.tight_layout()
plt.savefig(f'{dout}/mike_tumor_registry_number_records_bar_graph.png')
plt.close()

############################################
#### Pathology notes abstraction method ####
############################################

data = m_pt.iloc[m_pt['cancer site filter'].values=='yes',:][[
    'nmhc mrn',
    'pathology note title',
    'pathology note text',
    'cancer site', 
    'cancer histology', 
    'tumor size', 
    'cancer pathologic sbr grade', 
    'pathological tumor staging category', 
    'pathological nodes staging category', 
    'pathological metastasis staging category', 
    'estrogen receptor status', 
    'progesterone receptor status', 
    'her2 status', 
    'lymphovascular invasion', 
    'ki67', 
    'p53', 
    'number lymph nodes examined', 
    'number lymph nodes positive tumor',
    'has_tumor_size',
    'has_tumor_size_suggestions'
]].drop_duplicates()

total = len(data['nmhc mrn'].unique())

## Plot the number of available records for each variable
n_records = {}
for colname in data.columns[1:]: # cols except for 'nmhc mrn' 
    n = len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()) # number of records/patients with this data available
    n_records[colname] = n

plt.bar(n_records.keys(), n_records.values())
plt.xticks(rotation=30, ha='right')
plt.xlim(-1, len(n_records.keys()))
plt.hlines(y=total, linestyles='dashed', xmin=-1, xmax=len(n_records.keys())+1, colors='red')
plt.title('Data Availability from Pathol Reports')
plt.ylabel('Number of patients with records')
plt.tight_layout()
plt.savefig(f'{dout}/mike_pathol_abstraction_number_records_bar_graph.png')
plt.close()

## Inspect the note titles and texts of those missing tumor sizes 
mrn_hassize = data.iloc[~((pd.isnull(data['tumor size']).values) | (data['tumor size'].values=='not applicable')),:]['nmhc mrn'].unique()
notes_hassize = data.iloc[~((pd.isnull(data['tumor size']).values) | (data['tumor size'].values=='not applicable')),:]
nosize = data.iloc[(x not in mrn_hassize for x in data['nmhc mrn'].values),:]

import re

#### Get pattern "Tumor size: xx cm" ####
mrn=[]
for i in range(nosize.shape[0]):
    if len(re.findall(r'Tumor [Ss]ize[\s\S]{0,200}:.{0,50} [ xto\-<\.0-9]+ cm', nosize['pathology note text'].iloc[i]))>0:
        mrn.append(nosize['nmhc mrn'].iloc[i])

len(mrn); nosize = nosize.iloc[[x not in mrn for x in nosize['nmhc mrn']],:]

#### Get pattern "INFILTRATING TUMOR TYPE: ... SIZE: ..." ####
mrn=[]
for i in range(nosize.shape[0]):
    if len(re.findall(r'INFILTRATING TUMOR TYPE:[\s\S]{0,500} SIZE: [ \.Xx0-9]+ CM', nosize['pathology note text'].iloc[i]))>0: #
        mrn.append(nosize['nmhc mrn'].iloc[i])

len(mrn); nosize = nosize.iloc[[x not in mrn for x in nosize['nmhc mrn']],:]

#### Get pattern "measuring xx cm in greatest dimension" ####
mrn=[]
for i in range(nosize.shape[0]):
    if len(re.findall(r'measuring [\.0-9]+ cm in greatest dimension', nosize['pathology note text'].iloc[i]))>0: #
        mrn.append(nosize['nmhc mrn'].iloc[i])

len(mrn); nosize = nosize.iloc[[x not in mrn for x in nosize['nmhc mrn']],:]

#### Get pattern "Size/Extent of Lesion: ..." ####
mrn=[]
for i in range(nosize.shape[0]):
    # if len(re.findall(r'\nSize/Extent of Lesion:  \n Present on one slide (largest diameter): .{0,10} \n Present on more than one slide: .{0,10} \n Number of slides with .{0,10}: .{0,10} \n Number of slides examined: .{0,10} \n Largest single diameter: .{0,10} \n Overall estimate of the size of .{0,10}: .{0,10}', nosize['pathology note text'].iloc[i]))>0:
    if len(re.findall(r'\nSize/Extent of Lesion:', nosize['pathology note text'].iloc[i]))>0:
        mrn.append(nosize['nmhc mrn'].iloc[i])

len(mrn); nosize = nosize.iloc[[x not in mrn for x in nosize['nmhc mrn']],:]

#### Get pattern "measures ... cm in greatest dimension" ####
# mrn=[]
# for i in range(nosize.shape[0]):
#     # if len(re.findall(r'\nSize/Extent of Lesion:  \n Present on one slide (largest diameter): .{0,10} \n Present on more than one slide: .{0,10} \n Number of slides with .{0,10}: .{0,10} \n Number of slides examined: .{0,10} \n Largest single diameter: .{0,10} \n Overall estimate of the size of .{0,10}: .{0,10}', nosize['pathology note text'].iloc[i]))>0:
#     if len(re.findall(r'[Mm]easures [\.0-9] cm', nosize['pathology note text'].iloc[i]))>0: #  in greatest dimension
#         mrn.append(nosize['nmhc mrn'].iloc[i])

# len(mrn); nosize = nosize.iloc[[x not in mrn for x in nosize['nmhc mrn']],:]


# ## inspection code ##
# for i in range(30):
#     print(nosize['pathology note text'].iloc[i], '\n\n\n\n\n')

############################################
#### Identifying staging summary tables ####
############################################

# the three common patterns found
# np.sum(['\n Invasive Breast Cancer Staging Summary  \n' in x for x in data['pathology note text']])
# np.sum(['\nBREAST CANCER STAGING SUMMARY\n' in x for x in data['pathology note text']])
# np.sum(['\nINVASIVE BREAST CANCER SUMMARY\n' in x for x in data['pathology note text']])

table_patterns = [
    ## generic "breast cancer" 
    '\nBREAST CANCER STAGING SUMMARY\n',
    '\nLEFT BREAST CANCER STAGING SUMMARY\n',
    '\nRIGHT BREAST CANCER STAGING SUMMARY\n',
    '\nBreast Carcinoma Checklist:',
    ## DCIS 
    '\nDUCTAL CARCINOMA IN SITU SUMMARY', # '\n' , in case it looks like 'DUCTAL CARCINOMA IN SITU SUMMARY (Combined with 01-S-11-33962)'
    '\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n',
    '\n Ductal Carcinoma In Situ Staging Summary ', # ' \n'
    # '\n Ductal Carcinoma In Situ Staging Summary for Left Breast  \n',
    # '\n Ductal Carcinoma In Situ Staging Summary for the LEFT BREAST (part B)  \n',
    # '\n Ductal Carcinoma In Situ Staging Summary for the RIGHT BREAST (part C)  \n',
    # '\n Ductal Carcinoma In Situ Staging Summary - LEFT BREAST  \n',
    # '\n Ductal Carcinoma In Situ Staging Summary (Right breast) with focal microinvasion  \n',
    '\nIn Situ Breast Carcinoma Checklist:\n', # not necessarily DCIS 
    ## Invasive 
    # '\n Invasive Breast Cancer Staging Summary  \n',
    '\n Invasive Breast Cancer Staging Summary ',
    '\n Invasive Breast Cancer Staging Summary-',
    '\n Invasive Breast Cancer Staging Summary:',
    '\n Invasive Breast Cancer Staging Summary - RIGHT BREAST  \n',
    '\n Invasive Breast Cancer Staging Summary--RIGHT BREAST  \n',
    '\n Invasive Breast Cancer Staging Summary--LEFT BREAST  \n',
    '\n Invasive Breast Cancer POST-TREATMENT Staging Summary  \n',
    '\nInvasive Breast Carcinoma Checklist:', # '\n' , in case it looks like 'Invasive Breast Carcinoma Checklist (left breast):'
    '\nINVASIVE BREAST CANCER SUMMARY', # '\n' , in case it looks like 'INVASIVE BREAST CANCER SUMMARY (LEFT BREAST)' 
    '\nINVASIVE RIGHT BREAST CANCER SUMMARY\n',
    '\nINVASIVE LEFT BREAST CANCER SUMMARY\n',
    '\n Invasive Breast Cancer Staging Summary (Left breast)  \n',
    '\n Invasive Breast Cancer Staging Summary for Right Breast  \n',
    '\n Microinvasive Breast Cancer Staging Summary',
    ## postneoadjuvant 
    '\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary  \n',
    '\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ', # ' \n'
    # '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary - RIGHT BREAST',
    '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ',
    '\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n',
]

hastable=[]
notable=[]
for mrn in data['nmhc mrn'].unique():
    notes=data['pathology note text'].iloc[data['nmhc mrn'].values==mrn]
    if np.any([np.any([pattern in note for note in notes]) for pattern in table_patterns]):
        hastable.append(mrn)
    else:
        notable.append(mrn)

len(hastable); len(notable)

notable = data.iloc[[x in notable for x in data['nmhc mrn']],:]

# print(notable['pathology note text'].iloc[11])

# Show tables that have the work [Mm]astectomy
mast = notable.iloc[['astectomy' in x for x in notable['pathology note text']],:]
# Show tables that have the work [Ll]umpectomy
lump = notable.iloc[['umpectomy' in x for x in notable['pathology note text']],:]

