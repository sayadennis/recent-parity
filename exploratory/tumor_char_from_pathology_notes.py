import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import staging_summary_mining_functions

mike_dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/tumor_char'
# saya_dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/pathology'
dout = '/share/fsmresfiles/breast_cancer_pregnancy/data/exploratory'

# m_tr = pd.read_csv(f'{mike_dn}/Tumor_Registry_Tumor_Characterstics.tsv', sep='\t') # mike - tumor registry 
m_pt = pd.read_csv(f'{mike_dn}/Pathology_Cases_with_Surgeries_and_Abstractions.tsv', sep='\t') # mike - pathology cases abstractions 
# s_pt = pd.read_csv(f'{saya_dn}/isolated_sections_pathology.csv') # saya - pathology cases abstractions 

m_pt_sizes = pd.read_csv(f'{mike_dn}/Pathology_Cases_with_Surgeries_and_Abstractions_tumor_size_suggestions.tsv', sep='\t')
m_pt_sizes.rename({'procedure_occurrence_stable_identifier_value' : 'accession number'}, axis=1, inplace=True)

# Join the original pathology table and the size suggestion table 
m_pt = pd.merge(m_pt, m_pt_sizes, how='left', on='accession number')

############################################
#### Pathology notes abstraction method ####
############################################

## Fist record how many records are present per patient 
n_patient_pre = len(m_pt['nmhc mrn'].unique())
n_recs_pre = []
for nmhc in m_pt['nmhc mrn'].unique():
    obs=m_pt.iloc[m_pt['nmhc mrn'].values==nmhc,:]
    n_recs_pre.append(obs.shape[0])

## remove outside cases for simplicity for now 
m_pt = m_pt.iloc[m_pt['pathology case type'].values!='Outside Surgical Pathology',:]

## select path notes from 2010-2020
m_pt['accessioned year'] = [int(x.split('/')[-1]) for x in m_pt['accessioned date']]
m_pt = m_pt.iloc[[((x>=2010) & (x<=2020)) for x in m_pt['accessioned year']],:]

data = m_pt.iloc[m_pt['cancer site filter'].values=='yes',:][[
    'nmhc mrn',
    'accessioned date',
    'accessioned year',
    'pathology note title',
    'surgery name',
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

# How many patients are there? 
total = len(data['nmhc mrn'].unique())

############################################
#### Identifying staging summary tables ####
############################################

## select rows where path note includes a mastectomy or lumpectomy
# m_pt = m_pt.iloc[[(('astectomy' in x) | ('umpectomy' in x)) for x in m_pt['pathology note text']],:]

table_patterns = [
    ## generic "breast cancer" 
    '\nBREAST CANCER STAGING SUMMARY\n',
    '\nLEFT BREAST CANCER STAGING SUMMARY\n',
    '\nRIGHT BREAST CANCER STAGING SUMMARY\n',
    '\nBreast Carcinoma Checklist:',
    ## DCIS 
    '\nDUCTAL CARCINOMA IN SITU SUMMARY', # '\n' , in case it looks like 'DUCTAL CARCINOMA IN SITU SUMMARY (Combined with 01-S-11-33962)'
    '\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n',
    '\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary',
    '\n Ductal Carcinoma In Situ Staging Summary', # ' \n' ':' ','
    '\nIn Situ Breast Carcinoma Checklist:\n', # not necessarily DCIS but could be LCIS? 
    '\nIn-Situ Breast Carcinoma Checklist:\n',
    ## Invasive 
    '\n Invasive Breast Cancer Staging Summary', # ' ', '-', ':', ','
    '\nInvasive Breast Carcinoma Checklist', # '\n', ':'
    '\nINVASIVE BREAST CANCER SUMMARY', # '\n' , in case it looks like 'INVASIVE BREAST CANCER SUMMARY (LEFT BREAST)' 
    '\nINVASIVE RIGHT BREAST CANCER SUMMARY\n',
    '\nINVASIVE LEFT BREAST CANCER SUMMARY\n',
    '\n Microinvasive Breast Cancer Staging Summary',
    ## postneoadjuvant 
    '\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary', # '  \n'
    '\n POST NEOADJUUVANT Invasive Breast Cancer Staging Summary',
    '\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary',
    '\n Post-treatment Breast Cancer Staging Summary',
    '\n Post-Treatment Invasive Breast Cancer Staging Summary',
    '\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ', # ' \n'
    '\n Invasive Breast Cancer POST-TREATMENT Staging Summary',
    '\n Invasive Breast Cancer Post-Treatment Staging Summary',
    '\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary',
    # '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary - RIGHT BREAST',
    '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ',
    '\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n',
]

## select only the rows where there is a path table based on above patterns 
hastable = [np.any([len(re.findall(x, note_text))>0 for x in table_patterns]) for note_text in data['pathology note text']]

notable = data.iloc[~np.array(hastable),:]
data = data.iloc[hastable,:]

## Next record how many records are present per patient after selecting notes with table
n_patient_post = len(data['nmhc mrn'].unique())
n_recs_post = []
for nmhc in data['nmhc mrn'].unique():
    obs=data.iloc[data['nmhc mrn'].values==nmhc,:]
    n_recs_post.append(obs.shape[0])

plt.violinplot([n_recs_pre, n_recs_post], np.arange(2), showextrema=True, showmedians=True, bw_method='silverman')
plt.ylabel('Number of records')
plt.xticks([0,1], ['Pre', 'Post'])
plt.title('# path records per patient pre/post-selection via staging summary tables')
plt.tight_layout()
plt.savefig('recent_parity_n_path_records_violin.png')
plt.close()

# ## Show how many notes are present with each of the patterns (quick and dirty) 
# [np.sum([table_patterns[i] in x for x in data['pathology note text']]) for i in range(len(table_patterns))]

##################################################
#### Record tumor characteristics from tables ####
##################################################

###################
#### Histology ####
###################

colname='cancer histology'
data[colname]=None
print(f'Original number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

# table_pattern='\n Post-treatment Breast Cancer Staging Summary'
# test = data.iloc[[table_pattern in x for x in data['pathology note text']],:]

for i in data.index:
    for pattern in table_patterns:
        if pattern in data.loc[i,'pathology note text']:
            tum_grade=get_histology(note_text=data.loc[i,'pathology note text'], table_title=pattern)
            if tum_grade is not None:
                data.loc[i,colname] = tum_grade

print(f'New number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

###############
#### Grade ####
###############

colname='cancer pathologic sbr grade'
data[colname]=None
print(f'Original number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

for i in data.index:
    for pattern in table_patterns:
        if pattern in data.loc[i,'pathology note text']:
            tum_grade=get_grade(note_text=data.loc[i,'pathology note text'], table_title=pattern)
            if tum_grade is not None:
                data.loc[i,colname] = tum_grade

print(f'New number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

##############
#### Size ####
##############

colname='tumor size'
data[colname]=None
print(f'Original number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

for i in data.index:
    for pattern in table_patterns:
        if pattern in data.loc[i,'pathology note text']:
            tum_size=get_size(note_text=data.loc[i,'pathology note text'], table_title=pattern)
            if tum_size is not None:
                data.loc[i,'tumor size'] = tum_size

print(f'New number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

###############
#### Ki-67 ####
###############

colname='ki67'
data[colname]=None
print(f'Original number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

for i in data.index:
    for pattern in table_patterns:
        if pattern in data.loc[i,'pathology note text']:
            tum_grade=get_ki67(note_text=data.loc[i,'pathology note text'], table_title=pattern)
            if tum_grade is not None:
                data.loc[i,colname] = tum_grade

print(f'New number of patients with {colname} available: ', len(data['nmhc mrn'].iloc[~(pd.isnull(data[colname]).values | (data[colname]=='not applicable').values)].unique()))

############################################################
#### Clean up size, grade, histology, Ki-67 ####
############################################################

## Size (some of Mike's entries were numeric but strings, like '1.7')
for i in data.index:
    if type(data.loc[i,'tumor size'])==str:
        try:
            data.loc[i,'tumor size']=float(data.loc[i,'tumor size'])
        except: # most likely "not applicable"
            data.loc[i,'tumor size']=None

# ## Node staging - replace NX with None
# for i in data.index:
#     if data.loc[i,'pathological nodes staging category']=='pNX':
#         data.loc[i,'pathological nodes staging category']=None

## change grades
for i in data.index:
    if data.loc[i,'cancer pathologic sbr grade']=='Grade 1':
        data.loc[i,'cancer pathologic sbr grade']=1
    elif data.loc[i,'cancer pathologic sbr grade']=='Grade 2':
        data.loc[i,'cancer pathologic sbr grade']=2
    elif data.loc[i,'cancer pathologic sbr grade']=='Grade 3':
        data.loc[i,'cancer pathologic sbr grade']=3

#############################################################
#### Consolidate duplicate rows in variables of interest ####
#############################################################

for i in data.index:
    for j in data.columns:
        if data.loc[i,j]=='not applicable':
            data.loc[i,j]=None

data.shape
slim = data[
    ['nmhc mrn', 
    'cancer histology', 
    'cancer pathologic sbr grade', 
    'tumor size', 
    'estrogen receptor status', 
    'progesterone receptor status', 
    'her2 status', 
    'lymphovascular invasion', 
    'ki67', 
    'number lymph nodes examined', 
    'number lymph nodes positive tumor',
    'pathological tumor staging category', 
    'pathological nodes staging category']
].drop_duplicates(ignore_index=True)

def consolidate(df_slice): # slice of data frame. there are either 2, 3, or 4 rows (6/10/2021)
    consolidated = {}
    # get column names
    colnames = df_slice.columns
    # get NMHC MRN 
    nmhc = df_slice['nmhc mrn'].values[0]
    # iterate through columns
    isnan = pd.isnull(df_slice)
    for colname in colnames:
        nonnan_slice = df_slice.iloc[~isnan[colname].values,:]
        if np.all(isnan.loc[:,colname]): # all NaN
            continue
        elif np.sum(~isnan.loc[:,colname])==1: # only one is not NaN
            consolidated[colname] = list(df_slice.loc[:,colname].values)[np.where(~isnan.loc[:,colname])[0][0]]
        else: # multiple non-NaN values
            if np.all([nonnan_slice.loc[i,colname]==nonnan_slice.loc[nonnan_slice.index[0],colname] for i in nonnan_slice.index]): # if all non-NaN values match with each other 
                consolidated[colname] = df_slice.loc[df_slice.index[0],colname]
            elif colname=='tumor size': # take largest value
                consolidated[colname] = np.max(df_slice[colname].values)
            else:
                if colname=='cancer histology':
                    if 'IDC' in df_slice[colname].values:
                        consolidated[colname]='IDC'
                    elif 'ILC' in df_slice[colname].values:
                        consolidated[colname]='ILC'
                    elif 'DCIS' in df_slice[colname].values:
                        consolidated[colname]='DCIS'
                elif colname in ['estrogen receptor status', 'progesterone receptor status', 'her2 status']:
                    if 'positive' in df_slice[colname].values:
                        consolidated[colname]='positive'
                    else:
                        consolidated[colname]='negative'
                elif colname=='cancer pathologic sbr grade':
                    consolidated[colname]=max(df_slice[colname].values)
                elif colname=='ki67':
                    only_nums=[num for num in list(df_slice[colname].values) if isinstance(num, (int,float))]
                    if len(only_nums)>0:
                        consolidated[colname]=max(only_nums)
                    else:
                        consolidated[colname]=list(df_slice[colname].values)[0]
                elif colname=='pathological tumor staging category':
                    for t_cat in ['pT4', 'pT3', 'pT2', 'pT1mi', 'pT1', 'pT0', 'pTX']: # descending order of severity 
                        if np.any([t_cat in x for x in nonnan_slice[colname].values]):
                            consolidated[colname]=t_cat
                            break # stop after first instance of finding 
                elif colname=='pathological nodes staging category':
                    for n_cat in ['pN3', 'pN2', 'pN1', 'pN0', 'pNX']: # descending order of severity 
                        if np.any([n_cat in x for x in nonnan_slice[colname].values]):
                            consolidated[colname]=n_cat
                            break # stop after first instance of finding 
                elif colname in ['number lymph nodes examined', 'number lymph nodes positive tumor']:
                    consolidated[colname]=max([int(num) for num in nonnan_slice[colname].values if num.isnumeric()])
                else:
                    print(f'nonmatching values in MRN {nmhc} for {colname}: \n', df_slice[colname].values) # 
                    # continue
    return consolidated

singles = slim.iloc[[sum(slim['nmhc mrn']==x)==1 for x in slim['nmhc mrn']],:]

dup = slim.iloc[[sum(slim['nmhc mrn']==x)>1 for x in slim['nmhc mrn']],:]
dup_ids = dup['nmhc mrn'].unique()

dup_converted_to_singles = pd.DataFrame()

for dup_id in dup_ids:
    # test = consolidate(dup.iloc[dup['nmhc mrn'].values==dup_id,:])
    dup_converted_to_singles = dup_converted_to_singles.append(consolidate(dup.iloc[dup['nmhc mrn'].values==dup_id,:]), ignore_index=True)

singles_add=pd.concat((singles, dup_converted_to_singles))

# ###########################
# #### Plot distribution ####
# ###########################

# plot_items = [
#     'cancer histology',
#     'tumor size',
#     'cancer pathologic sbr grade', 
#     'pathological tumor staging category',
#     'pathological nodes staging category',
#     # 'pathological metastasis staging category', 
#     'estrogen receptor status',
#     'progesterone receptor status',
#     'her2 status'
# ]

# data = data[['nmhc mrn'] + plot_items]

# fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16,12))
# for i in range(len(plot_items)):
#     item = plot_items[i]
#     miss_rate = 100*pd.isnull(data[item]).sum()/data.shape[0]
#     if item=='tumor size':
#         axs[i//3, i%3].hist(data[item].iloc[~pd.isnull(data[item]).values], bins=30)
#         axs[i//3, i%3].set_title(f'{item} ({miss_rate:.1f}% missing)')
#         axs[i//3, i%3].set_xlabel('size (cm)')
#     else:
#         cts = data[item].value_counts()
#         if 'not applicable' in cts.index:
#             cts.drop('not applicable', axis=0, inplace=True)
#         if cts.shape[0]>5:
#             n_others = cts.iloc[5:].sum()
#             cts = cts.iloc[:5]
#             cts = pd.concat((cts, pd.Series(n_others, index=['others'])), axis=0)
#         axs[i//3, i%3].bar(np.arange(len(cts)), cts)
#         axs[i//3, i%3].set_xticks(np.arange(len(cts)))
#         axs[i//3, i%3].set_xticklabels(cts.index, rotation=30, ha='right') # , fontsize=14
#         item = re.sub('cancer', '', item)
#         item = re.sub('pathological', 'path', item)
#         axs[i//3, i%3].set_title(f'{item} ({miss_rate:.1f}% missing)')
#     # 
#     if i%3==0:
#         axs[i//3, i%3].set_ylabel('Counts') # , fontsize=14

# fig.suptitle('Current distributions of tumor characteristics')
# plt.tight_layout()
# fig.savefig('recent_parity_tum_char_distributions.png')
# plt.close()

###################################################################################################
#### Save spreadsheets separately for patients who have all info available and have duplicates ####
###################################################################################################

## Clean the no-table subset 
notable_slim = notable[
    ['nmhc mrn', 
    'cancer histology', 
    'cancer pathologic sbr grade', 
    'tumor size', 
    'estrogen receptor status', 
    'progesterone receptor status', 
    'her2 status', 
    'lymphovascular invasion', 
    'ki67', 
    'number lymph nodes examined', 
    'number lymph nodes positive tumor',
    'pathological tumor staging category', 
    'pathological nodes staging category']
].drop_duplicates(ignore_index=True)

for i in notable_slim.index:
    for j in notable_slim.columns:
        if notable_slim.loc[i,j]=='not applicable':
            notable_slim.loc[i,j]=None

for i in notable_slim.index:
    # first re-write histology 
    if notable_slim.loc[i,'cancer histology']=='intraductal carcinoma, noninfiltrating (8500/2)':
        notable_slim.loc[i,'cancer histology']='DCIS'
    elif notable_slim.loc[i,'cancer histology']=='infiltrating duct carcinoma (8500/3)':
        notable_slim.loc[i,'cancer histology']='IDC'
    else:
        notable_slim.loc[i,'cancer histology']=None
    # next rewrite grade 
    if notable_slim.loc[i,'cancer pathologic sbr grade']=='Grade 3':
        notable_slim.loc[i,'cancer pathologic sbr grade']=3
    elif notable_slim.loc[i,'cancer pathologic sbr grade']=='Grade 2':
        notable_slim.loc[i,'cancer pathologic sbr grade']=2
    elif notable_slim.loc[i,'cancer pathologic sbr grade']=='Grade 1':
        notable_slim.loc[i,'cancer pathologic sbr grade']=1
    else:
        notable_slim.loc[i,'cancer pathologic sbr grade']=None


## remove patients that already have clear records in the other table 
notable_slim = notable_slim.iloc[[(x not in data['nmhc mrn'].unique()) for x in notable_slim['nmhc mrn']],:]

## Load the earlier table 
combined = pd.read_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/recent_parity_combined_data_reviewed.csv')
combined = combined.iloc[combined['Review complete'].values=='x',]

## Join the previous results to the tumor characteristics 
singles_joined = pd.merge(combined, singles_add, how='left', left_on='EPIC_mrn', right_on='nmhc mrn')
dup_joined = pd.merge(combined, dup_converted_to_singles, how='left', left_on='EPIC_mrn', right_on='nmhc mrn')
notable_slim_joined = pd.merge(combined, notable_slim, how='left', left_on='EPIC_mrn', right_on='nmhc mrn')

## calculate missing rate for tumor characteristics 
singles_missing = pd.isnull(singles_joined.iloc[~pd.isnull(singles_joined['nmhc mrn']).values,:][plot_items]).sum().sum()/(len(plot_items)*singles_joined.iloc[~pd.isnull(singles_joined['nmhc mrn']).values,:].shape[0])
dup_missing = pd.isnull(dup_joined.iloc[~pd.isnull(dup_joined['nmhc mrn']).values,:][plot_items]).sum().sum()/(len(plot_items)*dup_joined.iloc[~pd.isnull(dup_joined['nmhc mrn']).values,:].shape[0])
notable_missing = pd.isnull(notable_slim_joined.iloc[~pd.isnull(notable_slim_joined['nmhc mrn']).values,:][plot_items]).sum().sum()/(len(plot_items)*notable_slim_joined.iloc[~pd.isnull(notable_slim_joined['nmhc mrn']).values,:].shape[0])

## Save tables 
singles_joined.to_csv(f'{dout}/singles.csv', index=False) # 905 rows, 905 patients 
dup_joined.to_csv(f'{dout}/multiple.csv', index=False) # 42 patients, 85 rows
notable_slim_joined.to_csv(f'{dout}/no_table.csv', index=False) # 114 patients, 162 rows 
