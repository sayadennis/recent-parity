import numpy as np
import pandas as pd

## Load data ## 
fin='/share/fsmresfiles/breast_cancer_pregnancy/data/20220530_redcap_import_data.csv'
data=pd.read_csv(fin)

## Demographic and gynecological history ##
varnames=[
    'age_at_diagnosis',
    'menarche',
    # 'parous',
    'number_pregnancies',
    'number_births',
    'age_at_first_pregnancy',
    'age_at_most_recent_pregnancy',
    # 'breastfed',
    'breastfeeding_duration',
    'recency_of_lactation'
]

varname_dict={
    'age_at_diagnosis' : 'Age at Diagnosis',
    'menarche' : 'Age at Menarche',
    # 'parous' : 'Parous', # categorical
    'number_pregnancies' : 'Number of Pregnancies',
    'number_births' : 'Number of Live Births',
    'age_at_first_pregnancy' : 'Age at First Pregnancy',
    'age_at_most_recent_pregnancy' : 'Age at Most Recent Pregnancy',
    # 'breastfed' : 'Breastfed', # categorical
    'breastfeeding_duration' : 'Duration of Breastfeeding (months)',
    'recency_of_lactation' : 'Recency of lactation (years)'
}

sumtab=pd.DataFrame(index=varname_dict.values(), columns=['mean','SD'])

for varname in varnames:
    sumtab.loc[varname_dict[varname],'mean']=np.nanmean(data[varname])
    sumtab.loc[varname_dict[varname],'SD']=np.nanstd(data[varname])

# proportion of parous 
print(data.parous.mean())

# proportion of records for breastfeeding
print(data.breastfed.sum()/data.shape[0])

sumtab.to_csv('demo_gyn_sumtab.csv')

## Genetic testing ## 
varnames=['brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']

for varname in ['any_patho_mutation', 'any_vus_mutation']:
    print(f'Proportion of presence in {varname}:', (data[varname]==1).sum()/data.shape[0])

num_encode={1:'Pathogenic', 2:'VUS', 3:'Negative', 4:'Not performed'}
sumtab=pd.DataFrame(index=[x.upper() for x in varnames], columns=['Pathogenic', 'VUS', 'Negative', 'Not performed'])
for varname in varnames:
    for i in np.arange(1,5):
        sumtab.loc[varname.upper(),num_encode[i]]=(data[varname].values==i).sum()

sumtab.to_csv('genetic_sumtab.csv')

## Tumor characteristics ## 
varnames=[
    'histology', 'tumor_size', 'histologic_grade', 'tumor_staging_category',
    'node_staging_category',
    'er_status', 'pr_status', 'her2_status', 'ki67']

# data.histology=data.histology.replace({'IDC':1, 'ILC':2, 'DCIS':3, 'Other':4}).astype('Int64')
# data.histologic_grade=data.histologic_grade.replace({'No histologic grade (DCIS)':1, '1':2, '2':3, '3':4, '.':None}).astype('Int64')
# data.tumor_staging_category=data.tumor_staging_category.replace({'pTx': 1, 'pTis': 2, 'pT1': 3, 'pT1mi': 4, 'pT2': 5, 'pT3': 6, 'pT4': 7, 'ypTx': 8, 'ypT0': 9, 'ypTis': 10, 'ypT1': 11, 'ypT1mi': 12, 'ypT2': 13, 'ypT3': 14, 'ypT4': 15}).astype('Int64')
# data.node_staging_category=data.node_staging_category.replace({'pNX': 1, 'pN0': 2, 'pN1': 3, 'pN1mi': 4, 'pN2': 5, 'pN3': 6, 'ypNX': 7, 'ypN0': 8, 'ypN1': 9, 'ypN1mi': 10, 'ypN2': 11, 'ypN3': 12}).astype('Int64')
# data.er_status=data.er_status.replace({'Positive':1, 'Negative':2}).astype('Int64')
# data.pr_status=data.pr_status.replace({'Positive':1, 'Negative':2}).astype('Int64')
# data.her2_status=data.her2_status.replace({'Positive':1, 'Negative':2, 'Not performed':3}).astype('Int64')
# data.ki67=data.ki67.replace({'Low':1, 'Intermediate':2, 'High':3, 'Not performed':4}).astype('Int64')

num_encode={
    'histology' : {'IDC':1, 'ILC':2, 'DCIS':3, 'Other':4},
    'histologic_grade' : {'No histologic grade (DCIS)':1, '1':2, '2':3, '3':4},
    'tumor_staging_category' : {'pTx': 1, 'pTis': 2, 'pT1': 3, 'pT1mi': 4, 'pT2': 5, 'pT3': 6, 'pT4': 7, 'ypTx': 8, 'ypT0': 9, 'ypTis': 10, 'ypT1': 11, 'ypT1mi': 12, 'ypT2': 13, 'ypT3': 14, 'ypT4': 15},
    'node_staging_category' : {'pNX': 1, 'pN0': 2, 'pN1': 3, 'pN1mi': 4, 'pN2': 5, 'pN3': 6, 'ypNX': 7, 'ypN0': 8, 'ypN1': 9, 'ypN1mi': 10, 'ypN2': 11, 'ypN3': 12},
    'er_status' : {'Positive':1, 'Negative':2},
    'pr_status' : {'Positive':1, 'Negative':2},
    'her2_status' : {'Positive':1, 'Negative':2, 'Not performed':3},
    'ki67' : {'Low':1, 'Intermediate':2, 'High':3, 'Not performed':4}
}

for varname in varnames:
    if varname in ['tumor_size', 'num_lymph_nodes_taken', 'num_ln_positive']:
        print(f'## {varname} ##')
        print('Mean', data[varname].mean(), '//', 'SD', data[varname].std())
        print()
    else:
        print(f'## {varname} ##')
        data[varname].replace(dict((v,k) for k,v in num_encode[varname].items())).value_counts().sort_values(ascending=False)/data[varname].replace(dict((v,k) for k,v in num_encode[varname].items())).value_counts().sum()
        print()

