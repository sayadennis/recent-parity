import os
import sys
import re
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import chisquare, ttest_ind, f_oneway
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Read data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-10_0938.csv'

data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

## Fill in parity category 
data['parity_category'] = None
for i in data.index:
    if data.loc[i,'parous']==0:
        data.loc[i,'parity_category'] = 'Nulliparous'
    elif data.loc[i,'years_since_pregnancy']<5:
        data.loc[i,'parity_category'] = '<5 years'
    elif data.loc[i,'years_since_pregnancy']<10:
        data.loc[i,'parity_category'] = '5-10 years'
    else:
        data.loc[i,'parity_category'] = '>=10 years'

## Fill in biomarker subtypes 
data['biomarker_subtypes'] = None
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        er = False
        pr = False
        print('Missing ER/PR at index:', i)
    her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    if ((er | pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'ER/PR+ HER2-'
    elif ((not er) & (not pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'Triple Negative'
    elif her2:
        data.loc[i,'biomarker_subtypes'] = 'HER2+'
    else:
        print('Unknown pattern at index:', i)

data['er_pr'] = None
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        er = False
        pr = False
        print('Missing ER/PR at index:', i)
    if (er | pr):
        data.loc[i,'er_pr'] = 'ER/PR+'
    else:
        data.loc[i,'er_pr'] = 'ER/PR-'

summary = pd.read_csv(f'{dn}/summary_tables/super_summary_table.csv')
summary = summary.iloc[summary.super_category.values=='tumor_char',:].drop('super_category', axis=1)

#################################################################
#### Assess the effects of potentially confounding variables ####
#################################################################

# one-way ANOVA for continuous variables 
# chi-square/Fisher’s Exact tests for categorical variables 

## Age at diagnosis
ages = {par_cat : data.iloc[data.parity_category.values==par_cat,:].age_at_diagnosis.dropna().values 
        for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']}
stat, p = f_oneway(ages['Nulliparous'], ages['<5 years'], ages['5-10 years'], ages['>=10 years'])

## Age at first birth 
ages = {par_cat : data.iloc[data.parity_category.values==par_cat,:].age_at_first_pregnancy.dropna().values 
        for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']}
stat, p = f_oneway(ages['<5 years'], ages['5-10 years'], ages['>=10 years'])

## Age at last birth 
ages = {par_cat : data.iloc[data.parity_category.values==par_cat,:].age_at_most_recent_pregnancy.dropna().values 
        for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']}
stat, p = f_oneway(ages['<5 years'], ages['5-10 years'], ages['>=10 years'])

## Tumor size 
sizes = {par_cat : data.iloc[data.parity_category.values==par_cat,:].tumor_size.dropna().values 
        for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']}
stat, p = f_oneway(sizes['Nulliparous'], sizes['<5 years'], sizes['5-10 years'], sizes['>=10 years'])

## Tumor grade 
crosstab = pd.crosstab(data['histologic_grade'], data['parity_category']).rename(
    dd['histologic_grade']['Choices, Calculations, OR Slider Labels'],
    axis=0
)
stat, p = chisquare(crosstab, axis=0)

## Positive lymph nodes (high missing for our data) 
# data.num_ln_positive>0

##########################################################
#### Assess the associations between tumor and parity ####
##########################################################

## ER

## PR

## ER or PR 

## HER2

## Triple Negative 

## Biological subtypes 

