import os
import sys
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import chisquare
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Read data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-10_0938.csv'

data = pd.read_csv(f'{dn}/{datadir}/{fn}')
data = data.iloc[data.nat.values==1,]

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

##########################
#### Define functions ####
##########################

# def generate_lrdata(df, feature=[], confounder=[], target=[], recency_thres=10):
#     lrdata = pd.DataFrame(None, index=None, columns=['parous', 'has_mutation', 'recent', 'age'])
#     for i in range(df.shape[0]):
#         pat_dict = {}
#         pat_dict['parous'] = [df['parous'].iloc[i]] # binary 0/1
#         pat_dict['recent'] = [int(df['years_since_pregnancy'].iloc[i] < recency_thres)]
#         pat_dict['ER/PR'] = [(df['er_status'].iloc[i]==1 | df['pr_status'].iloc[i]==1)]
#         pat_dict['target'] = [(df['rcb_category'].iloc[i]==1 | df['rcb_category'].iloc[i]==2)]
#         lrdata = pd.concat((lrdata, pd.DataFrame(pat_dict)), ignore_index=True)
#     # data for nulliparous vs. parous 
#     X_np = np.array(lrdata[['parous', 'age']], dtype=float)
#     X_np[:,1] = (X_np[:,1] - np.mean(X_np[:,1]))/np.std(X_np[:,1]) # standard scale
#     y_np = np.array(lrdata['has_mutation'], dtype=float)
#     # data for recency of parity --only select women who are parous
#     X_rec = np.array(lrdata.iloc[lrdata['parous'].values==1][['recent', 'age']], dtype=float)
#     X_rec[:,1] = (X_rec[:,1] - np.mean(X_rec[:,1]))/np.std(X_rec[:,1]) # standard scale
#     y_rec = np.array(lrdata.iloc[lrdata['parous'].values==1]['has_mutation'], dtype=float)
#     return X_np, y_np, X_rec, y_rec

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    oddsratio = {}
    for colnum in range(X.shape[1]):
        oddsratio[colnum] = []
    i = 0
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
        if ~np.all(y_bs==0):
            lrm = LogisticRegression(penalty='l2', solver='lbfgs')
            lrm.fit(X_bs, y_bs)
            for colnum in range(X.shape[1]):
                oddsratio[colnum].append(np.exp(lrm.coef_[0][colnum]))
        else:
            continue
    oddsratio = [np.mean(oddsratio[colnum]) for colnum in range(X.shape[1])]
    ci = []
    for colnum in range(X.shape[1]):
        # first get ci1
        p = ((1.0-alpha)/2.0) * 100
        lower = max(0.0, np.percentile(oddsratio[colnum], p))
        p = (alpha+((1.0-alpha)/2.0)) * 100
        upper = np.percentile(oddsratio[colnum], p)
        ci.append((lower, upper))
    return oddsratio, ci

##########################################################
#### Stratified analysis 1 - three biomarker subtypes ####
##########################################################

#### A. Chi-squared/Fisher's exact test ####
for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    print(f'#### {bm_cat} ####')
    crosstab = pd.crosstab(
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    stattab = pd.concat((
        pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
        pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    ))
    stattab.columns.name = None # remove column label for clarity 
    pd.DataFrame(
        np.array([chisquare(stattab).statistic, chisquare(stattab).pvalue]), 
        index=['chi-stat', 'p'], 
        columns=stattab.columns
    )
    print('\n\n')

#### B. Odds-ratio and confidence intervals calculated via logistic regression ####
for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    print(f'#### {bm_cat} ####')
    sub_data = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
    for par_query in ['<5 years', '5-10 years', '>=10 years']: # parity query category 
        X = sub_data.iloc[
            ((sub_data.parity_category.values=='Nulliparous') | (sub_data.parity_category.values==par_query))
        ].parity_category.map(
            {'Nulliparous' : 0, par_query : 1}
        ).to_numpy().reshape(-1,1)
        y = sub_data.iloc[
            ((sub_data.parity_category.values=='Nulliparous') | (sub_data.parity_category.values==par_query))
        ].rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']).map(
            {'0' : 1, '1' : 1, '2' : 0, '3' : 0}
        ).fillna(0).to_numpy()
        # perform LR analysis 
        or_mean, ci = get_oddsratio_ci(X, y)


#########################################################
#### Stratified analysis 2 - ER/PR positive vs. rest ####
#########################################################

#### A. Chi-squared/Fisher's exact test ####
for erpr_cat in ['ER/PR+', 'ER/PR-']:
    print(f'#### {erpr_cat} ####')
    crosstab = pd.crosstab(
        data.iloc[data.er_pr.values==erpr_cat,:].rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.er_pr.values==erpr_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    stattab = pd.concat((
        pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
        pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    ))
    stattab.columns.name = None # remove column label for clarity 
    pd.DataFrame(
        np.array([chisquare(stattab).statistic, chisquare(stattab).pvalue]), 
        index=['chi-stat', 'p'], 
        columns=stattab.columns
    )
    print('\n\n')

#### B. Odds-ratio and confidence intervals calculated via logistic regression ####

# #####################################
# #### Overall logistic regression ####
# #####################################

# """
# X_np, X_rec - input matrix for LR
# - "np" means parous/nulliparous comparison, and "rec" means recent/non-recent comparison 
# - input columns: 
#     1) gynecological history
#     2) confounder (1-2 columns, depending on how to categorize confounder)

# y_np, y_red - target array for LR. 0 means poor response and 1 means good response. 
# - poor vs. good response is defined by the RCB category. 0-1 is good and 2-3 is bad. 
# """

# #### A. with three biological subtypes ####

# #### B. with ER/PR positive vs. rest ####


