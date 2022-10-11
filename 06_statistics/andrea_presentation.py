import os
import sys
import re
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import chisquare, ttest_ind
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

outdir = 'stat_results/nac_only'

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

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    oddsratio = {}
    for colnum in range(X.shape[1]):
        oddsratio[colnum] = []
    # 
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
        if ~np.all(y_bs==0):
            lrm = LogisticRegression(penalty='l2', solver='lbfgs')
            lrm.fit(X_bs, y_bs)
            for colnum in range(X.shape[1]):
                oddsratio[colnum].append(np.exp(lrm.coef_[0][colnum]))
        else:
            continue
    mean_or = [np.mean(oddsratio[colnum]) for colnum in range(X.shape[1])]
    ci = []
    for colnum in range(X.shape[1]):
        # first get ci1
        p = ((1.0-alpha)/2.0) * 100
        lower = max(0.0, np.percentile(oddsratio[colnum], p))
        p = (alpha+((1.0-alpha)/2.0)) * 100
        upper = np.percentile(oddsratio[colnum], p)
        ci.append((lower, upper))
    return mean_or, ci

##########################################################
#### Stratified analysis 1 - three biomarker subtypes ####
##########################################################

#### A. Chi-squared/Fisher's exact test ####
for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    crosstab = pd.crosstab(
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    stattab = pd.concat((
        pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
        pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    ))
    stattab.columns.name = None # remove column label for clarity 
    fn_bm_cat = re.sub(r"\s|[^A-Za-z]+", "_", bm_cat).lower()
    pd.DataFrame(
        np.array([chisquare(stattab).statistic, chisquare(stattab).pvalue]), 
        index=['chi-stat', 'p'], 
        columns=stattab.columns
    ).to_csv(f'{dn}/{outdir}/chi_test_biomarksub_{fn_bm_cat}.csv')

#### B. Odds-ratio and confidence intervals calculated via logistic regression ####
record_ors = pd.DataFrame(
    None, 
    index=['<5 years', '5-10 years', '>=10 years'],
    columns=['ER/PR+ HER2-', 'Triple Negative', 'HER2+']
)

for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
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
        record_ors.loc[par_query, bm_cat] = f'{or_mean[0]:.4f} (95% CI {ci[0][0]:.4f}-{ci[0][1]:.4f})'

record_ors.to_csv(f'{dn}/{outdir}/logistic_oddsratio_biomarkercategory.csv', index=True)

#### C. Student's t-test with continuous RCB ####
record_t = pd.DataFrame(
    None, 
    index=['<5 years', '5-10 years', '>=10 years'],
    columns=['ER/PR+ HER2-', 'Triple Negative', 'HER2+']
)

for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    sub_data = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
    for par_query in ['<5 years', '5-10 years', '>=10 years']: # parity query category 
        ref_dist = sub_data.iloc[sub_data.parity_category.values=='Nulliparous',:].rcb.values
        alt_dist = sub_data.iloc[sub_data.parity_category.values==par_query,:].rcb.values
        t_result = ttest_ind(alt_dist, ref_dist, nan_policy='omit')
        t = t_result.statistic
        p = t_result.pvalue
        record_t.loc[par_query, bm_cat] = f'{t:.4f} (p={p:.4f})'

record_t.to_csv(f'{dn}/{outdir}/ttest_biomarkercategory.csv', index=True)

#########################################################
#### Stratified analysis 2 - ER/PR positive vs. rest ####
#########################################################

#### A. Chi-squared/Fisher's exact test ####
for erpr_cat in ['ER/PR+', 'ER/PR-']:
    crosstab = pd.crosstab(
        data.iloc[data.er_pr.values==erpr_cat,:].rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.er_pr.values==erpr_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    stattab = pd.concat((
        pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
        pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    ))
    stattab.columns.name = None # remove column label for clarity 
    if '+' in erpr_cat:
        fn_bm_cat = 'pos'
    elif '-' in erpr_cat:
        fn_bm_cat = 'neg'
    pd.DataFrame(
        np.array([chisquare(stattab).statistic, chisquare(stattab).pvalue]), 
        index=['chi-stat', 'p'], 
        columns=stattab.columns
    ).to_csv(f'{dn}/{outdir}/chi_test_erpr_{fn_bm_cat}.csv')

#### B. Odds-ratio and confidence intervals calculated via logistic regression ####
record_ors = pd.DataFrame(
    None, 
    index=['<5 years', '5-10 years', '>=10 years'],
    columns=['ER/PR+', 'ER/PR-']
)

for erpr_cat in ['ER/PR+', 'ER/PR-']:
    sub_data = data.iloc[data.er_pr.values==erpr_cat,:]
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
        record_ors.loc[par_query, erpr_cat] = f'{or_mean[0]:.4f} (95% CI {ci[0][0]:.4f}-{ci[0][1]:.4f})'

record_ors.to_csv(f'{dn}/{outdir}/logistic_oddsratio_erpr_category.csv', index=True)

#### C. Student's t-test with continuous RCB ####
record_t = pd.DataFrame(
    None, 
    index=['<5 years', '5-10 years', '>=10 years'],
    columns=['ER/PR+', 'ER/PR-']
)

for erpr_cat in ['ER/PR+', 'ER/PR-']:
    sub_data = data.iloc[data.er_pr.values==erpr_cat,:]
    for par_query in ['<5 years', '5-10 years', '>=10 years']: # parity query category 
        ref_dist = sub_data.iloc[sub_data.parity_category.values=='Nulliparous',:].rcb.values
        alt_dist = sub_data.iloc[sub_data.parity_category.values==par_query,:].rcb.values
        t_result = ttest_ind(alt_dist, ref_dist, nan_policy='omit')
        t = t_result.statistic
        p = t_result.pvalue
        record_t.loc[par_query, erpr_cat] = f'{t:.4f} (p={p:.4f})'

record_t.to_csv(f'{dn}/{outdir}/ttest_erpr_category.csv', index=True)

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
