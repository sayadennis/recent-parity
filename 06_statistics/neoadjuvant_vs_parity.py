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

datestring = datetime.now().strftime("%Y%m%d")

###################
#### Read data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'

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
    # ER/PR status
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        er = False
        pr = False
        print('Missing ER/PR at index:', i)
    # HER2 status
    try:
        her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    except:
        her2 = False
        print('Missing HER2 at index:', i)
    # fill in subtypes
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
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].rcb_category.map(dd['response_nat']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    # stattab = pd.concat((
    #     pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
    #     pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    # ))
    stattab = crosstab.drop('Unknown', axis=0)
    stattab.columns.name = None # remove column label for clarity 
    fn_bm_cat = re.sub(r"\s|[^A-Za-z]+", "_", bm_cat).lower()
    pd.DataFrame(
        np.array([chisquare(stattab).statistic, chisquare(stattab).pvalue]), 
        index=['chi-stat', 'p'], 
        columns=stattab.columns
    ).to_csv(f'{dn}/{outdir}/{datestring}_chi_test_biomarksub_{fn_bm_cat}.csv')

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
        ].rcb_category.map(dd['response_nat']['Choices, Calculations, OR Slider Labels']).map(
            {'No residual disease' : 1, 'Residual disease' : 0, 'Unknown' : 0}
        ).fillna(0).to_numpy()
        # perform LR analysis 
        or_mean, ci = get_oddsratio_ci(X, y)
        record_ors.loc[par_query, bm_cat] = f'{or_mean[0]:.4f} (95% CI {ci[0][0]:.4f}-{ci[0][1]:.4f})'

record_ors.to_csv(f'{dn}/{outdir}/{datestring}_logistic_oddsratio_biomarkercategory.csv', index=True)

#### C. t-test comparing RCB score ####
record_t = pd.DataFrame(
    None, 
    index=['parity', 'recency (thres=5)', 'recency (thres=10)'],
    columns=['ER/PR+ HER2-', 'Triple Negative', 'HER2+']
)

for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    ## Parous vs. Nulliparous 
    parous = data.iloc[data.biomarker_subtypes.values==bm_cat,:].parous
    t, p = ttest_ind(
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].rcb.iloc[parous.values==1].values,
        data.iloc[data.biomarker_subtypes.values==bm_cat,:].rcb.iloc[parous.values==0].values,
        nan_policy='omit'
    )
    record_t.loc['parity', bm_cat] = f'{t:.4f} (p={p:.6f})'
    ## Recency 
    for thres in [5,10]:
        recent = data.iloc[((data.parous.values==1) & (data.biomarker_subtypes.values==bm_cat)),:].years_since_pregnancy < thres
        t, p = ttest_ind(
            data.iloc[((data.parous.values==1) & (data.biomarker_subtypes.values==bm_cat)),:].rcb.iloc[recent.values],
            data.iloc[((data.parous.values==1) & (data.biomarker_subtypes.values==bm_cat)),:].rcb.iloc[recent.values],
            nan_policy='omit'
        )
        record_t.loc[f'recency (thres={thres})', bm_cat] = f'{t:.4f} (p={p:.6f})'

record_t.to_csv(f'{dn}/{outdir}/{datestring}_t_test_biomarkercategory.csv')

#########################################################
#### Stratified analysis 2 - ER/PR positive vs. rest ####
#########################################################

#### A. Chi-squared/Fisher's exact test ####
record_chi = pd.DataFrame(
    None, 
    index=['<5 years', '5-10 years', '>=10 years'],
    columns=['ER/PR+', 'ER/PR-']
)

for erpr_cat in ['ER/PR+', 'ER/PR-']:
    crosstab = pd.crosstab(
        data.iloc[data.er_pr.values==erpr_cat,:].rcb_category.map(dd['response_nat']['Choices, Calculations, OR Slider Labels']),
        data.iloc[data.er_pr.values==erpr_cat,:].parity_category
    )[['Nulliparous', '<5 years', '5-10 years', '>=10 years']]
    # stattab = pd.concat((
    #     pd.DataFrame(crosstab.loc[['0','1']].sum(axis=0), columns=['Positive']).transpose(),
    #     pd.DataFrame(crosstab.loc[['2','3']].sum(axis=0), columns=['Negative']).transpose()
    # ))
    stattab = crosstab.drop('Unknown', axis=0)
    stattab.columns.name = None # remove column label for clarity 
    if '+' in erpr_cat:
        fn_bm_cat = 'pos'
    elif '-' in erpr_cat:
        fn_bm_cat = 'neg'
    #
    for par_query in ['<5 years', '5-10 years', '>=10 years']:
        f_obs = stattab.loc[:,['Nulliparous', par_query]]
        f_exp = np.array([f_obs.sum(axis=0)/2, f_obs.sum(axis=0)/2])
        chi, p = chisquare(f_obs=f_obs, f_exp=f_exp, axis=None)
        record_chi.loc[par_query, erpr_cat] = f'{chi:.2f} (p={p:.6f})'

record_chi.to_csv(f'{dn}/{outdir}/{datestring}_chi_test_erpr.csv')

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
        ].rcb_category.map(dd['response_nat']['Choices, Calculations, OR Slider Labels']).map(
            {'No residual disease' : 1, 'Residual disease' : 0, 'Unknown' : 0}
        ).fillna(0).to_numpy()
        # perform LR analysis 
        or_mean, ci = get_oddsratio_ci(X, y)
        record_ors.loc[par_query, erpr_cat] = f'{or_mean[0]:.4f} (95% CI {ci[0][0]:.4f}-{ci[0][1]:.4f})'

record_ors.to_csv(f'{dn}/{outdir}/{datestring}_logistic_oddsratio_erpr_category.csv', index=True)

#### C. t-test comparing RCB score ####
record_t = pd.DataFrame(
    None, 
    index=['parity', 'recency (thres=5)', 'recency (thres=10)'],
    columns=['ER/PR+', 'ER/PR-']
)

for erpr_cat in ['ER/PR+', 'ER/PR-']:
    ## Parous vs. Nulliparous 
    parous = data.iloc[data.er_pr.values==erpr_cat,:].parous
    t, p = ttest_ind(
        data.iloc[data.er_pr.values==erpr_cat,:].rcb.iloc[parous.values==1].values,
        data.iloc[data.er_pr.values==erpr_cat,:].rcb.iloc[parous.values==0].values,
        nan_policy='omit'
    )
    record_t.loc['parity', erpr_cat] = f'{t:.4f} (p={p:.6f})'
    ## Recency 
    for thres in [5,10]:
        recent = data.iloc[((data.parous.values==1) & (data.er_pr.values==erpr_cat)),:].years_since_pregnancy < thres
        t, p = ttest_ind(
            data.iloc[((data.parous.values==1) & (data.er_pr.values==erpr_cat)),:].rcb.iloc[recent.values],
            data.iloc[((data.parous.values==1) & (data.er_pr.values==erpr_cat)),:].rcb.iloc[recent.values],
            nan_policy='omit'
        )
        record_t.loc[f'recency (thres={thres})', erpr_cat] = f'{t:.4f} (p={p:.6f})'

record_t.to_csv(f'{dn}/{outdir}/{datestring}_t_test_erpr_category.csv')
