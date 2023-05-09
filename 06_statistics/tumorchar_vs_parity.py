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
dout = '/share/fsmresfiles/breast_cancer_pregnancy/stat_results'

###################
#### Read data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'

data = pd.read_csv(f'{dn}/{datadir}/{fn}')

# remove exluded patients
data = data.iloc[(data.exclude_demo.values!=1) & (data.exclude_tum.values!=1),:]

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

#########################################
#### Fill in useful categorical info ####
#########################################

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
data['er_pr_positive'] = None
data['her2_positive'] = None
data['her2_or_tn'] = None
data['triple_negative'] = None
data['grade_severe_3'] = (data['histologic_grade'].map(dd['histologic_grade']['Choices, Calculations, OR Slider Labels'])=='3').astype(int)
data['grade_severe_2or3'] = ((data['histologic_grade'].map(dd['histologic_grade']['Choices, Calculations, OR Slider Labels'])=='3') | (data['histologic_grade'].map(dd['histologic_grade']['Choices, Calculations, OR Slider Labels'])=='2')).astype(int)

missing_ix = []
for i in data.index:
    ## Get ER data and set to None if missing 
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
    except:
        er = None
    ## Researt for PR 
    try:
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        pr = None
    ## Researt for HER2 
    try:
        her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    except:
        her2 = None
    ## Now record information in the main dataframe 
    if np.all(~pd.isnull([er, pr, her2])): # if all information is available
        if (np.any([er, pr]) & (not her2)):
            data.loc[i,'biomarker_subtypes'] = 'ER/PR+ HER2-'
            data.loc[i,'triple_negative'] = 0
        elif ((not er) & (not pr) & (not her2)):
            data.loc[i,'biomarker_subtypes'] = 'Triple Negative'
            data.loc[i,'triple_negative'] = 1
        elif her2:
            data.loc[i,'biomarker_subtypes'] = 'HER2+'
            data.loc[i,'triple_negative'] = 0
        else:
            print('Unknown pattern at index:', i)
    # 
    if np.any(~pd.isnull([er, pr])): # if one of ER/PR data is available
        if np.any([er, pr]):
            data.loc[i,'er_pr_positive'] = 1
        else:
            data.loc[i,'er_pr_positive'] = 0
    # 
    if ~pd.isnull(her2): # If HER2 data is available 
        if her2:
            data.loc[i,'her2_positive'] = 1
        else:
            data.loc[i,'her2_positive'] = 0

data['her2_or_tn'] = (data['her2_positive'] | data['triple_negative']).astype(int)

#################################################################
#### Assess the effects of potentially confounding variables ####
#################################################################

# one-way ANOVA for continuous variables 
# chi-square/Fisher’s Exact tests for categorical variables 

## Tumor size 
sizes = {par_cat : data.iloc[data.parity_category.values==par_cat,:].tumor_size.dropna().values 
        for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']}
stat, p = f_oneway(sizes['Nulliparous'], sizes['<5 years'], sizes['5-10 years'], sizes['>=10 years'])

print('#### Statistical test for tumor size ####')
print(f'One-way ANOVA f={stat:.4f} (p={p:.4f})\n')

## Tumor grade ### We need to determine how to test this 
crosstab = pd.crosstab(data['histologic_grade'], data['parity_category']).rename(
    dd['histologic_grade']['Choices, Calculations, OR Slider Labels'],
    axis=0
)
stat, p = chisquare(crosstab, axis=0)

print('#### Statistical test for tumor histologic grade ####')
print(crosstab)
for i, parity_category in enumerate(crosstab.columns):
    print(f'Chi-squared test for {parity_category} chi={stat[i]:.2e} (p={p[i]:.2e})')

##########################################################
#### Assess the associations between tumor and parity ####
##########################################################

feature_names=['her2_positive', 'triple_negative', 'er_pr_positive', 'her2_or_tn', 'grade_severe_3', 'grade_severe_2or3']

def generate_lrdata(df, feature_name, recency_thres=10):
    lrdata = pd.DataFrame(None, index=None, columns=['parous', feature_name, 'recent', 'age', 'fam_hx'])
    for i in range(df.shape[0]):
        pat_dict = {}
        pat_dict['parous'] = [df['parous'].iloc[i]] # binary 0/1
        pat_dict[feature_name] = [df[feature_name].iloc[i]]
        pat_dict['recent'] = [int(df['years_since_pregnancy'].iloc[i] < recency_thres)]
        pat_dict['age'] = [df['age_at_diagnosis'].iloc[i]]
        pat_dict['fam_hx'] = [df['fam_hx'].iloc[i]]
        if np.any(pd.isnull(list(pat_dict.values()))):
            continue
        lrdata = pd.concat((lrdata, pd.DataFrame(pat_dict)), ignore_index=True)
    # data for nulliparous vs. parous 
    X_np = np.array(lrdata[['parous', 'age', 'fam_hx']], dtype=float)
    X_np[:,1] = (X_np[:,1] - np.mean(X_np[:,1]))/np.std(X_np[:,1]) # standard scale
    y_np = np.array(lrdata[feature_name], dtype=float)
    # data for recency of parity --only select women who are parous
    X_rec = np.array(lrdata.iloc[lrdata['parous'].values==1][['recent', 'age', 'fam_hx']], dtype=float)
    X_rec[:,1] = (X_rec[:,1] - np.mean(X_rec[:,1]))/np.std(X_rec[:,1]) # standard scale
    y_rec = np.array(lrdata.iloc[lrdata['parous'].values==1][feature_name], dtype=float)
    return X_np, y_np, X_rec, y_rec

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    or1, or2, or3 = [], [], []
    i = 0
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
        if ~np.all(y_bs==0):
            lrm = LogisticRegression(penalty='l2', solver='lbfgs')
            lrm.fit(X_bs, y_bs)
            or1.append(np.exp(lrm.coef_[0][0]))
            or2.append(np.exp(lrm.coef_[0][1]))
            or3.append(np.exp(lrm.coef_[0][2]))
        else:
            continue
    oddsratios = [np.mean(or1), np.mean(or2), np.mean(or3)]
    # first get ci1
    ci_lower = ((1.0-alpha)/2.0) * 100
    ci_higher = (alpha+((1.0-alpha)/2.0)) * 100
    ci, pvals = [], []
    for bs_sample in [or1, or2, or3]:
        lower = max(0.0, np.percentile(bs_sample, ci_lower))
        upper = np.percentile(bs_sample, ci_higher)
        ci.append((lower, upper))
        pvals.append(np.min([(np.array(bs_sample)<1).mean(), (np.array(bs_sample)>1).mean()])*2)
    return oddsratios, ci, pvals

##########################
#### Perform Analysis ####
##########################

results_parity = pd.DataFrame(index=feature_names, columns=['varname', 'or', 'low', 'high'])
results_recency5 = pd.DataFrame(index=feature_names, columns=['varname', 'or', 'low', 'high'])
results_recency10 = pd.DataFrame(index=feature_names, columns=['varname', 'or', 'low', 'high'])

for feature_name in feature_names:
    print(f'######## Results for {feature_name} ########')
    X_np, y_np, X_rec, y_rec = generate_lrdata(data, feature_name=feature_name, recency_thres=10)
    #
    # if np.all(y_np==0):
    if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
        print('# Cannot perform parous vs. nulliparous comparison due to lack of data\n')
    else:
        try:
            X_np_nonan = np.delete(X_np, np.where(np.isnan(X_np))[0], axis=0)
            y_np_nonan = np.delete(y_np, np.where(np.isnan(X_np))[0])
            oddsratios, cis, pvals = get_oddsratio_ci(X_np_nonan, y_np_nonan)
            results_parity.loc[feature_name,'varname'] = feature_name
            results_parity.loc[feature_name,'or'] = oddsratios[0]
            results_parity.loc[feature_name,'low'] = cis[0][0]
            results_parity.loc[feature_name,'high'] = cis[0][1]
            results_parity.loc[feature_name,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
            results_parity.loc[feature_name,'pval'] = pvals[0]
            # print('\n#### parous vs nulliparous ####')
            # print(f'Odds ratio for parity: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')
    # 
    # if np.all(y_rec==0):
    if ((np.sum((X_rec[:,0]==0) & (y_rec==1))==0) | (np.sum((X_rec[:,0]==1) & (y_rec==1))==0)):
        print('# Cannot perform recent vs. non-recent comparison with 10-year cutoff due to lack of data\n')
    else:
        try:
            X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
            y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
            oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
            results_recency10.loc[feature_name,'varname'] = feature_name
            results_recency10.loc[feature_name,'or'] = oddsratios[0]
            results_recency10.loc[feature_name,'low'] = cis[0][0]
            results_recency10.loc[feature_name,'high'] = cis[0][1]
            results_recency10.loc[feature_name,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
            results_recency10.loc[feature_name,'pval'] = pvals[0]
            # print('\n#### recent vs non-recent (recency threshold 10 years) ####')
            # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')
    # 
    X_np, y_np, X_rec, y_rec = generate_lrdata(data, feature_name=feature_name, recency_thres=5)
    # if np.all(y_rec==0):
    if ((np.sum((X_rec[:,0]==0) & (y_rec==1))==0) | (np.sum((X_rec[:,0]==1) & (y_rec==1))==0)):
        print('# Cannot perform recent vs. non-recent comparison with 5-year cutoff due to lack of data\n')
    else:
        try:
            X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
            y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
            oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
            results_recency5.loc[feature_name,'varname'] = feature_name
            results_recency5.loc[feature_name,'or'] = oddsratios[0]
            results_recency5.loc[feature_name,'low'] = cis[0][0]
            results_recency5.loc[feature_name,'high'] = cis[0][1]
            results_recency5.loc[feature_name,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
            results_recency5.loc[feature_name,'pval'] = pvals[0]
            # print('\n#### recent vs non-recent (recency threshold 5 years) ####')
            # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')

datestring = datetime.now().strftime("%Y%m%d")

results_parity.to_csv(f'{dout}/{datestring}_tumchar_vs_parity.csv', index=False)
results_recency10.to_csv(f'{dout}/{datestring}_tumchar_vs_recencyparity10.csv', index=False)
results_recency5.to_csv(f'{dout}/{datestring}_tumchar_vs_recencyparity5.csv', index=False)
