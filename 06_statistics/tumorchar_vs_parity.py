import os
import sys
import re
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import chisquare, ttest_ind, f_oneway
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler
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

data['nulliparous'] = np.invert(data['parous'].astype(bool)).astype(float)
data['parity <5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 5)
data['parity >=5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 5)
data['parity <10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 10)
data['parity >=10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 10)

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

feature_names=['er_pr_positive', 'her2_positive', 'triple_negative', 'grade_severe_2or3', 'her2_or_tn', 'grade_severe_3']

def generate_lrdata(df, parity_ref, parity_comp, feature_name):
    lrdata = pd.DataFrame(None, index=df.index, columns=[parity_comp, feature_name, 'age', 'fam_hx'])
    lrdata[feature_name] = df[feature_name].astype(float)
    lrdata['age'] = df['age_at_diagnosis'].astype(float)
    lrdata['fam_hx'] = df['fam_hx'].astype(float)
    for i in df.index:
        lrdata.loc[i,parity_comp] = 1 if df.loc[i,parity_comp]==1 else 0 if df.loc[i,parity_ref]==1 else np.nan
    # drop any rows with NaN
    lrdata.dropna(inplace=True, axis=0)
    # separate X and y
    X = lrdata[[parity_comp,'age','fam_hx']]
    y = lrdata[feature_name]
    # standard scale age
    scaler = StandardScaler()
    X['age'] = scaler.fit_transform(X['age'].values.reshape(-1,1))
    return X, y

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

datestring = datetime.now().strftime("%Y%m%d")

parity_comparisons = {
    'Parous vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parous'},
    '<5 vs. >=5 years' : {'ref' : 'parity >=5 years', 'comp' : 'parity <5 years'},
    '<10 vs. >=10 years' : {'ref' : 'parity >=10 years', 'comp' : 'parity <10 years'},
    '<5 vs. >=10 years' : {'ref' : 'parity >=10 years', 'comp' : 'parity <5 years'},
    '<5 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity <5 years'},
    '<10 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity <10 years'},
    '>=5 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity >=5 years'},
    '>=10 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity >=10 years'},
}

for parity_comparison in parity_comparisons.keys():
    print(f'#### {parity_comparison} ####')
    ref = parity_comparisons[parity_comparison]['ref']
    comp = parity_comparisons[parity_comparison]['comp']
    results = pd.DataFrame(columns=['Variable of interest', 'Category', comp, ref, 'OR (95% CI)', 'p-value'])
    for feature_name in feature_names:
        X, y = generate_lrdata(data, parity_ref=ref, parity_comp=comp, feature_name=feature_name)
        if ((np.sum((X[comp]==0) & (y==1))==0) | (np.sum((X[comp]==1) & (y==1))==0)):
            print(f'# Cannot perform comparison for {feature_name} due to lack of data\n')
        else:
            crosstab = pd.crosstab(X[comp], y)
            totals = crosstab.sum(axis=0)
            for i in crosstab.index:
                for j in crosstab.columns:
                    cts = crosstab.loc[i,j]
                    pcts = 100 * cts/totals[j]
                    crosstab.loc[i,j] = f'{cts} ({pcts:.1f}%)'
            oddsratios, cis, pvals = get_oddsratio_ci(X, y)
            or_formatted = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
            pval_formatted = f'{pvals[0]:.4f}'
            results = pd.concat((
                results, 
                pd.DataFrame({
                    'Variable of interest' : [feature_name, None],
                    'Category' : ['Yes', 'No'],
                    parity_comparisons[parity_comparison]['comp'] : [crosstab.loc[1,1], crosstab.loc[1,0]],
                    parity_comparisons[parity_comparison]['ref'] : [crosstab.loc[0,1], crosstab.loc[0,0]],
                    'OR (95% CI)' : [or_formatted, None],
                    'p-value' : [pval_formatted, None]
                })
            ))
    fout = datestring + '_tumchar_' + parity_comparison.replace(' ', '_') + '.csv'
    results.to_csv(f'{dout}/{fout}', index=False)
