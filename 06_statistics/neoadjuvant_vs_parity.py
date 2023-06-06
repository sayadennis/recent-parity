import os
import sys
import re
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import chisquare, ttest_ind
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

data['nulliparous'] = np.invert(data['parous'].astype(bool)).astype(float)
data['parity <5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 5)
data['parity >=5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 5)
data['parity <10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 10)
data['parity >=10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 10)

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

data['pCR'] = (data.response_nat.map(dd['response_nat']['Choices, Calculations, OR Slider Labels'])=='No residual disease').astype(int)
data['positive_response'] = ((
    (data.rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels'])=='0').values |
    (data.rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels'])=='1').values
)).astype(int)

##########################
#### Define functions ####
##########################

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
    ci, pvals = [], []
    ci_lower = ((1.0-alpha)/2.0) * 100
    ci_higher = (alpha+((1.0-alpha)/2.0)) * 100
    for colnum in range(X.shape[1]):
        # first get ci1
        lower = max(0.0, np.percentile(oddsratio[colnum], ci_lower))
        upper = np.percentile(oddsratio[colnum], ci_higher)
        ci.append((lower, upper))
        pvals.append(np.min([(np.array(oddsratio[colnum])<1).mean(), (np.array(oddsratio[colnum])>1).mean()])*2)
    return mean_or, ci, pvals

##########################
#### Perform Analysis ####
##########################

datestring = datetime.now().strftime("%Y%m%d")

parity_comparisons = {
    'Parous vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parous'},
    '<5 vs. >=5 years' : {'ref' : 'parity <5 years', 'comp' : 'parity >=5 years'},
    '<10 vs. >=10 years' : {'ref' : 'parity <10 years', 'comp' : 'parity >=10 years'},
    '<5 vs. >=10 years' : {'ref' : 'parity <5 years', 'comp' : 'parity >=10 years'},
    '<5 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity <5 years'},
    '<10 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity <10 years'},
    '>=5 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity >=5 years'},
    '>=10 years vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parity >=10 years'},
}

stratifications = data.biomarker_subtypes.unique()
feature_names = ['positive_response', 'pCR']

for parity_comparison in parity_comparisons.keys():
    print(f'#### {parity_comparison} ####')
    ref = parity_comparisons[parity_comparison]['ref']
    comp = parity_comparisons[parity_comparison]['comp']
    results = pd.DataFrame(columns=['Variable of interest', 'Stratification', 'Category', ref, comp, 'OR (95% CI)', 'p-value'])
    for feature_name in feature_names:
        for stratification in stratifications:
            subdata = data.iloc[data.biomarker_subtypes.values==stratification,:]
            X, y = generate_lrdata(subdata, parity_ref=ref, parity_comp=comp, feature_name=feature_name)
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
                        'Stratification' : [stratification, None],
                        'Category' : ['Yes', 'No'],
                        parity_comparisons[parity_comparison]['ref'] : [crosstab.loc[1,1], crosstab.loc[1,0]],
                        parity_comparisons[parity_comparison]['comp'] : [crosstab.loc[0,1], crosstab.loc[0,0]],
                        'OR (95% CI)' : [or_formatted, None],
                        'p-value' : [pval_formatted, None]
                    })
                ))
    fout = datestring + '_nat_' + parity_comparison.replace(' ', '_') + '.csv'
    results.to_csv(f'{dout}/{fout}', index=False)

#############################
#### t-test of RCB score ####
#############################

print('\n## Parous vs. Nulliparous Comparison ##')
# Non-stratified 
t, p = ttest_ind(data.rcb.iloc[data.parous.values==1], data.rcb.iloc[data.parous.values==0])
print(f'Overall: t={t:.2f} (p={p:.4f})')

# Stratified 
for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    subdata = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
    t, p = ttest_ind(subdata.rcb.iloc[subdata.parous.values==1], subdata.rcb.iloc[subdata.parous.values==0])
    print(f'{bm_cat}: t={t:.2f} (p={p:.4f})')

print('\n\n## Recent Parity vs. Non-Recent Parity Comparison ##')
for thres in [5,10]:
    print(f'\n# {thres}-year recency threshold')
    t, p = ttest_ind(
        data.rcb.iloc[data.years_since_pregnancy.values<thres], 
        data.rcb.iloc[data.years_since_pregnancy.values>=thres]
    )
    print(f'Overall: t={t:.2f} (p={p:.4f})')
    # Stratified 
    for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
        subdata = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
        t, p = ttest_ind(
            subdata.rcb.iloc[subdata.years_since_pregnancy.values<thres], 
            subdata.rcb.iloc[subdata.years_since_pregnancy.values>=thres]
        )
        print(f'{bm_cat}: t={t:.2f} (p={p:.4f})')
