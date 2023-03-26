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

data['pCR'] = (data.response_nat.map(dd['response_nat']['Choices, Calculations, OR Slider Labels'])=='No residual disease').astype(int)
data['positive'] = ((
    (data.rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels'])=='0').values |
    (data.rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels'])=='1').values
)).astype(int)

##########################
#### Define functions ####
##########################

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

##############################################################################
#### Stratified analysis for pCR and positive response (RCB category 0-1) ####
##############################################################################

results = {
    'ER/PR+ HER2-' : {
        'pCR' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
        'positive' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
    }, 
    'Triple Negative' : {
        'pCR' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
        'positive' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
    }, 
    'HER2+' : {
        'pCR' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
        'positive' : pd.DataFrame(index=['parity', 'recency (thres=5)', 'recency (thres=10)'],columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value']),
    },
}

for bm_cat in ['ER/PR+ HER2-', 'Triple Negative', 'HER2+']:
    for outcome in ['pCR', 'positive']:
        subdata = data.iloc[data.biomarker_subtypes.values==bm_cat,:]
        X_np, y_np, X_rec, y_rec = generate_lrdata(subdata, feature_name=outcome, recency_thres=5)
        if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
            print('# Cannot perform parous vs. nulliparous comparison due to lack of mutation carriers\n')
        else:
            try:
                X_np_nonan = np.delete(X_np, np.where(np.isnan(X_np))[0], axis=0)
                y_np_nonan = np.delete(y_np, np.where(np.isnan(X_np))[0])
                oddsratios, cis, pvals = get_oddsratio_ci(X_np_nonan, y_np_nonan)
                results[bm_cat][outcome].loc['parity','varname'] = 'parity'
                results[bm_cat][outcome].loc['parity','or'] = oddsratios[0]
                results[bm_cat][outcome].loc['parity','low'] = cis[0][0]
                results[bm_cat][outcome].loc['parity','high'] = cis[0][1]
                results[bm_cat][outcome].loc['parity','formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
                results[bm_cat][outcome].loc['parity','p-value'] = pvals[0]
                # print('\n#### parous vs nulliparous ####')
                # print(f'Odds ratio for parity: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
                # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
            except:
                print(f'Could not calculate odds ratio for {bm_cat} and {outcome}.\n')
        #
        if ((np.sum((X_rec[:,0]==0) & (y_rec==1))==0) | (np.sum((X_rec[:,0]==1) & (y_rec==1))==0)):
            print('# Cannot perform recent vs. non-recent comparison with 5-year cutoff due to lack of data\n')
        else:
            try:
                X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
                y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
                oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
                results[bm_cat][outcome].loc['recency (thres=5)','varname'] = 'recency (thres=5)'
                results[bm_cat][outcome].loc['recency (thres=5)','or'] = oddsratios[0]
                results[bm_cat][outcome].loc['recency (thres=5)','low'] = cis[0][0]
                results[bm_cat][outcome].loc['recency (thres=5)','high'] = cis[0][1]
                results[bm_cat][outcome].loc['recency (thres=5)','formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
                results[bm_cat][outcome].loc['recency (thres=5)','p-value'] = pvals[0]
                # print('\n#### recent vs non-recent (recency threshold 10 years) ####')
                # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
                # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
            except:
                print(f'Could not calculate odds ratio for {bm_cat} and {outcome}.\n')
        # 
        X_np, y_np, X_rec, y_rec = generate_lrdata(subdata, feature_name=outcome, recency_thres=10)
        # if np.all(y_rec==0):
        if ((np.sum((X_rec[:,0]==0) & (y_rec==1))==0) | (np.sum((X_rec[:,0]==1) & (y_rec==1))==0)):
            print('# Cannot perform recent vs. non-recent comparison with 10-year cutoff due to lack of data\n')
        else:
            try:
                X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
                y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
                oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
                results[bm_cat][outcome].loc['recency (thres=10)','varname'] = 'recency (thres=10)'
                results[bm_cat][outcome].loc['recency (thres=10)','or'] = oddsratios[0]
                results[bm_cat][outcome].loc['recency (thres=10)','low'] = cis[0][0]
                results[bm_cat][outcome].loc['recency (thres=10)','high'] = cis[0][1]
                results[bm_cat][outcome].loc['recency (thres=10)','formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
                results[bm_cat][outcome].loc['recency (thres=10)','p-value'] = pvals[0]
                # print('\n#### recent vs non-recent (recency threshold 5 years) ####')
                # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
                # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
            except:
                print(f'Could not calculate odds ratio for {bm_cat} and {outcome}.\n')

datestring = datetime.now().strftime("%Y%m%d")

results['ER/PR+ HER2-']['pCR'].to_csv(f'{dout}/{datestring}_nat_erpr_pCR.csv')
results['ER/PR+ HER2-']['positive'].to_csv(f'{dout}/{datestring}_nat_erpr_positive_response.csv')
results['HER2+']['pCR'].to_csv(f'{dout}/{datestring}_nat_her2_pCR.csv')
results['HER2+']['positive'].to_csv(f'{dout}/{datestring}_nat_her2_positive_response.csv')
results['Triple Negative']['pCR'].to_csv(f'{dout}/{datestring}_nat_TN_pCR.csv')
results['Triple Negative']['positive'].to_csv(f'{dout}/{datestring}_nat_TN_positive_response.csv')

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
