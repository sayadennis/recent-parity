import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

#######################################
#### Read and lightly process data ####
#######################################

data = pd.read_csv(f'{dn}/tkhr_pat_wo_gentest_added/recent_parity_combined_data_reviewed_tkhr_added.csv')

data['YearsSinceLastPreg']=data['AgeBCDx']-data['AgeLastPreg']

## Remove rows where 
data=data.iloc[~(data.YearsSinceLastPreg.values<0),:]

######################################
#### Logistic regression analysis ####
######################################

genes=['HasPathoMut', 'BRCA1_patho', 'BRCA2_patho', 'PALB2_patho', 'TP53_patho', 'PTEN_patho', 'CDH1_patho', 'STK11_patho', 'CHEK2_patho', 'ATM_patho']

def generate_lrdata(df, genelist=genes, recency_thres=10):
    data = pd.DataFrame(None, index=None, columns=['parous', 'has_mutation', 'recent', 'age'])
    for i in range(df.shape[0]):
        pat_dict = {}
        pat_dict['parous'] = [int(df['Para'].iloc[i]>0)]
        pat_dict['has_mutation'] = [int(np.any(df[genelist].iloc[i]==1))] # binary 0/1 
        agedx = df['AgeBCDx'].iloc[i]
        agelastpreg = df['AgeLastPreg'].iloc[i]
        pat_dict['recent'] = [int(agedx - agelastpreg < recency_thres)]
        pat_dict['age'] = [df['AgeBCDx'].iloc[i]]
        data = pd.concat((data, pd.DataFrame(pat_dict)), ignore_index=True)
    # 
    X_np = np.array(data[['parous', 'age']], dtype=int)
    X_np[:,1] = (X_np[:,1] - np.mean(X_np[:,1]))/np.std(X_np[:,1]) # standard scale
    y_np = np.array(data['has_mutation'], dtype=int)
    X_rec = np.array(data.iloc[data['parous'].values==1][['recent', 'age']], dtype=int)
    X_rec[:,1] = (X_rec[:,1] - np.mean(X_rec[:,1]))/np.std(X_rec[:,1]) # standard scale
    y_rec = np.array(data.iloc[data['parous'].values==1]['has_mutation'], dtype=int)
    return X_np, y_np, X_rec, y_rec

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    or1, or2 = [], []
    i = 0
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i)
        if ~np.all(y_bs==0):
            lrm = LogisticRegression(penalty='l2', solver='lbfgs')
            lrm.fit(X_bs, y_bs)
            or1.append(np.exp(lrm.coef_[0][0]))
            or2.append(np.exp(lrm.coef_[0][1]))
        else:
            continue
    oddsratio = (np.mean(or1), np.mean(or2))
    ci = ()
    # first get ci1
    p = ((1.0-alpha)/2.0) * 100
    lower = max(0.0, np.percentile(or1, p))
    p = (alpha+((1.0-alpha)/2.0)) * 100
    upper = np.percentile(or1, p)
    ci1 = (lower, upper)
    # next get ci2
    p = ((1.0-alpha)/2.0) * 100
    lower = max(0.0, np.percentile(or2, p))
    p = (alpha+((1.0-alpha)/2.0)) * 100
    upper = np.percentile(or2, p)
    ci2 = (lower, upper)
    ci = (ci1, ci2)
    return oddsratio, ci

#####################################
#### First include all mutations ####
#####################################
print('######## Results for all mutations combined ########')
X_np, y_np, X_rec, y_rec = generate_lrdata(data, recency_thres=10)

print('#### parous vs nulliparous ####')
# lrm = LogisticRegression(penalty='none', solver='lbfgs')
# lrm.fit(X_np, y_np)
if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
    print('# Cannot perform parous vs. nulliparous comparison due to lack of mutation carriers\n')
else:
    # lrm = LogisticRegression(penalty='none', solver='lbfgs')
    # lrm.fit(X_np, y_np)
    try:
        oddsratios, cis = get_oddsratio_ci(X_np, y_np)
        print(f'Odds ratio for parity: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
        print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
    except:
        print('Could not calculate odds ratio.\n')


print('#### recent vs non-recent (recency threshold 10 years) ####')
# lrm = LogisticRegression(penalty='none', solver='lbfgs')
# lrm.fit(X_rec, y_rec)
if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
    print('# Cannot perform recent vs. non-recent comparison with 10-year cutoff due to lack of mutation carriers in certain categories\n')
else:
    # lrm = LogisticRegression(penalty='none', solver='lbfgs')
    # lrm.fit(X_rec, y_rec)
    try:
        oddsratios, cis = get_oddsratio_ci(X_rec, y_rec)
        print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
        print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
    except:
        print('Could not calculate odds ratio.\n')


print('#### recent vs non-recent (recency threshold 5 years) ####')
X_np, y_np, X_rec, y_rec = generate_lrdata(data, recency_thres=5)
# lrm = LogisticRegression(penalty='none', solver='lbfgs')
# lrm.fit(X_rec, y_rec)
if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
    print('# Cannot perform recent vs. non-recent comparison with 5-year cutoff due to lack of mutation carriers in certain categories\n')
else:
    # lrm = LogisticRegression(penalty='none', solver='lbfgs')
    # lrm.fit(X_rec, y_rec)
    try:
        oddsratios, cis = get_oddsratio_ci(X_rec, y_rec)
        print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
        print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
    except:
        print('Could not calculate odds ratio.\n')


#########################################
#### Next perform test for each gene ####
#########################################

for gene in genes:
    print(f'######## Results for {gene} ########')
    X_np, y_np, X_rec, y_rec = generate_lrdata(data, genelist=[gene], recency_thres=10)
    #
    # if np.all(y_np==0):
    if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
        print('# Cannot perform parous vs. nulliparous comparison due to lack of mutation carriers\n')
    else:
        # lrm = LogisticRegression(penalty='none', solver='lbfgs')
        # lrm.fit(X_np, y_np)
        try:
            oddsratios, cis = get_oddsratio_ci(X_np, y_np)
            print('\n#### parous vs nulliparous ####')
            print(f'Odds ratio for parity: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')
    # 
    # if np.all(y_rec==0):
    if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
        print('# Cannot perform recent vs. non-recent comparison with 10-year cutoff due to lack of mutation carriers in certain categories\n')
    else:
        # lrm = LogisticRegression(penalty='none', solver='lbfgs')
        # lrm.fit(X_rec, y_rec)
        try:
            oddsratios, cis = get_oddsratio_ci(X_rec, y_rec)
            print('\n#### recent vs non-recent (recency threshold 10 years) ####')
            print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')
    # 
    X_np, y_np, X_rec, y_rec = generate_lrdata(data, genelist=[gene], recency_thres=5)
    # if np.all(y_rec==0):
    if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
        print('# Cannot perform recent vs. non-recent comparison with 5-year cutoff due to lack of mutation carriers in certain categories\n')
    else:
        # lrm = LogisticRegression(penalty='none', solver='lbfgs')
        # lrm.fit(X_rec, y_rec)
        try:
            oddsratios, cis = get_oddsratio_ci(X_rec, y_rec)
            print('\n#### recent vs non-recent (recency threshold 5 years) ####')
            print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')

