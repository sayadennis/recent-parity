import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/06_exported_from_redcap'

#######################################
#### Read and lightly process data ####
#######################################

data = pd.read_csv(f'{dn}/FrequencyAndResultsO_DATA_2022-09-26_1703.csv')

######################################
#### Logistic regression analysis ####
######################################

genes=['any_patho_mutation', 'brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']

def generate_lrdata(df, genelist=genes, recency_thres=10):
    lrdata = pd.DataFrame(None, index=None, columns=['parous', 'has_mutation', 'recent', 'age'])
    for i in range(df.shape[0]):
        pat_dict = {}
        pat_dict['parous'] = [df['parous'].iloc[i]] # binary 0/1
        pat_dict['has_mutation'] = [float(df[genelist].iloc[i]==1)] # {'Present':1, 'Absent':2}
        pat_dict['recent'] = [int(df['years_since_pregnancy'].iloc[i] < recency_thres)]
        pat_dict['age'] = [df['age_at_diagnosis'].iloc[i]]
        lrdata = pd.concat((lrdata, pd.DataFrame(pat_dict)), ignore_index=True)
    # data for nulliparous vs. parous 
    X_np = np.array(lrdata[['parous', 'age']], dtype=float)
    X_np[:,1] = (X_np[:,1] - np.mean(X_np[:,1]))/np.std(X_np[:,1]) # standard scale
    y_np = np.array(lrdata['has_mutation'], dtype=float)
    # data for recency of parity --only select women who are parous
    X_rec = np.array(lrdata.iloc[lrdata['parous'].values==1][['recent', 'age']], dtype=float)
    X_rec[:,1] = (X_rec[:,1] - np.mean(X_rec[:,1]))/np.std(X_rec[:,1]) # standard scale
    y_rec = np.array(lrdata.iloc[lrdata['parous'].values==1]['has_mutation'], dtype=float)
    return X_np, y_np, X_rec, y_rec

def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    or1, or2 = [], []
    i = 0
    for i in range(rep):
        X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
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


######################
#### Perform test ####
######################

for gene in genes:
    print(f'######## Results for {gene} ########')
    X_np, y_np, X_rec, y_rec = generate_lrdata(data, genelist=[gene], recency_thres=10)
    #
    # if np.all(y_np==0):
    if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
        print('# Cannot perform parous vs. nulliparous comparison due to lack of mutation carriers\n')
    else:
        try:
            X_np_nonan=np.delete(X_np, np.where(np.isnan(X_np))[0], axis=0)
            y_np_nonan=np.delete(y_np, np.where(np.isnan(X_np))[0])
            oddsratios, cis = get_oddsratio_ci(X_np_nonan, y_np_nonan)
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
        try:
            X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
            y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
            oddsratios, cis = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
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
        try:
            X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
            y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
            oddsratios, cis = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
            print('\n#### recent vs non-recent (recency threshold 5 years) ####')
            print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
            print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
        except:
            print('Could not calculate odds ratio.\n')

