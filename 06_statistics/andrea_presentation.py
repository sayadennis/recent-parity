import os
import sys
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Read data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-07_1648.csv'

data = pd.read_csv(f'{dn}/{datadir}/{fn}')
data = data.iloc[data.nat.values==1,]

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

##########################
#### Define functions ####
##########################

def generate_lrdata(df, feature=[], confounder=[], target=[], recency_thres=10):
    lrdata = pd.DataFrame(None, index=None, columns=['parous', 'has_mutation', 'recent', 'age'])
    for i in range(df.shape[0]):
        pat_dict = {}
        pat_dict['parous'] = [df['parous'].iloc[i]] # binary 0/1
        pat_dict['recent'] = [int(df['years_since_pregnancy'].iloc[i] < recency_thres)]
        pat_dict['ER/PR'] = [(df['er_status'].iloc[i]==1 | df['pr_status'].iloc[i]==1)]
        pat_dict['target'] = [(df['rcb_category'].iloc[i]==1 | df['rcb_category'].iloc[i]==2)]
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

# def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
#     or1, or2 = [], []
#     i = 0
#     for i in range(rep):
#         X_bs, y_bs = resample(X, y, random_state=i) # create bootstrap (bs) sample
#         if ~np.all(y_bs==0):
#             lrm = LogisticRegression(penalty='l2', solver='lbfgs')
#             lrm.fit(X_bs, y_bs)
#             or1.append(np.exp(lrm.coef_[0][0]))
#             or2.append(np.exp(lrm.coef_[0][1]))
#         else:
#             continue
#     oddsratio = (np.mean(or1), np.mean(or2))
#     ci = ()
#     # first get ci1
#     p = ((1.0-alpha)/2.0) * 100
#     lower = max(0.0, np.percentile(or1, p))
#     p = (alpha+((1.0-alpha)/2.0)) * 100
#     upper = np.percentile(or1, p)
#     ci1 = (lower, upper)
#     # next get ci2
#     p = ((1.0-alpha)/2.0) * 100
#     lower = max(0.0, np.percentile(or2, p))
#     p = (alpha+((1.0-alpha)/2.0)) * 100
#     upper = np.percentile(or2, p)
#     ci2 = (lower, upper)
#     ci = (ci1, ci2)
#     return oddsratio, ci

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


# ###########################################################
# #### Stratified analysis 1 - three biological subtypes ####
# ###########################################################

# #### A. Chi-squared/Fisher's exact test ####

# #### B. Odds-ratio and confidence intervals calculated via logistic regression ####

# #########################################################
# #### Stratified analysis 2 - ER/PR positive vs. rest ####
# #########################################################

# #### A. Chi-squared/Fisher's exact test ####

# #### B. Odds-ratio and confidence intervals calculated via logistic regression ####

