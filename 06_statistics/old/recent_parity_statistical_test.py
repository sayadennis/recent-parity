import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

data = pd.read_csv(f'{dn}/recent_parity_combined_data_reviewed.csv')
reviewed = data.iloc[data['Review complete'].values=='x',:]

print('\nTotal number of patients included: %s\n' % reviewed.shape[0])

def calc_or(ct, recency_cutoff=None):
    if ct.shape[1]==2:
        if ((ct[0,1]==0) | (ct[1,0]==0) | (ct[0,0]==0)):
            oddsratio=0
        else:
            oddsratio = (ct[1,1]/ct[0,1])/(ct[1,0]/ct[0,0])
    elif ct.shape[1]==3:
        if recency_cutoff==5:
            if ((ct[0,1]+ct[0,2]==0) | (ct[1,0]==0) | (ct[0,0]==0)):
                oddsratio=0
            else:
                oddsratio = ((ct[1,1]+ct[1,2])/(ct[0,1]+ct[0,2])) / (ct[1,0]/ct[0,0])
        elif recency_cutoff==10:
            if ((ct[0,2]==0) | ((ct[1,0]+ct[1,1])==0) | ((ct[0,0]+ct[0,1])==0)):
                oddsratio=0
            else:
                oddsratio = (ct[1,2]/ct[0,2]) / ((ct[1,0]+ct[1,1])/(ct[0,0]+ct[0,1]))
    return oddsratio

def evaluate_recent_parity(df, genes, include_vus=False):
    #### First evaluate parous vs. nulliparous ####
    ct = pd.DataFrame(0, index=['no mutation', 'mutation'], columns=['nulliparous', 'parous'])
    # count patients 
    for i in range(df.shape[0]):
        parous = int(df['para'].iloc[i]>0) # binary 0/1 
        if include_vus:
            hasmut = int(
                np.any(df[genes].iloc[i]=='positive') | 
                np.any(df[genes].iloc[i]=='positive ') | 
                np.any(df[genes].iloc[i]=='variant') | 
                np.any(df[genes].iloc[i]=='variant ')
            ) # binary 0/1 
        else:
            hasmut = int(np.any(df[genes].iloc[i]=='positive') | np.any(df[genes].iloc[i]=='positive ')) # binary 0/1 
        ct.iloc[hasmut, parous] += 1
    # calculate odds ratio 
    oddsratio = calc_or(np.array(ct))
    print('## Results for parous vs. nulliparous ##')
    print('# Observed frequency:')
    print(ct)
    print('')
    print(f'# Odds-ratio: {oddsratio:.4f}')
    if np.all(ct>=5):
        # perform chi-square contingency test 
        chi2, p, dof, expected = chi2_contingency(ct)
        # report results
        print(f'# Chi-square test p-value: {p:.4f}\n\n')
        # print(pd.DataFrame(np.array([[chi2, p, dof]]), columns=['chi2stat', 'p-value', 'DOF']))
    else:
        _, p = fisher_exact(ct)
        print(f'# Fisher\'s exact test p-value: {p:.4f}\n\n')
        # print(pd.DataFrame(np.array([[oddsratio, p]]), columns=['odds ratio', 'p-value']))
        # print('Not enough patients to perform test:')
    df = df.iloc[~pd.isnull(df['ageoflastpreg']).values,:]
    ctdict_parnulli = {}
    ctdict_parnulli['no mut & nulliparous'] = ct.iloc[0,0]
    ctdict_parnulli['no mut & parous'] = ct.iloc[0,1]
    ctdict_parnulli['mutation & nulliparous'] = ct.iloc[1,0]
    ctdict_parnulli['mutation & parous'] = ct.iloc[1,1]
    ## Next perform a chi-square test based on recency of pregnancy ## 
    ct = pd.DataFrame(0, index=['no mutation', 'mutation'], columns=['10+', '5-9', '<5'])
    # count patients 
    for i in range(df.shape[0]):
        agedx = df['age_at_diagnosis'].iloc[i]
        agelastpreg = df['ageoflastpreg'].iloc[i]
        if (agedx-agelastpreg) < 5:
            recency = 2
        elif (agedx-agelastpreg) < 10:
            recency = 1
        else:
            recency = 0
        if include_vus:
            hasmut = int(
                np.any(df[genes].iloc[i]=='positive') | 
                np.any(df[genes].iloc[i]=='positive ') | 
                np.any(df[genes].iloc[i]=='variant') | 
                np.any(df[genes].iloc[i]=='variant ')
            ) # binary 0/1 
        else:
            hasmut = int(np.any(df[genes].iloc[i]=='positive') | np.any(df[genes].iloc[i]=='positive ')) # binary 0/1 
        ct.iloc[hasmut, recency] += 1
    ctdict_recency = {}
    ctdict_recency['no mut & 10+'] = ct.iloc[0,0]
    ctdict_recency['no mut & 5-9'] = ct.iloc[0,1]
    ctdict_recency['no mut & <5'] = ct.iloc[0,2]
    ctdict_recency['mutation & 10+'] = ct.iloc[1,0]
    ctdict_recency['mutation & 5-9'] = ct.iloc[1,1]
    ctdict_recency['mutation & <5'] = ct.iloc[1,2]
    # Calculate odds ratio
    ctarray = np.array(ct)
    oddsratio_5 = calc_or(ctarray, recency_cutoff=5)
    oddsratio_10 = calc_or(ctarray, recency_cutoff=10)
    print('## Results for recency of pregnancy ##')
    print('# Observed frequency:')
    print(ct)
    print('')
    print(f'# Odds-ratio with 5 year recency cutoff: {oddsratio_5:.4f}')
    print(f'# Odds-ratio with 10 year recency cutoff: {oddsratio_10:.4f}')
    # get two tables with different recency thresholds 
    thres5 = np.array(np.array([ctarray[:,0], ctarray[:,1] + ctarray[:,2]]).T)
    thres10 = np.array(np.array([ctarray[:,0] + ctarray[:,1], ctarray[:,2]]).T)
    # test for the 5 year cutoff 
    if np.all(thres5>5):
        # perform chi-square contingency test 
        chi2, p, dof, expected = chi2_contingency(thres5)
        # report results 
        print(f'# Chi-square test p-value with 5 year recency cutoff: {p:.4f}')
        # print(pd.DataFrame(np.array([[chi2, p, dof]]), columns=['chi2stat', 'p-value', 'DOF']))
    else:
        _, p = fisher_exact(thres5)
        print(f'# Fisher\'s exact test p-value with 5 year recency cutoff: {p:.4f}')
    # test for the 10 year cutoff 
    if np.all(thres10>5):
        # perform chi-square contingency test 
        chi2, p, dof, expected = chi2_contingency(thres10)
        # report results 
        print(f'# Chi-square test p-value with 10 year recency cutoff: {p:.4f}\n\n')
        # print(pd.DataFrame(np.array([[chi2, p, dof]]), columns=['chi2stat', 'p-value', 'DOF']))
    else:
        _, p = fisher_exact(thres10)
        print(f'# Fisher\'s exact test p-value with 10 year recency cutoff: {p:.4f}\n\n')
    return ctdict_parnulli, ctdict_recency

#############################################
#### Combine all mutation results in one ####
#############################################
print('######## Results for all mutations combined ########')

genes = ['general', 'brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

ctdict_parnulli, ctdict_recency = evaluate_recent_parity(reviewed, genes)

# save all counts to a separate table
allct_parnulli = pd.DataFrame(None, columns=['no mut & nulliparous', 'no mut & parous', 'mutation & nulliparous', 'mutation & parous'])
allct_recency = pd.DataFrame(None, columns=['no mut & <5', 'no mut & 5-9', 'no mut & 10+', 'mutation & <5', 'mutation & 5-9', 'mutation & 10+'])

allct_parnulli = pd.concat((allct_parnulli, pd.DataFrame(ctdict_parnulli, index=['general'])), sort=True)
allct_recency = pd.concat((allct_recency, pd.DataFrame(ctdict_recency, index=['general'])), sort=True)

##################################
#### Look at individual genes ####
##################################

for gene in genes[2:]:
    gene_list = [gene]
    print(f'######## Results for {gene} ########')
    ctdict_parnulli, ctdict_recency = evaluate_recent_parity(reviewed, gene_list)
    allct_parnulli = pd.concat((allct_parnulli, pd.DataFrame(ctdict_parnulli, index=[gene])), sort=True)
    allct_recency = pd.concat((allct_recency, pd.DataFrame(ctdict_recency, index=[gene])), sort=True)

allct_parnulli.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/all_patient_counts_parous_vs_nulliparous.csv', header=True, index=True)
allct_recency.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/all_patient_counts_recency_of_parity.csv', header=True, index=True)

#############################################################################################
#### Combine all mutation results in one with Variant of Uncertain Sifnificance included ####
#############################################################################################
print('######## Results all mutations combined (including VUS) ########')

genes = ['general', 'brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

ctdict_parnulli, ctdict_recency = evaluate_recent_parity(reviewed, genes, include_vus=True)

# save all counts to a separate table
allct_parnulli = pd.DataFrame(None, columns=['no mut & nulliparous', 'no mut & parous', 'mutation & nulliparous', 'mutation & parous'])
allct_recency = pd.DataFrame(None, columns=['no mut & <5', 'no mut & 5-9', 'no mut & 10+', 'mutation & <5', 'mutation & 5-9', 'mutation & 10+'])

allct_parnulli = pd.concat((allct_parnulli, pd.DataFrame(ctdict_parnulli, index=['general'])), sort=True)
allct_recency = pd.concat((allct_recency, pd.DataFrame(ctdict_recency, index=['general'])), sort=True)

##################################################################################
#### Look at individual genes with Variant of Uncertain Sifnificance included ####
##################################################################################

for gene in genes[2:]:
    print('')
    gene_list = [gene]
    print(f'######## Results for {gene} (including VUS) ########')
    ctdict_parnulli, ctdict_recency = evaluate_recent_parity(reviewed, gene_list, include_vus=True)
    allct_parnulli = pd.concat((allct_parnulli, pd.DataFrame(ctdict_parnulli, index=[gene])), sort=True)
    allct_recency = pd.concat((allct_recency, pd.DataFrame(ctdict_recency, index=[gene])), sort=True)

allct_parnulli.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/all_patient_counts_parous_vs_nulliparous_withVUS.csv', header=True, index=True)
allct_recency.to_csv('/share/fsmresfiles/breast_cancer_pregnancy/data/all_patient_counts_recency_of_parity_withVUS.csv', header=True, index=True)
