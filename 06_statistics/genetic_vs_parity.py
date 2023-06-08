import os
import sys
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
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

# mark family history as none for patients who are missing it 
data.fam_hx = data.fam_hx.fillna(0.)

#########################################
#### Fill in useful categorical info ####
#########################################

data['nulliparous'] = np.invert(data['parous'].astype(bool)).astype(float)
data['parity <5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 5)
data['parity >=5 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 5)
data['parity <10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] < 10)
data['parity >=10 years'] = np.where(np.isnan(data['years_since_pregnancy']), np.nan, data['years_since_pregnancy'] >= 10)

genes=['any_patho_mutation', 'brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']

data['any_patho_mutation'] = (data['any_patho_mutation'].map(dd['any_patho_mutation']['Choices, Calculations, OR Slider Labels'])=='Pathogenic').astype(int)

for feature_name in genes[1:]:
    data[feature_name] = (data[feature_name].map(dd[feature_name]['Choices, Calculations, OR Slider Labels'])=='Pathogenic').astype(int)

genes = [item for item in genes if data[item].sum(axis=0)>=2]

###################################
#### Define analysis functions ####
###################################

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
    results = pd.DataFrame(columns=['Variable of interest', 'Category', ref, comp, 'OR (95% CI)', 'p-value'])
    for feature_name in genes:
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
                    parity_comparisons[parity_comparison]['ref'] : [crosstab.loc[1,1], crosstab.loc[1,0]],
                    parity_comparisons[parity_comparison]['comp'] : [crosstab.loc[0,1], crosstab.loc[0,0]],
                    'OR (95% CI)' : [or_formatted, None],
                    'p-value' : [pval_formatted, None]
                })
            ))
    fout = datestring + '_genetic_' + parity_comparison.replace(' ', '_') + '.csv'
    results.to_csv(f'{dout}/{fout}', index=False)


# results_parity = pd.DataFrame(index=genes, columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value'])
# results_recency10 = pd.DataFrame(index=genes, columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value'])
# results_recency5 = pd.DataFrame(index=genes, columns=['varname', 'or', 'low', 'high', 'formatted', 'p-value'])

# for gene in genes:
#     print(f'######## Results for {gene} ########')
#     X_np, y_np, X_rec, y_rec = generate_lrdata(data, genelist=[gene], recency_thres=10)
#     #
#     # if np.all(y_np==0):
#     if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
#         print('# Cannot perform parous vs. nulliparous comparison due to lack of mutation carriers\n')
#     else:
#         try:
#             X_np_nonan = np.delete(X_np, np.where(np.isnan(X_np))[0], axis=0)
#             y_np_nonan = np.delete(y_np, np.where(np.isnan(X_np))[0])
#             oddsratios, cis, pvals = get_oddsratio_ci(X_np_nonan, y_np_nonan)
#             results_parity.loc[gene,'varname'] = gene
#             results_parity.loc[gene,'or'] = oddsratios[0]
#             results_parity.loc[gene,'low'] = cis[0][0]
#             results_parity.loc[gene,'high'] = cis[0][1]
#             results_parity.loc[gene,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
#             results_parity.loc[gene,'p-value'] = pvals[0]
#             # print('\n#### parous vs nulliparous ####')
#             # print(f'Odds ratio for parity: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
#             # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
#         except:
#             print(f'Could not calculate odds ratio for {gene} in parous vs. nulliparous.\n')
#     # 
#     # if np.all(y_rec==0):
#     if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
#         print('# Cannot perform recent vs. non-recent comparison with 10-year cutoff due to lack of mutation carriers in certain categories\n')
#     else:
#         try:
#             X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
#             y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
#             oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
#             results_recency10.loc[gene,'varname'] = gene
#             results_recency10.loc[gene,'or'] = oddsratios[0]
#             results_recency10.loc[gene,'low'] = cis[0][0]
#             results_recency10.loc[gene,'high'] = cis[0][1]
#             results_recency10.loc[gene,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
#             results_recency10.loc[gene,'p-value'] = pvals[0]
#             # print('\n#### recent vs non-recent (recency threshold 10 years) ####')
#             # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
#             # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
#         except:
#             print(f'Could not calculate odds ratio for {gene} in recency (10 years).\n')
#     # 
#     X_np, y_np, X_rec, y_rec = generate_lrdata(data, genelist=[gene], recency_thres=5)
#     # if np.all(y_rec==0):
#     if ((np.sum((X_np[:,0]==0) & (y_np==1))==0) | (np.sum((X_np[:,0]==1) & (y_np==1))==0)):
#         print('# Cannot perform recent vs. non-recent comparison with 5-year cutoff due to lack of mutation carriers in certain categories\n')
#     else:
#         try:
#             X_rec_nonan=np.delete(X_rec, np.where(np.isnan(X_rec))[0], axis=0)
#             y_rec_nonan=np.delete(y_rec, np.where(np.isnan(X_rec))[0])
#             oddsratios, cis, pvals = get_oddsratio_ci(X_rec_nonan, y_rec_nonan)
#             results_recency5.loc[gene,'varname'] = gene
#             results_recency5.loc[gene,'or'] = oddsratios[0]
#             results_recency5.loc[gene,'low'] = cis[0][0]
#             results_recency5.loc[gene,'high'] = cis[0][1]
#             results_recency5.loc[gene,'formatted'] = f'{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})'
#             results_recency5.loc[gene,'p-value'] = pvals[0]
#             # print('\n#### recent vs non-recent (recency threshold 5 years) ####')
#             # print(f'Odds ratio for recency: {oddsratios[0]:.4f} (95% CIs {cis[0][0]:.4f}-{cis[0][1]:.4f})')
#             # print(f'Odds ratio for age: {oddsratios[1]:.4f} (95% CIs {cis[1][0]:.4f}-{cis[1][1]:.4f})\n')
#         except:
#             print(f'Could not calculate odds ratio for {gene} in recency (5 years).\n')

# datestring = datetime.now().strftime("%Y%m%d")

# results_parity.to_csv(f'{dout}/{datestring}_mutation_vs_parity.csv', index=False)
# results_recency10.to_csv(f'{dout}/{datestring}_mutation_vs_recencyparity10.csv', index=False)
# results_recency5.to_csv(f'{dout}/{datestring}_mutation_vs_recencyparity5.csv', index=False)
