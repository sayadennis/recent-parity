from datetime import datetime
import pickle

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import chi2_contingency
from scipy.stats import ttest_ind

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

## Load data & data dictionary
datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# remove exluded patients
data = data.iloc[(data.exclude_demo.values!=1) & (data.exclude_tum.values!=1),:]

## Fill in biomarker subtypes 
data['biomarker_subtypes'] = None
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
    except:
        er = np.nan
    try:
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        pr = np.nan
    if not pd.isnull(data.loc[i,'her2_status']):
        her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    else:
        her2 = np.nan
    #   
    if ((er==True) | (pr==True)):
        erpr = True
    else:
        erpr = False
    #
    if (np.all(np.isnan([er, pr])) | np.isnan(her2)):
        data.loc[i,'biomarker_subtypes'] = None
    elif ((erpr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'ER/PR+ HER2-'
    elif ((not erpr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'Triple Negative'
    elif her2:
        data.loc[i,'biomarker_subtypes'] = 'HER2+'
    else:
        print('Unknown pattern at index:', i)

data['er_pr'] = None
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
    except:
        er = np.nan
    try:
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        pr = np.nan
    if ((er==True) | (pr==True)):
        data.loc[i,'er_pr'] = 'ER/PR+'
    elif np.all(np.isnan([er, pr])):
        data.loc[i,'er_pr'] = None
    else:
        data.loc[i,'er_pr'] = 'ER/PR-'

#####################
#### NAC summary ####
#####################

## Biologic subtype
nac_biomarker_crosstab = pd.crosstab(data.biomarker_subtypes, data.nat)
nac_biomarker_crosstab.columns = nac_biomarker_crosstab.columns.map({0: 'No NAC', 1: 'Had NAC'})
nac_biomarker_crosstab.columns.name = None; nac_biomarker_crosstab.index.name = None
nac_biomarker_crosstab_pct = nac_biomarker_crosstab.copy()
for ix in nac_biomarker_crosstab.index:
    for colname in nac_biomarker_crosstab.columns:
        ct = nac_biomarker_crosstab.loc[ix,colname]
        pct = 100 * nac_biomarker_crosstab.loc[ix,colname] / nac_biomarker_crosstab.loc[ix,:].sum()
        nac_biomarker_crosstab_pct.loc[ix,colname] = f'{ct} ({pct:.1f}%)'

nac_biomarker_crosstab_pct.to_csv(f'{dn}/summary_tables/nac_by_biomarker.csv', index=True, header=True)

chi2, p, dof, expected = chi2_contingency(nac_biomarker_crosstab)
fisher_exact(nac_biomarker_crosstab)

## Tumor size
fig, ax = plt.subplots(figsize=(8,5))

nonac_patho = ax.violinplot(data.iloc[data.nat.values!=1,:]['tumor_size'].dropna().values, positions=[1], showmedians=True)
nac_clin = ax.violinplot(data.iloc[data.nat.values==1,:]['clin_tumor_size'].dropna().values, positions=[2], showmedians=True)
nac_patho = ax.violinplot(data.iloc[data.nat.values==1,:]['tumor_size'].dropna().values, positions=[3], showmedians=True)

# Customize the violin plot
ax.set_xticks([1, 2, 3])  # Position of the columns on the x-axis
ax.set_xticklabels(["No NAC\n(pathological)", "NAC\n(clinical)", "NAC\n(pathological)"])
ax.set_ylabel("Tumor size (cm)")
ax.set_title("Tumor size distributions by NAC category")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
fig.savefig(f'{dn}/plots/tumor_size_by_nac.png')
plt.close()

t, p = ttest_ind(
    data.iloc[data.nat.values!=1,:]['tumor_size'].dropna().values,
    data.iloc[data.nat.values==1,:]['clin_tumor_size'].dropna().values
)
print(f'Tumor size comparison between no NAC (pathological size) and NAC (clinical size): t={t:.2f} (p={p:.2e})')

t, p = ttest_ind(
    data.iloc[data.nat.values==1,:]['clin_tumor_size'].dropna().values,
    data.iloc[data.nat.values==1,:]['tumor_size'].dropna().values
)
print(f'Tumor size comparison between NAC (clinical size) and NAC (pathological size): t={t:.2f} (p={p:.2e})')

## Nodal involvement

# Distribution among NAC patients
nac_clinical = data.iloc[data.nat.values==1,:].clin_node_stag_cat.map(dd['clin_node_stag_cat']['Choices, Calculations, OR Slider Labels']).value_counts()
nac_patho = data.iloc[data.nat.values==1,:].node_staging_category.map(dd['node_staging_category']['Choices, Calculations, OR Slider Labels']).value_counts()
# Distribution among all patients
nonac_patho = data.iloc[data.nat.values!=1,:].node_staging_category.map(dd['node_staging_category']['Choices, Calculations, OR Slider Labels']).value_counts()

nodal_involvement = pd.DataFrame(index=['N0', 'N1-3', 'NX'], columns=['No NAC (patho)', 'Had NAC (clinical)', 'Had NAC (patho)'])

for label, df in zip(['No NAC (patho)', 'Had NAC (clinical)', 'Had NAC (patho)'], [nonac_patho, nac_clinical, nac_patho]):
    nodal_involvement.loc['N0', label] = df.loc[['N0' in x for x in df.index]].sum()
    nodal_involvement.loc['N1-3', label] = df.iloc[[('N1' in x) | ('N2' in x) | ('N3' in x) for x in df.index]].sum()
    nodal_involvement.loc['NX', label] = df.iloc[['NX' in x for x in df.index]].sum()

nodal_involvement_pct = nodal_involvement.copy()
for i in nodal_involvement.index:
    for j in nodal_involvement.columns:
        ct = nodal_involvement.loc[i,j]
        pct = 100 * ct / nodal_involvement.loc[:,j].sum(axis=0)
        nodal_involvement_pct.loc[i,j] = f'{ct} ({pct:.1f}%)'

for label in ['No NAC (patho)', 'Had NAC (clinical)', 'Had NAC (patho)']:
    total_ct = nodal_involvement.loc[:,label].sum(axis=0)
    nodal_involvement_pct.loc['Total', label] = f'{total_ct} (100.0%)'

nodal_involvement_pct.to_csv(f'{dn}/summary_tables/nodal_involvement_nac_vs_not.csv')

## Age
fig, ax = plt.subplots(figsize=(6,4))

nonac_age = ax.violinplot(data.iloc[data.nat.values!=1,:]['age_at_diagnosis'].dropna().values, positions=[1], showmedians=True)
nac_age = ax.violinplot(data.iloc[data.nat.values==1,:]['age_at_diagnosis'].dropna().values, positions=[2], showmedians=True)

# Customize the violin plot
ax.set_xticks([1, 2])  # Position of the columns on the x-axis
ax.set_xticklabels(["No NAC", "NAC"])
ax.set_ylabel("Age at diagnosis")
ax.set_title("Age distributions by NAC category")

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
fig.savefig(f'{dn}/plots/age_by_nac.png')
plt.close()

t, p = ttest_ind(
    data.iloc[data.nat.values!=1,:]['age_at_diagnosis'].dropna().values,
    data.iloc[data.nat.values==1,:]['age_at_diagnosis'].dropna().values
)
print(f'Age comparison: t={t:.2f} (p={p:.2e})')

# Age separated by NAC status and biologic subtype! 
fig, ax = plt.subplots(figsize=(6,4))

pos = [1,2,4,5,7,8]
violins = ax.violinplot(
    [data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='ER/PR+ HER2-'),:]['age_at_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='ER/PR+ HER2-'),:]['age_at_diagnosis'].dropna().values,
    data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='HER2+'),:]['age_at_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='HER2+'),:]['age_at_diagnosis'].dropna().values,
    data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='Triple Negative'),:]['age_at_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='Triple Negative'),:]['age_at_diagnosis'].dropna().values],
    pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [0,2,4]:
    violins['bodies'][i].set_facecolor('blue')

for i in [1,3,5]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([1.5, 4.5, 7.5])
ax.set_xticklabels(['ER/PR+ HER2-', 'HER2+', 'Triple Negative'], ha='right', rotation=30)
ax.set_ylabel('Age at diagnosis')
ax.set_title('Diagnosis age by NAC status and biologic subtype')
ax.legend(['No NAC', 'Had NAC'], loc='lower right')

plt.tight_layout()
fig.savefig(f'{dn}/plots/age_by_nac_and_biomarker.png')
plt.close()

# Year of diagnosis by NAC status and biologic subtype

fig, ax = plt.subplots(figsize=(6,4))

data = data.iloc[data.year_of_diagnosis.values>=2010,:]

pos = [1,2,4,5,7,8]
violins = ax.violinplot(
    [data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='ER/PR+ HER2-'),:]['year_of_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='ER/PR+ HER2-'),:]['year_of_diagnosis'].dropna().values,
    data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='HER2+'),:]['year_of_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='HER2+'),:]['year_of_diagnosis'].dropna().values,
    data.iloc[(data.nat.values!=1) & (data.biomarker_subtypes.values=='Triple Negative'),:]['year_of_diagnosis'].dropna().values,
    data.iloc[(data.nat.values==1) & (data.biomarker_subtypes.values=='Triple Negative'),:]['year_of_diagnosis'].dropna().values],
    pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [0,2,4]:
    violins['bodies'][i].set_facecolor('blue')

for i in [1,3,5]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([1.5, 4.5, 7.5])
ax.set_xticklabels(['ER/PR+ HER2-', 'HER2+', 'Triple Negative'], ha='right', rotation=30)
ax.set_ylabel('Year of diagnosis')
ax.set_title('Diagnosis year by NAC status and biologic subtype')
ax.legend(['No NAC', 'Had NAC'], loc='lower right')

plt.tight_layout()
fig.savefig(f'{dn}/plots/year_dx_by_nac_and_biomarker.png')
plt.close()

