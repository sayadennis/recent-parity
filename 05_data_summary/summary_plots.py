import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

#####################################
#### Load data & data dictionary ####
#####################################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-10_0938.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

plotdir = 'plots'

nat_only = True

if nat_only:
    data = data.iloc[data.nat.values==1,]
    plotdir = 'plots/nac_only'

#####################
#### Demographis ####
#####################

## Age at diagnosis 
fig, ax = plt.subplots()
ax.hist(data.age_at_diagnosis)
ax.set_xlabel('Age')
ax.set_ylabel('Number of patients')
ax.set_title('Distribution of Cancer Diagnosis Age')
plt.savefig(f'{dn}/{plotdir}/age_at_diagnosis.png')

###############################
#### Gynecological history ####
###############################


###################################
#### For Andrea's presentation ####
###################################

# Andrea's table: tumor characteristics separated by parity information 
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

## Scatter plot: clinical vs pathological tumor size 
fig, ax = plt.subplots(figsize=(8,8))
for par_cat in data.parity_category.unique():
    par_data = data.iloc[data.parity_category.values==par_cat,:]
    ax.scatter(par_data.clin_tumor_size, par_data.tumor_size, label=par_cat, alpha=0.5) # .map(colors)
ax.set_xlabel('Clinical size (cm)')
ax.set_ylabel('Pathological size (cm)')
ax.legend()
fig.savefig(f'{dn}/{plotdir}/clin_path_size_scatter.png')
plt.close()

## Heatmap: clinical vs pathological node staging category
node = pd.crosstab(
    data.clin_node_stag_cat.map(dd['clin_node_stag_cat']['Choices, Calculations, OR Slider Labels']), 
    data.node_staging_category.map(dd['node_staging_category']['Choices, Calculations, OR Slider Labels']), 
)
node.loc[
    ['NX','N0','N1','N2','N3'],
    ['ypNX','ypN0','ypN1','ypN1mi','ypN2','ypN3']
]

fig, ax = plt.subplots()
im = ax.imshow(node)

ax.set_xticks(np.arange(node.shape[1]))
ax.set_xticklabels(node.columns)
ax.set_yticks(np.arange(node.shape[0]))
ax.set_yticklabels(node.index)

ax.set_ylabel('Clinical node category')
ax.set_xlabel('Pathological node category')

for i in range(node.shape[0]):
    for j in range(node.shape[1]):
        text = ax.text(j, i, node.iloc[i, j],
                       ha="center", va="center", color="w")

fig.tight_layout()
fig.savefig(f'{dn}/{plotdir}/clin_path_node_heatmap.png')
plt.close()

