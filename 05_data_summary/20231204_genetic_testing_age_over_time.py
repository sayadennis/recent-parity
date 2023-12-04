from datetime import datetime
import pickle

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

## Load data & data dictionary
datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# remove exluded patients
data = data.iloc[(data.exclude_demo.values!=1) & (data.exclude_tum.values!=1),:]

##################
## Total counts ##
##################

data = data.iloc[data.year_of_diagnosis.values>=2010,:]

# Record ages as array
ages = {}
for i in range(2010,2021):
    ages[i] = data.iloc[data.year_of_diagnosis.values==i,:].age_at_diagnosis.values

thres = 40  # age threshold

categories = [f'>{thres} years old', f'<={thres} years old']
age_category_ct = {
    i : [sum(ages[i]>thres), sum(ages[i]<=thres)] for i in range(2010,2021)
}

age_category_pct = {
    i : [100 * sum(ages[i]>thres)/len(ages[i]), 100 * sum(ages[i]<=thres)/len(ages[i])] for i in range(2010,2021)
}

fig, ax = plt.subplots(1, 2, figsize=(10,5))

ax[0].stackplot(age_category_ct.keys(), np.array(list(age_category_ct.values())).T, labels=categories, alpha=0.7)
ax[0].set_xlabel('Dx year')
ax[0].set_ylabel('Number of patients')
#ax[0].set_title('Number of patients')
ax[0].legend(loc='upper left')

ax[1].stackplot(age_category_pct.keys(), np.array(list(age_category_pct.values())).T, labels=categories, alpha=0.7)
ax[1].set_xlabel('Dx year')
ax[1].set_ylabel('Percentage')
#ax[1].set_title(f'Proportion of patients') 
ax[1].legend(loc='upper left')
ax[1].set_ylim(0,100)

fig.suptitle(f'Temporal changes in age distributions of patients with genetic testing')

plt.tight_layout()
fig.savefig(f'{dn}/plots/age_genetic_testing_over_time_{thres}yo.png')
plt.close()

