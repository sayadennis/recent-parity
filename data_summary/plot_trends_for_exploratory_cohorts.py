import os
import numpy as np
import pandas as pd
import re
from datetime import datetime
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

din='/share/fsmresfiles/breast_cancer_pregnancy/explore_cohort'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

## Load data 
colnames=['ir_id', 'Epic_MRN', 'Powerchart_MRN', 'first_name', 'last_name', 'dob', 'earliest_diagnosis']

bc = pd.read_csv(f'{din}/cohort_by_breastcanceronly_withdiagyear.csv', header=None, names=colnames)
prov = pd.read_csv(f'{din}/cohort_by_prov_name_withdiagyear.csv', header=None, names=colnames)
genn = pd.read_csv(f'{din}/cohort_by_genetic_notes_withdiagyear.csv', header=None, names=colnames)

## Clean data
for df in [bc, prov, genn]:
    # add year of diagnosis
    df['year_of_dx'] = [datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f').year for x in df['earliest_diagnosis'].values]
    # add age of diagnosis
    df['age_of_dx'] = [(datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f') - datetime.strptime(y, '%Y-%m-%d %H:%M:%S.%f')).days // 365 for (x,y) in zip(df['earliest_diagnosis'].values, df['dob'].values)]
    # reformat dob to exclude time (only year, month, and day)
    df['dob'] = [x.split()[0] for x in df['dob'].values]
    # eliminate unnecessary columns and remove duplicates
    df.drop(['Epic_MRN', 'Powerchart_MRN', 'earliest_diagnosis'], axis=1, inplace=True)
    df.drop_duplicates(ignore_index=True, inplace=True)

## Plot Venn diagram of three groups 
olsize = 0
for ir_id in prov['ir_id'].unique():
    if ir_id in genn['ir_id'].values:
        olsize += 1

venn2(subsets=(len(prov['ir_id'].unique())-olsize, len(genn['ir_id'].unique())-olsize, olsize), set_labels=['Provider encounter', 'Genetic notes'])
plt.title('Overlap of cohort identified by two methods')
plt.savefig(f'{dout}/prov_genn_cohort_overlap_venn.png')
plt.close()

## Plot histograms of diagnosis years for the three groups 
bins_list = np.arange(2010, 2022)
fig, axs = plt.subplots(1,3, figsize=(8,4))
fig.suptitle('Histograms of diagnosis year', fontsize=15)
axs[0].hist(bc['year_of_dx'], bins=bins_list)
axs[0].set_title('All breast cancer')
axs[1].hist(prov['year_of_dx'], bins=bins_list)
axs[1].set_title('Provider encounter')
axs[2].hist(genn['year_of_dx'], bins=bins_list)
axs[2].set_title('Genetic notes')
fig.savefig(f'{dout}/histogram_year_of_dx.png')
plt.close()

## Plot histograms of diagnosis ages for the three groups
fig, axs = plt.subplots(1,3, figsize=(8,4))
fig.suptitle('Histograms of age of diagnosis', fontsize=15)
axs[0].hist(bc['age_of_dx'])
axs[0].set_title('All breast cancer')
axs[1].hist(prov['age_of_dx'])
axs[1].set_title('Provider encounter')
axs[2].hist(genn['age_of_dx'])
axs[2].set_title('Genetic notes')
fig.savefig(f'{dout}/histogram_age_of_dx.png')
plt.close()

