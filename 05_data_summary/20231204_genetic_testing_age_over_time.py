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

########################################
#### Get some counts for evaluation ####
########################################

## Diagnosis years -- Theresa let me know that patients diagnosed in 2015 or later is more likely to have received panel testing
num_dx = {i : len(ages[i]) for i in range(2010,2021)}

post_2015 = sum([num_dx[i] for i in range(2015,2021)])
pre_2015 = sum([num_dx[i] for i in range(2010,2015)])

print(f'Of our whole REDCap cohort, {100*post_2015/(pre_2015+post_2015):.1f}% of our cohort was diagnosed/tested in 2015 or later.')

## When we count
exec(open('/home/srd6051/recent_parity/recent-parity/02_process_notes/20231128_count_panel_testing.py').read())

irid_wo_genenames = mention.iloc[np.all(mention.drop('Panel', axis=1, inplace=False)==0, axis=1).values,:].index
no_mention_genename = struc.iloc[[x in irid_wo_genenames for x in struc.ir_id.values],:]

likely_bart = no_mention_genename.iloc[['BRACAnalysis' in x for x in no_mention_genename.note_text],:].ir_id.unique()
nobart = no_mention_genename.iloc[['BRACAnalysis' not in x for x in no_mention_genename.note_text],:].ir_id.unique()
nobart_notes = struc.iloc[[x in nobart for x in struc.ir_id],:].note_text

print(f'{len(likely_bart)} ({100 * (len(likely_bart))/mention.shape[0]:.1f}%) patients had mention of BRACAnalysis (BART?) in their genetic counseling notes.')

possibly_no_panel = len(likely_bart) + len(nobart)
print(f'Among patients whose genetic counsel notes are available, at most {possibly_no_panel} ({100 * possibly_no_panel/mention.shape[0]:.1f}%) patients had non-panel testing.')

mention.index = [irid_to_redcapid[x] for x in mention.index]  # change `mention` index from MRN to redcap ID 
mined = redcap.iloc[[x in mention.index for x in redcap.record_id.values],:]  # get the rows of redcap data that were mined from EDW (NOT manually entered later by Takahiro)
post_2015 = sum(mined.year_of_diagnosis>=2015)
pre_2015 = sum(mined.year_of_diagnosis<2015)

print(f'Among the {mined.shape[0]} patients whose data were obtained via text mining, {100*post_2015/(pre_2015+post_2015):.1f}% of our cohort was diagnosed/tested in 2015 or later.')

