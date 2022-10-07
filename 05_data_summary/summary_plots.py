import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

#####################################
#### Load data & data dictionary ####
#####################################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-09-26_1703.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

plotdir = 'plots'

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


