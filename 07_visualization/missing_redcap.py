import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-01-17_1220.csv'
redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

tumchar = np.array(list(dd.keys()))[[(dd[x]['Form Name']=='tumor_characteristics') & (dd[x]['Required Field?']=='y') for x in dd.keys()]]
nat = np.array(list(dd.keys()))[[(dd[x]['Form Name']=='treatment') & (dd[x]['Required Field?']=='y') for x in dd.keys()]]

fig, ax = plt.subplots()
ax.bar(np.arange(len(tumchar)), [pd.isnull(redcap[x]).sum() for x in tumchar])
ax.set_xticks(np.arange(len(tumchar)))
ax.set_xticklabels([dd[x]['Field Label'] for x in tumchar], rotation=45, ha='right')
ax.set_ylabel('Number of patients')
# ax.set_xlabel('Data Field')
ax.set_title('Number of Patients Currently Missing Tumor Characteristics Data')

# ax[1].bar(np.arange(len(nat)), [pd.isnull(redcap[x]).sum() for x in nat])
# ax[1].set_xticks(np.arange(len(nat)))
# ax[1].set_xticklabels([dd[x]['Field Label'] for x in nat], rotation=45, ha='right')
# ax[1].set_xlabel('Data Field')
# ax[1].set_title('Neoadjuvant Therapy')

# fig.suptitle('Number of Patients Currently Missing Data')
plt.tight_layout()

fig.savefig(f'{dn}/plots/redcap_missing_tumchar_nat_summary.png')
plt.close()
