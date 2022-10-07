import numpy as np
import pandas as pd
import pickle

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

## Load data & data dictionary

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-07_1648.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# Select only patients with neoadjuvant therapy if necessary
sumtabdir = 'summary_tables'

nat_only = False

if nat_only:
    data = data.iloc[data.nat.values==1,]
    sumtabdir = 'summary_tables/nac_only'

##########################################
#### Demographic & Gynecological info ####
##########################################

# Parity distribution - 4 categories 
pd.DataFrame(
    [np.sum(data.parous.values==0),
     np.sum(data.years_since_pregnancy<5),
     np.sum((data.years_since_pregnancy>=5) & (data.years_since_pregnancy<10)),
     np.sum(data.years_since_pregnancy>=10)],
    index=['nulliparous', 'less than 5', 'less than 10', '10 or more'],
    columns=['Number of patients']
)

# Demographic & Gynecological history table 
gyn_demo_table = pd.DataFrame(columns=['Mean (±SD)'])
ignore = ['parous', 'breastfed', 'breastfeeding_duration', 'recency_of_lactation']

for key in dd.keys():
    if key=='age_at_diagnosis':
        mean = np.nanmean(data[key])
        std = np.nanstd(data[key])
        gyn_demo_table = pd.concat((gyn_demo_table, pd.DataFrame(
            f'{mean:.2f} (±{std:.2f})',
            index=[dd[key]['Field Label']],
            columns=['Mean (±SD)']
        )))
    elif dd[key]['Form Name']=='parity_information':
        if key not in ignore:
            mean = np.nanmean(data[key])
            std = np.nanstd(data[key])
            gyn_demo_table = pd.concat((gyn_demo_table, pd.DataFrame(
                f'{mean:.2f} (±{std:.2f})',
                index=[dd[key]['Field Label']],
                columns=['Mean (±SD)']
            )))

gyn_demo_table.to_csv(f'{dn}/{sumtabdir}/gyn_demo.csv')

###################
#### Treatment ####
###################

# Presence of NAC
data.nat.map({1: 'Yes', 0: 'No'}).value_counts().sort_values()

# NAC regimen distribution
df = pd.DataFrame(
    data.nat_reg.map(dd['nat_reg']['Choices, Calculations, OR Slider Labels']).value_counts().values,
    index=data.nat_reg.map(dd['nat_reg']['Choices, Calculations, OR Slider Labels']).value_counts().index,
    columns=['Number of patients']
).loc[['Anthracycline/FEC/EC/AC', 'Taxane/DOC/PAC/TC', 'Both', 'Other']]
print(df)
df.to_csv(f'{dn}/{sumtabdir}/nac_regimen.csv')

for key in dd.keys():
    if dd[key]['Form Name']=='treatment':
        if isinstance(dd[key]['Choices, Calculations, OR Slider Labels'], dict):
            cts = pd.DataFrame(
                data[key].map(dd[key]['Choices, Calculations, OR Slider Labels']).value_counts().values,
                index=data[key].map(dd[key]['Choices, Calculations, OR Slider Labels']).value_counts().index,
                columns=['counts']
            )
            if cts.shape[0]>1:
                print(f'''#### {dd[key]['Field Label']} ####''')
                print(cts)
                print()


##############################
#### Genetic test results ####
##############################

# Common mutations
genetic = pd.DataFrame()

gyn_category_bool = {
    'Nulliparous'       : data.parous.values==0,
    'Less than 5 years' : data.years_since_pregnancy.values<5,
    '5-10 years'        : ((data.years_since_pregnancy.values>=5) & (data.years_since_pregnancy.values<10)),
    '10 or more years'  : data.years_since_pregnancy.values>=10
}

for key in dd.keys():
    if dd[key]['Form Name']=='geneticsfam_hx':
        if dd[key]['Field Type']=='radio':
            if key in ['any_patho_mutation', 'any_vus_mutation']:
                cts = {}
                cts['Total'] = np.sum(data[key].map(dd[key]['Choices, Calculations, OR Slider Labels'])=='Present')
                for gyn_cat in gyn_category_bool.keys():
                    cts[gyn_cat] = np.sum(data.iloc[gyn_category_bool[gyn_cat]][key].map(
                        dd[key]['Choices, Calculations, OR Slider Labels']
                    )=='Present')
            else:
                cts = {}
                cts['Total'] = np.sum(data[key].map(dd[key]['Choices, Calculations, OR Slider Labels'])=='Pathogenic')
                for gyn_cat in gyn_category_bool.keys():
                    cts[gyn_cat] = np.sum(data.iloc[gyn_category_bool[gyn_cat]][key].map(
                        dd[key]['Choices, Calculations, OR Slider Labels']
                    )=='Pathogenic')
            genetic = pd.concat((genetic, pd.DataFrame(cts, index=[key.upper()])))

# Re-order columns 
genetic = genetic[['Nulliparous', 'Less than 5 years', '5-10 years', '10 or more years', 'Total']]
genetic.to_csv(f'{dn}/{sumtabdir}/genetic_mutation_counts.csv')

###############################
#### Tumor characteristics ####
###############################

# Biomarker subtypes
status = {}
for biomarker in ['er', 'pr', 'her2']:
    status[biomarker] = (data[f'{biomarker}_status'].map(
        dd[f'{biomarker}_status']['Choices, Calculations, OR Slider Labels']
    )=='Positive').values

df = pd.DataFrame(
    [np.sum(status['er'] & status['pr'] & ~status['her2']),
     np.sum(~status['er'] & ~status['pr'] & ~status['her2']),
     np.sum(status['her2'])],
    index=['ER/PR+ HER2-', 'Triple Negative', 'HER2+'],
    columns=['counts']
)
print(df)
df.to_csv(f'{dn}/{sumtabdir}/biomarker_subtypes.csv')

for key in dd.keys():
    if dd[key]['Form Name']=='tumor_characteristics':
        if ((dd[key]['Field Type']=='radio') | (dd[key]['Field Type']=='dropdown')):
            cts = data[key].map(dd[key]['Choices, Calculations, OR Slider Labels']).value_counts()
            if len(cts)>1:
                print(f'''#### {dd[key]['Field Label']} ####''')
                df = pd.DataFrame(
                    cts.sort_index().values,
                    index=cts.sort_index().index,
                    columns=['counts']
                )
                print(df)
                print('')
                df.to_csv(f'{dn}/{sumtabdir}/{key}.csv')
