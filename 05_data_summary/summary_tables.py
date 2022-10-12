import numpy as np
import pandas as pd
import pickle

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

## Load data & data dictionary
datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-10_0938.csv'
data = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# Select only patients with neoadjuvant therapy if necessary
sumtabdir = 'summary_tables'

nat_only = True

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

data['biomarker_subtypes'] = None
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        er = False
        pr = False
        print('Missing ER/PR at index:', i)
    her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    if ((er | pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'ER/PR+ HER2-'
    elif ((not er) & (not pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'Triple Negative'
    elif her2:
        data.loc[i,'biomarker_subtypes'] = 'HER2+'
    else:
        print('Unknown pattern at index:', i)

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


if nat_only:
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

    data['er_pr'] = None
    for i in data.index:
        try:
            er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
            pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
        except:
            er = False
            pr = False
            print('Missing ER/PR at index:', i)
        if (er | pr):
            data.loc[i,'er_pr'] = 'ER/PR+'
        else:
            data.loc[i,'er_pr'] = 'ER/PR-'

    #### Age ####
    show_tab = pd.DataFrame(
        index=['Mean age (±SD)'],
        columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    )
    for par_cat in show_tab.columns:
        mean = data.iloc[data.parity_category.values==par_cat,:].age_at_diagnosis.mean()
        std = data.iloc[data.parity_category.values==par_cat,:].age_at_diagnosis.std()
        show_tab.loc['Mean age (±SD)',par_cat] = f'{mean:.2f} (±{std:.2f})'
    # 
    show_tab.to_csv(f'{dn}/{sumtabdir}/demo_age_vs_parity.csv')

    #### Clinical T staging category ####
    show_tab = pd.crosstab(data.clin_tumor_stag_cat, data.parity_category).rename(
        dd['clin_tumor_stag_cat']['Choices, Calculations, OR Slider Labels'],
        axis=0
    )
    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ].to_csv(f'{dn}/{sumtabdir}/counts_clintumor_vs_parity.csv')

    #### Clinical N staging category ####
    show_tab = pd.crosstab(data.clin_node_stag_cat, data.parity_category).rename(
        dd['clin_node_stag_cat']['Choices, Calculations, OR Slider Labels'],
        axis=0
    )
    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    ].to_csv(f'{dn}/{sumtabdir}/counts_clinnode_vs_parity.csv')

    #### Biomarker subtypes ####
    show_tab = pd.crosstab(data.biomarker_subtypes, data.parity_category)
    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ].to_csv(f'{dn}/{sumtabdir}/counts_biomarkersubtype_vs_parity.csv')

    #### Treatment regimen ####
    show_tab = pd.crosstab(
        data.nat_reg.map(dd['nat_reg']['Choices, Calculations, OR Slider Labels']), 
        data.parity_category
    )
    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ].to_csv(f'{dn}/{sumtabdir}/counts_natregimen_vs_parity.csv')

    #### RCB ####
    show_tab = pd.crosstab(
        data.rcb_category.map(dd['rcb_category']['Choices, Calculations, OR Slider Labels']), 
        data.parity_category
    )
    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ].to_csv(f'{dn}/{sumtabdir}/counts_rcbcategory_vs_parity.csv')

    #### RCB and ER/PR status ####
    show_tab = pd.DataFrame(data.groupby(['er_pr', 'rcb_category', 'parity_category']).size(), columns=['count'])
    show_tab = show_tab.reset_index(level='parity_category').pivot(columns='parity_category')
    # clean up column names 
    levels = show_tab.columns.levels
    labels = show_tab.columns.labels
    show_tab.columns = levels[1][labels[1]]
    show_tab.columns.name = None
    show_tab.to_csv(f'{dn}/{sumtabdir}/counts_rcbcategory_vs_erpr_vs_parity.csv')

    #### Genetic status - do all genes ####
    show_tab = pd.DataFrame(columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years'])
    for gene in ['brca1', 'brca2', 'palb2', 'tp53', 'chek2', 'pten', 'cdh1', 'stk11', 'atm']:
        try:
            show_tab = pd.concat((
                show_tab, 
                pd.DataFrame(pd.crosstab(
                    data[gene].map(dd[gene]['Choices, Calculations, OR Slider Labels']), 
                    data.parity_category
                ).loc['Pathogenic',:].rename(gene)).transpose()
            ))
        except:
            print('No mutations in ', gene)
            show_tab.loc[gene,:] = 0

    show_tab.loc[
        :,['Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ].to_csv(f'{dn}/{sumtabdir}/counts_genetic_vs_parity.csv')
