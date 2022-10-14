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

nat_only = False

if nat_only:
    data = data.iloc[data.nat.values==1,]
    sumtabdir = 'summary_tables/nac_only'

#########################################
#### Fill in useful categorical info ####
#########################################

## Fill in parity category 
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

## Fill in biomarker subtypes 
data['biomarker_subtypes'] = None
missing_ix = []
for i in data.index:
    try:
        er = dd['er_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'er_status']]=='Positive'
        pr = dd['pr_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'pr_status']]=='Positive'
    except:
        er = False
        pr = False
        if not i in missing_ix:
            missing_ix.append(i)
    if not pd.isnull(data.loc[i,'her2_status']):
        her2 = dd['her2_status']['Choices, Calculations, OR Slider Labels'][data.loc[i,'her2_status']]=='Positive'
    else:
        her2 = False
        if not i in missing_ix:
            missing_ix.append(i)
    if ((er | pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'ER/PR+ HER2-'
    elif ((not er) & (not pr) & (not her2)):
        data.loc[i,'biomarker_subtypes'] = 'Triple Negative'
    elif her2:
        data.loc[i,'biomarker_subtypes'] = 'HER2+'
    else:
        print('Unknown pattern at index:', i)

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

super_table = pd.DataFrame(
    columns=['super_category', 'main_category', 'sub_category', 
             'Nulliparous', '<5 years', '5-10 years', '>=10 years',
             'Total']
)

#####################
#### Demographic ####
#####################

#### Age ####
main_table = pd.DataFrame(
    index=['Mean age (±SD)'],
    columns=super_table.columns
)
for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    mean = data.iloc[data.parity_category.values==par_cat,:].age_at_diagnosis.mean()
    std = data.iloc[data.parity_category.values==par_cat,:].age_at_diagnosis.std()
    main_table.loc['Mean age (±SD)',par_cat] = f'{mean:.2f} (±{std:.2f})'

main_table['super_category'] = 'demo'
main_table['main_category'] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((
    super_table,
    main_table
))

###############################
#### Gynecological history ####
###############################

main_table = pd.DataFrame(
    index=['Mean Gravida (±SD)', 'Mean Para (±SD)'],
    columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years']
)

for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    # Gravida 
    mean = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.mean()
    std = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.std()
    main_table.loc['Mean Gravida (±SD)',par_cat] = f'{mean:.2f} (±{std:.2f})'
    # Para 
    mean = data.iloc[data.parity_category.values==par_cat,:].number_births.mean()
    std = data.iloc[data.parity_category.values==par_cat,:].number_births.std()
    main_table.loc['Mean Para (±SD)',par_cat] = f'{mean:.2f} (±{std:.2f})'

main_table['super_category'] = 'gynecological'
main_table['main_category'] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((
    super_table,
    main_table
))

##############################
#### Genetic test results ####
##############################

# Common mutations
main_table = pd.DataFrame(
    columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years']
)

for key in dd.keys():
    if dd[key]['Form Name']=='geneticsfam_hx':
        if dd[key]['Field Type']=='radio':
            if key in ['any_patho_mutation', 'any_vus_mutation']:
                cts = {}
                for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
                    cts[par_cat] = np.sum(data.iloc[data.parity_category.values==par_cat,:][key].map(
                        dd[key]['Choices, Calculations, OR Slider Labels']
                    )=='Present')
            else:
                cts = {}
                for par_cat in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
                    cts[par_cat] = np.sum(data.iloc[data.parity_category.values==par_cat,:][key].map(
                        dd[key]['Choices, Calculations, OR Slider Labels']
                    )=='Pathogenic')
            main_table = pd.concat((main_table, pd.DataFrame(cts, index=[key.upper()])))

main_table['Total'] = main_table[
    ['Nulliparous', '<5 years', '5-10 years', '>=10 years']
].sum(axis=1).astype(int)

main_table['super_category'] = 'genetic'
main_table['main_category'] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((
    super_table,
    main_table
))

###############################
#### Tumor characteristics ####
###############################

main_table = pd.DataFrame(
    columns=['main_category', 'sub_category', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
)

for colname in ['biomarker_subtypes', 'clin_tumor_stag_cat', 'clin_node_stag_cat']:    
    if colname!='biomarker_subtypes':
        subtable = pd.crosstab(data[colname], data['parity_category']).rename(
            dd[colname]['Choices, Calculations, OR Slider Labels'],
            axis=0
        ).reset_index()
    else:
        subtable = pd.crosstab(data[colname], data['parity_category']).reset_index()
    #
    subtable = subtable.rename({colname : 'sub_category'}, axis=1)
    subtable.columns.name = None
    subtable['main_category'] = colname
    main_table = pd.concat((
        main_table,
        subtable
    ))

main_table['Total'] = main_table[
    ['Nulliparous', '<5 years', '5-10 years', '>=10 years']
].sum(axis=1).astype(int)

main_table['super_category'] = 'tumor_char'

super_table = pd.concat((
    super_table,
    main_table
))

###################
#### Treatment ####
###################

main_table = pd.DataFrame(
    columns=['main_category', 'sub_category', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
)

for colname in ['nat', 'nat_reg', 'her2_therapy', 'rcb', 'rcb_category']:
    if colname=='nat':
        sub_table = pd.crosstab(data[colname], data['parity_category']).rename(
            {1: 'Yes', 0: 'No'}, 
            axis=0
        )
    elif colname=='rcb':
        continue
        ## The below needs to NOT be part of this loop. Fix ## 
        # sub_table = pd.DataFrame(
        #     index=['Mean RCB (±SD)'], 
        #     columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years']
        # )
        # for par_cat in sub_table.columns:
        #     mean = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.mean()
        #     std = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.std()
        # sub_table.loc['Mean RCB (±SD)',par_cat] = f'{mean} (±{std})'
    else:
        sub_table = pd.crosstab(data[colname], data['parity_category']).rename(
            dd[colname]['Choices, Calculations, OR Slider Labels'],
            axis=0
        )
    sub_table = sub_table.reset_index()
    sub_table.columns.name = None
    sub_table = sub_table.rename({colname : 'sub_category'}, axis=1)
    sub_table['main_category'] = colname
    main_table = pd.concat((main_table, sub_table))


main_table['Total'] = main_table[
    ['Nulliparous', '<5 years', '5-10 years', '>=10 years']
].sum(axis=1).astype(int)

main_table['super_category'] = 'treatment'
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((
    super_table,
    main_table
))

super_table.to_csv(f'{dn}/{sumtabdir}/super_summary_table.csv', index=False)
