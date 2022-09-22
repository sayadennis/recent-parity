from math import comb
import os
import sys
import numpy as np
import pandas as pd

drev = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing' # reviewed data 
dori = '/share/fsmresfiles/breast_cancer_pregnancy' # original data (from SABCS - from EDW and Erica's chart review)
dout = '/share/fsmresfiles/breast_cancer_pregnancy/data/tkhr_pat_wo_gentest_added'

###################
#### Read data ####
###################

# New data, reviewed by Takahiro 
rev = pd.read_csv(f'{drev}/04_reviewed/sayas_study_final.tsv', sep='\t')
# Corresponding processed (unreviewed) data of the above, to get basic ID info that Takahiro removed from spreadsheet
pro = pd.read_csv(f'{drev}/03_processed/recent_parity_combined_data.csv')

# SABCS data, reviewed by Erica 
ori = pd.read_csv(f'{dori}/data/recent_parity_combined_data_reviewed.csv')
ori = ori.iloc[ori['Review complete'].values=='x',:]

data_dict = pd.read_csv(f'{dori}/data/recent_parity_data_dictionary.csv')

combined = pd.DataFrame(None, columns=data_dict.ColumnName)

################################################################
#### Transform Takahiro's reviewed table to standard format ####
################################################################

## First add necessary identification information like MRNs 

rev['ir_id'] = None; rev['EPIC_mrn'] = None; rev['Powerchart_mrn'] = None; rev['birth_date'] = None;

for i in rev.index:
    first_name = rev.loc[i,'Name']
    last_name = rev.loc[i,'Name.1']
    processed_match = pro.iloc[((pro['first_name'].values==first_name) & (pro['last_name'].values==last_name)),:]
    rev.loc[i,'ir_id'] = processed_match['ir_id'].values[0]
    rev.loc[i,'EPIC_mrn'] = processed_match['EPIC_mrn'].values[0]
    rev.loc[i,'Powerchart_mrn'] = processed_match['Powerchart_mrn'].values[0]
    rev.loc[i,'birth_date'] = processed_match['birth_date'].values[0]

## 
rev.rename({
    'Name' : 'FirstName',
    'Name.1' : 'LastName',
    'birth_date' : 'DOB',
    'The age at BC diagnosis' : 'AgeBCDx',
    'The year of BC diagnosis' : 'TimeBCDx',
    'gravida' : 'Gravida',
    'para' : 'Para',
    'First menarche' : 'Menarche',
    'Age of first preg' : 'AgeFirstPreg',
    'Age of last preg' : 'AgeLastPreg',
    'Last nursed(years ago)' : 'LastNursed',
    'duration of breasting (months)' : 'DurationBF',
    'gene test (yes/no)' : 'HasGenTest',
    'genetic mutation(yes/no)' : 'Has mutation (yes/no)',
    'surgery(Tm, Bcs, Only biopcy)' : 'Surgery',
    'Tumor location [R/L/Bilateral]' : 'TumorLocation',
    'Tumor size (cm) [US or MRI]' : 'ClinTumorSize',
    'Lymphovascular Invasion' : 'LymphInv',
    'Axillary Lymph Nodes Examined' : 'NumLymphExamined',
    'Number of Positive Versus Total' : 'NumPosLymph',
    'Histologic type (DCIS,IDC,ILC or other)' : 'Histology', 
    'Tumor Grade' : 'PathGrade',
    'Pathological Tumor size(cm)' : 'PathTumorSize',
    'NAC' : 'NACT',
    'NAC regimen' : 'NACTReg',
    'response(CR,PR,SD,PD)': 'Response',
    'Adjuvant therapy (Chemo and/or hormone)' : 'Adjuvant therapy (chemo and/or hormone)',
    'Distant Free Survival (month)' : 'DRFS', 
    'OS(month)' : 'OS', 
    'Metastatic site' : 'MetaSite',
    'Others' : 'Note'
}, axis=1, inplace=True)

# Drop columns -- I think these are columns from my algorthim that Takahiro kept in the spreadsheet 
rev.drop(['ER.1', 'PR.1', 'HER2.1', 'Ki67', 'p53.1', 'Patient No.'], axis=1, inplace=True) # check code 

# Reorder columns
rev = rev[[
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'FirstName', 'LastName', 'DOB', 'AgeBCDx', 'TimeBCDx', 'Gravida', 'Para', 'Menarche', 'AgeFirstPreg', 'AgeLastPreg', 'LastNursed', 'DurationBF',
       'HasGenTest', 'Has mutation (yes/no)', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm', 'others', 'VUS', 
       'TumorLocation', 'Histology', 'ClinTumorSize', 'cTNM', 'Surgery', 'PathTumorSize', 'PathGrade', 'cTNM or pTNM', 
       'LymphInv', 'NumLymphExamined', 'NumPosLymph', 'ER', 'PR', 'HER2', 'Ki-67', 'p53', 'NACT', 'NACTReg',
       'HER2 therapy (No/Trastuzumab/More)', 'yTNM', 'Response', 'DRFS', 'OS', 'MetaSite', 'Note'
]]

## Loop through rows and write data to the "combined" dataframe 
combined = pd.concat([combined, rev])

# Enter values for pathogenic mutations 
for i in combined.index:
    # first record VUS 
    vus=combined.loc[i,'VUS'] # this will be gene name 
    for gene in ['brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']:
        if gene in vus.lower():
            combined.loc[i,gene]='VUS'
        # next record pathogenic 
        if combined.loc[i,gene]=='no':
            if combined.loc[i,gene]!='VUS':
                combined.loc[i,gene]='Negative'
        elif (('variant' in str(combined.loc[i,gene]).lower()) | ('vus' in str(combined.loc[i,gene]).lower())):
            continue
        else:
            combined.loc[i,gene]='Pathogenic'
            # print(f'Unrecognized pattern in {rev_gene.upper()} at index {i}')

# combined.drop(['brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm'], axis=1, inplace=True)

# Enter values for VUS
for i in combined.index:
    if combined.loc[i,'VUS']=='no':
        combined.loc[i,'HasVUS']=0
    else:
        combined.loc[i,'HasVUS']=1

# Clean up pathological tumor/node staging category 
combined['pTNM']=None
for i in combined.index:
    if str(combined.loc[i,'cTNM or pTNM']).startswith('p'):
        combined.loc[i,'pTNM']=combined.loc[i,'cTNM or pTNM']
    elif str(combined.loc[i,'cTNM or pTNM']).startswith('yp'):
        combined.loc[i,'pTNM']=combined.loc[i,'cTNM or pTNM']

# Clean up tumor location 
for i in combined.index:
    if 'Bilateral' in str(combined.loc[i,'TumorLocation']):
        combined.loc[i,'TumorLocation'] = 'B'
    elif combined.loc[i,'TumorLocation']=='x':
        combined.loc[i,'TumorLocation'] = np.nan

# Clean up HasGenTest column (for Andrea's review)
for i in combined.index:
    if combined.loc[i,'HasGenTest'].lower().strip()=='yes':
        combined.loc[i,'HasGenTest']=1
    elif combined.loc[i,'HasGenTest'].lower().strip()=='no':
        combined.loc[i,'HasGenTest']=0
    else:
        entry=combined.loc[i,'HasGenTest']
        print(f'Unrecognized pattern at index {i}: {entry}')

# Clean up AgeLastPreg column 
for i in combined.index:
    if combined.loc[i,'AgeLastPreg'].lower().strip()=='no':
        combined.loc[i,'HasGenTest']=np.nan

# Clean Gravida and Para
for colname in ['Gravida', 'Para']:
    for i in combined.index:
        if combined.loc[i,colname].isnumeric():
            combined.loc[i,colname]=int(combined.loc[i,colname])
        else:
            combined.loc[i,colname]=np.nan

# Clean AgeFirstPreg and AgeLastPreg
for i in combined.index:
    for colname in ['AgeFirstPreg', 'AgeLastPreg']:
        if combined.loc[i,colname].isnumeric():
            combined.loc[i,colname]=int(combined.loc[i,colname])
        elif combined.loc[i,colname]=='34(now)':
            combined.loc[i,colname]=34
        else:
            combined.loc[i,colname]=np.nan
    # Once both columns are formatted, fill in "AgeLastPreg" with first pregnancy age if empty 
    if ((~pd.isnull(combined.loc[i,'AgeFirstPreg'])) & pd.isnull(combined.loc[i,'AgeLastPreg'])):
        combined.loc[i,'AgeLastPreg'] = combined.loc[i,'AgeFirstPreg']

##############################################################################
#### Transform the data table used for SABCS into the new standard format ####
##############################################################################

ori.rename({
    'first_name' : 'FirstName', 
    'last_name' : 'LastName',
    'birth_date' : 'DOB', 
    'gravida' : 'Gravida', 
    'para' : 'Para', 
    'menarche' : 'Menarche', 
    'diagnosis_year' : 'TimeBCDx',
    'age_at_diagnosis' : 'AgeBCDx', 
    'agefirstpreg' : 'AgeFirstPreg', 
    'ageoflastpreg' : 'AgeLastPreg', 
    'lastnursed' : 'LastNursed',
    'duration' : 'DurationBF', 
    'general' : 'HasPathoMut', 
    # 'brca1' : 'BRCA1_patho', 
    # 'brca2' : 'BRCA2_patho', 
    # 'palb2' : 'PALB2_patho', 
    # 'tp53' : 'TP53_patho',
    # 'pten' : 'PTEN_patho', 
    # 'cdh1' : 'CDH1_patho', 
    # 'stk11' : 'STK11_patho', 
    # 'chek2' : 'CHEK2_patho', 
    # 'atm' : 'ATM_patho', 
}, axis=1, inplace=True)

## For each gene, clean up
for colname in ['brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']:
    for i in ori.index:
        if ori.loc[i,colname]=='positive':
            ori.loc[i,colname]='Pathogenic'
        elif ori.loc[i,colname]=='negative':
            ori.loc[i,colname]='Negative'
        elif ori.loc[i,colname]=='variant':
            ori.loc[i,colname]='VUS'
        else:
            ori.loc[i,colname]='Not performed'

ori['region_mutation']=None
for i in ori.index:
    for gene in ['brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']:
        if not pd.isnull(ori.loc[i,f'{gene}_mutregion']):
            if pd.isnull(ori.loc[i,'region_mutation']):
                ori.loc[i,'region_mutation']=gene.upper()+'-'+ori.loc[i,f'{gene}_mutregion']+'; '
            else:
                ori.loc[i,'region_mutation']=ori.loc[i,'region_mutation']+gene.upper()+'-'+ori.loc[i,f'{gene}_mutregion']+'; '

## Change DOB format 


###################################################
#### Concatenate and data and clean duplicates ####
###################################################

combined = combined.iloc[rev['HasGenTest'].values=='yes',]

combined = pd.concat([ori, combined])
# combined = combined[data_dict.ColumnName]

# clean duplicates 
combined.drop_duplicates(subset='ir_id', keep='last', inplace=True, ignore_index=True) # inplace=True

dup_ix = combined.iloc[[sum(combined['EPIC_mrn']==x)>1 for x in combined['EPIC_mrn']],:].index # weird ones where EPIC MRN are the same but names and all info are different... 
combined.loc[dup_ix,'EPIC_mrn']=None

# save combined table, with only necessary columns (?)
combined.to_csv(f'{dout}/recent_parity_combined_data_reviewed_tkhr_added.csv', index=False, header=True)
