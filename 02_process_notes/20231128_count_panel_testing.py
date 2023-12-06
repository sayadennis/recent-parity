import os
import sys
import re
import pickle
from datetime import datetime

import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/01_ssms_raw'

with open(os.path.join(dn, 'column_names_note_type_detail.txt'), 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())

##################################
#### Genetic counseling notes #### 
##################################

## Read in SQL query results
struc = pd.read_csv(
    os.path.join(dn, 'genetic/sql_query_results_notes_genetic_counsel.csv'), 
    header=None, names=colnames
)

## Store isolated sections in this dataframe 
separated_sections = pd.DataFrame(None, columns=['ir_id', 'test_results']) # , 'interpretation' 'created_year'

patterns = [
    "Test results:",
    "Interpretation:",
    "Risk Assessment"
]

texts = pd.DataFrame('', index=struc['ir_id'].unique(), columns=patterns[:-1])

cts = 0
## Loop through patient ir_id's and isolate sections 
for i in struc['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(struc['note_text'].iloc[[x==i for x in struc['ir_id']]].values)
    patient_notelist = []
    # subsetyears = list(struc['created_year'].iloc[[x==i for x in struc['ir_id']]].values)
    for j in range(len(subsetlines)):
        if re.search("Test results:", subsetlines[j]):
            sliced_string = re.split('Test results:', subsetlines[j], 1)[1]
            sliced_string = re.split('Interpretation', sliced_string, 1)[0]
            if np.any([x in sliced_string for x in ['BRCA', 'BRCA1', 'BRCA2', 'PALB2', 'TP53', 'PTEN', 'CDH1', 'STK11', 'CHEK2', 'ATM']]):
                texts.loc[i,'Test results:'] = sliced_string
        elif re.search(r"GENETIC RESULT ADDENDUM    ADDENDUM: .{30,200} Test results ", subsetlines[j]):
            sliced_string = re.split(r"GENETIC RESULT ADDENDUM    ADDENDUM: .{30,200} Test results ", subsetlines[j], 1)[1]
            sliced_string = re.split("Interpretation", sliced_string, 1)[0]
            if len(sliced_string) > len(texts.loc[i,'Test results:']):
                texts.loc[i,'Test results:'] = sliced_string
        elif re.search("1. Genetic testing for the following genes was offered:", subsetlines[j]):
            sliced_string = re.split("1. Genetic testing for the following genes was offered:", subsetlines[j], 1)[1]
            sliced_string = re.split("2. ", sliced_string, 1)[0]
            if len(sliced_string) > len(texts.loc[i,'Test results:']):
                texts.loc[i,'Test results:'] = sliced_string
        elif re.search("1. Genetic testing for hereditary breast and ovarian cancer syndrome was offered", subsetlines[j]):
            sliced_string = re.split("1. Genetic testing for hereditary breast and ovarian cancer syndrome was offered", subsetlines[j], 1)[1]
            sliced_string = re.split("2. ", sliced_string, 1)[0]
            if ("Gene Dx for BRCA1/2 Sequencing" in sliced_string) | ("Gene Dx for BRCA1/2 NGS" in sliced_string):
                if len('Gene Dx for BRCA1, BRCA2') > len(texts.loc[i,'Test results:']):
                    texts.loc[i,'Test results:'] = 'Gene Dx for BRCA1, BRCA2'
            elif ("GeneDx for the High Risk Breast Cancer Panel" in sliced_string):
                if len('GeneDx for the High Risk Breast Cancer Panel: BRCA1, BRCA2, CDH1, PTEN, STK11 and TP53?') > len(texts.loc[i,'Test results:']):
                    texts.loc[i,'Test results:'] = 'GeneDx for the High Risk Breast Cancer Panel: BRCA1, BRCA2, CDH1, PTEN, STK11 and TP53?'
            elif "GeneDx labs for sequencing and deletion/duplication analysis of the BRCA1 and BRCA2 genes" in sliced_string:
                if len("GeneDx labs for sequencing and deletion/duplication analysis of the BRCA1 and BRCA2 genes") > len(texts.loc[i,'Test results:']):
                    texts.loc[i,'Test results:'] = "GeneDx for BRCA1, BRCA2"
            elif "Gene Dx for NGS of the BRCA genes" in sliced_string:
                texts.loc[i,'Test results:'] = "GeneDx for BRCA1, BRCA2"
        elif re.search(r"1\. Test name/genes analyzed: BRCAplus panel \(ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, PTEN, TP53\)", subsetlines[j]):
            texts.loc[i,'Test results:'] = "BRCAplus panel (ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, PTEN, TP53)"
        elif re.search(r"Invitae(?: Breast)? STAT Panel \(BRCA1, BRCA2, CDH1, PALB2, PTEN, STK11, TP53\)", subsetlines[j]):
            texts.loc[i,'Test results:'] = "Invitae Breast STAT Panel (BRCA1, BRCA2, CDH1, PALB2, PTEN, STK11, TP53)"
        elif re.search(r"1\. Test name/genes analyzed: Invitae STAT breast cancer panel \(7 genes\)", subsetlines[j]):
            texts.loc[i,'Test results:'] = "Invitae Breast STAT Panel (BRCA1, BRCA2, CDH1, PALB2, PTEN, STK11, TP53)"
        elif re.search(r"1\. Test name/genes analyzed: Invitae STAT add-on panel \(ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, PTEN, STK11, TP53\)", subsetlines[j]):
            texts.loc[i,'Test results:'] = "Invitae STAT add-on panel (ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, PTEN, STK11, TP53)"
        elif re.search(r"1\. Test name/genes analyzed: Invitae STAT add on panel \(9 genes total\)", subsetlines[j]):
            texts.loc[i,'Test results:'] = "Invitae STAT add-on panel (ATM, BRCA1, BRCA2, CDH1, CHEK2, PALB2, PTEN, STK11, TP53)"
        elif re.search(r"1. Test name/genes analyzed:", subsetlines[j]):
            sliced_string = re.split(r"1. Test name/genes analyzed:", subsetlines[j], 1)[1]
            sliced_string = re.split(r"2. ", sliced_string, 1)[0]
            texts.loc[i,'Test results:'] = sliced_string

genelist = ['BRCA', 'BRCA1', 'BRCA2', 'PALB2', 'TP53', 'PTEN', 'CDH1', 'STK11', 'CHEK2', 'ATM']

mention = pd.DataFrame(0, index=struc['ir_id'].unique(), columns=genelist)

for i in texts.index:
    for genename in genelist:
        if genename in texts.loc[i,'Test results:']:
            mention.loc[i,genename] = 1

##################################
#### Get REDCap data to align ####
##################################

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2023-03-17_0949.csv'

redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

# remove exluded patients
redcap = redcap.iloc[(redcap.exclude_demo.values!=1) & (redcap.exclude_tum.values!=1),:]

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)

# mark family history as none for patients who are missing it 
redcap.fam_hx = redcap.fam_hx.fillna(0.)

redcap = redcap.iloc[(redcap.exclude_demo.values!=1) & (redcap.exclude_tum.values!=1),:]

include_ir_ids = []
irid_to_redcapid = {}  # mapping from ir_id to redcap ID

struc.birth_date = [x.split()[0] for x in struc.birth_date]

for i in redcap.index:
    epic_mrn = redcap.loc[i,'epic_mrn']
    dob = redcap.loc[i,'dob']
    subsetlines = struc.iloc[(struc.EPIC_mrn.values==epic_mrn) & (struc.birth_date.values==dob)]
    if subsetlines.shape[0]>0:
        ir_id = subsetlines['ir_id'].iloc[0]
        redcap_id = redcap.iloc[redcap.epic_mrn.values==epic_mrn,:]['record_id'].values[0]
        include_ir_ids.append(ir_id)
        irid_to_redcapid[ir_id] = redcap_id

mention = mention.loc[include_ir_ids,:]

mention.iloc[mention.BRCA.values==1,:].sum(axis=0)/606

# Estimating what proportion of patients had panel testing: of the patients that have genes mentioned in their genetic testing results section, what proportion has mention of any other gene?  
np.any(mention.iloc[:,3:], axis=1).sum()/606


mention['Panel'] = np.any(mention[['PALB2', 'TP53', 'PTEN', 'CDH1', 'STK11', 'CHEK2', 'ATM']], axis=1).astype(int)
print('Proportion of patients with non-BRCA gene mention in results section, amonth those who have BRCA genes mentioned:', mention.iloc[np.any(mention[['BRCA', 'BRCA1', 'BRCA2']], axis=1).values,:]['Panel'].sum())

