import os
import sys
import re
import numpy as np
import pandas as pd
from datetime import datetime

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

# with open(f'{dn}/column_names_note_type_detail.txt', 'r') as f:
#     lines = f.readlines()

# colnames = []
# for line in lines:
#     colnames.append(line.strip())

colnames = [ # the alternative colnames for when focusing on note types etc. 
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date',
    'department_key', 'department_name', 'department_external_name', 'department_specialty', 'location_name',
    'created_year', 'updated_year', 'signed_year', 'source_system_table', 'note_type', 'note_status', 'note_text'
]

#########################
#### Pathology notes ####
#########################

## Read in SQL query results
data = pd.read_csv(
    f'{dn}/pathology/notes_pathology.csv', 
    header=None, names=colnames
)

data['template'] = None

for i in range(data.shape[0]):
    if len(re.findall(r'ANTIGEN RESULT/INTENSITY OF STAINING INTERPRETATION', re.sub(r'\s+', ' ', data.iloc[i,-2])))>0:
        data.iloc[i,-1] = 'ANTIGEN RESULT'
    elif 'Breast Cancer Staging Summary' in data.iloc[i,-2]:
        data.iloc[i,-1] = 'Breast Cancer Staging Summary'
    elif 'BREAST CANCER STAGING SUMMARY' in data.iloc[i,-2]:
        data.iloc[i,-1] = 'BREAST CANCER STAGING SUMMARY'
    else:
        data.iloc[i,-1] = None


#######################################
#### Breast Cancer Staging Summary ####
#######################################

bcs = data.iloc[data['template'].values=='Breast Cancer Staging Summary',:]

def find_and_slice(string, begin, end):
    """
    string : string to be sliced
    begin : REGEX indicating where the slice begins
    end : REGEX expression indicating where the slice ends
    """
    if len(re.findall(begin, string))>0: # if a match exists for beginning of slice
        if len(re.findall(end, string))>0: # if a match exists for end of slice
            if re.split(begin, string)[1]: # if the string after beginning pattern is not None ?? - is this step still necessary after fixing the parentheses problem? 
                string_slice = re.split(end, re.split(begin, string)[1])[0].strip()
            else: # if string after beginning pattern is None
                string_slice = re.split(end, re.split(begin, string)[2])[0].strip()
        else:
            string_slice = re.split(begin, string)[1].strip()
    else: # if no match exists for beginning of slice 
        string_slice = None
    return string_slice


def get_tumor_char(notestring):
    # notestring = bcs['note_text'].iloc[i]
    # first eliminate the part preceding the staging summary 
    staging_sum = re.sub(r'\s+', ' ', notestring.split('Breast Cancer Staging Summary')[1])
    # now collect the entries
    tum_char_dict = {
        'Specimen Submitted'                        : find_and_slice(staging_sum, begin='Specimen Submitted:', end=r'Specimen Dimensions:|Tumor Size:'),
        'Specimen Dimensions'                       : find_and_slice(staging_sum, begin='Specimen Dimensions:', end=r'Tumor [Ss]ize:|Tumor Size \(greatest dimension of each\):'),
        'Tumor Size'                                : find_and_slice(staging_sum, begin=r'Tumor [Ss]ize:|Tumor Size \(greatest dimension of each\):', end=r'Histologic Type:|Tumor Type:|Final Diagnosis|Lymphovascular Invasion:'),
        'Histologic Type'                           : find_and_slice(staging_sum, begin=r'Histologic Type:|Tumor Type:', end=r'Grade:|Grade \(cumulative\):'),
        'Grade (cumulative)'                        : find_and_slice(staging_sum, begin=r'Grade:|Grade \(cumulative\):', end=r'Tubule Formation:|Lymph[o\-]vascular Invasion:|Lymphatic Vascular Invasion:'),
        'Tubule Formation'                          : find_and_slice(staging_sum, begin='Tubule Formation:', end=r'Nuclear Pleomorphism:|Nuclear Atypia:'),
        'Nuclear Pleomorphism/Atypia'               : find_and_slice(staging_sum, begin=r'Nuclear Pleomorphism:|Nuclear Atypia:', end='Mitotic Rate:'),
        'Mitotic Rate'                              : find_and_slice(staging_sum, begin='Mitotic Rate:', end=r'Lymph[o\-]vascular Invasion:|Lymph-vascular Space Invasion:|Grade \(cumulative\):'),
        'Lymphovascular Invasion'                   : find_and_slice(staging_sum, begin=r'Lymph[o\-]vascular Invasion:|Lymphatic Vascular Invasion:', end=r'DCIS as Extensive Intraductal Component:|DCIS as Extensive Component:|LCIS:'),
        'DCIS as Extensive Intraductal Component'   : find_and_slice(staging_sum, begin='DCIS as Extensive Intraductal Component:', end='DCIS Measurement/Proportion:'),
        'DCIS Measurement/Proportion'               : find_and_slice(staging_sum, begin='DCIS Measurement/Proportion:', end=r'LCIS:|Lobular neoplasia Focal ALH Calcifications:'),
        'LCIS'                                      : find_and_slice(staging_sum, begin=r'LCIS:', end=r'Calcifications:'),
        'Calcifications'                            : find_and_slice(staging_sum, begin=r'Calcifications:|Lobular neoplasia Focal ALH Calcifications:', end=r'Locations of|Margins of Excision:'), # its 'Locations of Calcification:', but the "Calcification" part will be cut off with the 'begin' slice
        'Locations of Calcifications'               : find_and_slice(staging_sum, begin='Locations of Calcifications:', end=r'Margins of Excision:?|# Axillary Lymph Nodes Examined:'),
        'Margins of Excision'                       : find_and_slice(staging_sum, begin=r'Margins of Excision:|Negative Distance to Margin:', end=r'Invasive Cancer:|Invasive Cancer and DCIS:|DCIS:'),
        'Invasive Cancer'                           : find_and_slice(staging_sum, begin='Invasive Cancer:', end=r'Distance to Margin:|Axillary Lymph Nodes Number of Positive Versus Total:'),
        'Distance to Margin'                        : find_and_slice(staging_sum, begin='Distance to Margin:', end=r'Margin Widely Free \(more than 0.5 cm\)|DCIS:'),
        'Margin Widely Free (more than 0.5 cm)'     : find_and_slice(staging_sum, begin=r'Margin Widely Free \(more than 0.5 cm\)|Negative Distance to Margin:', end=r'DCIS:|Axillary Lymph Nodes Number of Positive Versus Total:'),
        'DCIS'                                      : find_and_slice(staging_sum, begin='DCIS:', end=r'Distance to Margin:|# Axillary Lymph Nodes Examined:'), # why is there another distance to margin here? ,
        '# Axillary Lymph Nodes Examined'           : find_and_slice(staging_sum, begin='# Axillary Lymph Nodes Examined:', end=r'Number of Positive Versus Total:|Breast Tumor Markers:'),
        'Number of Positive Versus Total'           : find_and_slice(staging_sum, begin=r'Number of Positive Versus Total:|Axillary Lymph Nodes Number of Positive Versus Total:', end=r'Size of Largest Metastasis:|Breast Tumor Markers:'),
        'Size of Largest Metastasis'                : find_and_slice(staging_sum, begin='Size of Largest Metastasis:', end='Extranodal Extension:'),
        'Extranodal Extension'                      : find_and_slice(staging_sum, begin='Extranodal Extension:', end='Breast Tumor Markers:'),
        'Breast Tumor Markers'                      : find_and_slice(staging_sum, begin='Breast Tumor Markers:', end=r'ER:|TNM Staging:'),
        'ER'                                        : find_and_slice(staging_sum, begin='ER:', end='PR:'),
        'PR'                                        : find_and_slice(staging_sum, begin='PR:', end='HER2:'),
        'HER2'                                      : find_and_slice(staging_sum, begin='HER2:', end=r'Ki-67|Diagnosis'),
        'Ki-67'                                     : find_and_slice(staging_sum, begin='Ki-67', end=r'p53:?|Labs:|TNM Staging:'),
        'p53'                                       : find_and_slice(staging_sum, begin=r'p53:?', end=r'TNM Staging:|IMPRESSION AND RECOMMENDATION:|Final Diagnosis|ASSESSMENT/PLAN:|IMPRESSION:|RS of|Breast Biopsy Case Number:'),
        'TNM'                                       : find_and_slice(staging_sum, begin='TNM Staging:', end=r'Grading of invasive carcinoma is based on|Labs:|Final Diagnosis|BREAST CANCER TUMOR MARKERS FOR INVASIVE CARCINOMA:|Past Medical History:|[1-9][1-9]?/[1-9][1-9]?/[0-9][0-9]|Oncotype|PMH/PSH:|SUBJECTIVE:|ASSESSMENT:|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|MX PE:'),
        'Oncotype RS'                               : find_and_slice(staging_sum, begin=r'Oncotype RS|Oncotype testing', end='Impression and Plan'),
        'Diagnosis'                                 : find_and_slice(staging_sum, begin=r'Joint Commission on Cancer \(AJCC, 7th Edition, 2010\). Diagnosis|Final Diagnosis', end=r'Impression and Plan|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|IMPRESSION:'),
        'Impression and Plan'                       : find_and_slice(staging_sum, begin=r'Impression and Plan|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|IMPRESSION:', end='Seen and discussed with')
    }
    return tum_char_dict

tum_char = pd.DataFrame()

for i in range(bcs.shape[0]):
    tum_char_dict = get_tumor_char(bcs['note_text'].iloc[i])
    for varname in ['ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date']:
        tum_char_dict[varname] = bcs[varname].iloc[i]
    tum_char = tum_char.append(tum_char_dict, ignore_index=True)


## Save a subset of columns temporarily for exploratory purposes
columns_to_save = [
    'ir_id',
    'EPIC_mrn',
    'Powerchart_mrn',
    'first_name',
    'last_name',
    'birth_date',
    'Specimen Submitted',
    'Specimen Dimensions',
    'Tumor Size',
    'Histologic Type',
    'Grade (cumulative)',
    'Tubule Formation',
    'Nuclear Pleomorphism/Atypia',
    'Mitotic Rate',
    'Lymphovascular Invasion',
    'DCIS as Extensive Intraductal Component',
    'DCIS Measurement/Proportion',
    'LCIS',
    'Calcifications',
    'Locations of Calcifications',
    'Margins of Excision',
    'Invasive Cancer',
    'Distance to Margin',
    'Margin Widely Free (more than 0.5 cm)',
    'DCIS',
    '# Axillary Lymph Nodes Examined',
    'Number of Positive Versus Total',
    'Size of Largest Metastasis',
    'Extranodal Extension',
    'Breast Tumor Markers',
    'ER',
    'PR',
    'HER2',
    'Ki-67',
    'p53',
    'TNM',
]
tum_char[columns_to_save].drop_duplicates(ignore_index=True).to_csv(f'{dn}/pathology/isolated_sections_pathology.csv', index=False, header=True)

# varname = 'Diagnosis'
# for i in range(500):
#     print(tum_char[varname].iloc[i])

########################
#### ANTIGEN RESULT ####
########################

ant = data.iloc[data['template'].values=='ANTIGEN RESULT',:]

def get_tumor_char(notestring):
    # notestring = bcs['note_text'].iloc[i]
    # first eliminate the part preceding the staging summary 
    notestring = re.sub(r'\s+', ' ', notestring).split('ANTIGEN RESULT/INTENSITY OF STAINING INTERPRETATION')[1]
    # now collect the entries
    tum_char_dict = {
        'ER'                                        : find_and_slice(notestring, begin='Estrogen Receptor', end='Progesterone Receptor'),
        'PR'                                        : find_and_slice(notestring, begin='Progesterone Receptor', end='HER-2/neu Score'),
        'HER2'                                      : find_and_slice(notestring, begin='HER-2/neu Score', end=r'Ki-67'),
        'Ki-67'                                     : find_and_slice(notestring, begin='Ki-67', end=r'p53'),
        'p53'                                       : find_and_slice(notestring, begin=r'p53:?', end=r'Family History:|1\.|Impression/Plan:|[1-9][1-9]?/[1-9][1-9]?/[0-9][0-9]|Staging Scans:|PMH/PSH:|TNM Staging:|Invasive Breast Cancer Staging Summary|Assessment and Plan:|RS|Breast Biopsy Case Number:|REFERENCE RANGES:'),
        'Genetic Testing'                           : find_and_slice(notestring, begin=r'Genetic testing:', end=r'Ashkenazi Jewish heritage:|Time since diagnosis:|TREATMENT SIDE EFFECTS:')
    }
    return tum_char_dict

tum_char = pd.DataFrame()

for i in range(ant.shape[0]):
    tum_char_dict = get_tumor_char(ant['note_text'].iloc[i])
    for varname in ['ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date']:
        tum_char_dict[varname] = ant[varname].iloc[i]
    tum_char = tum_char.append(tum_char_dict, ignore_index=True)

tum_char

