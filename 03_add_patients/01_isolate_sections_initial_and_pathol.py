import os
import sys
import re
import numpy as np
import pandas as pd
from datetime import datetime

din = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/01_raw'
dout = '/share/fsmresfiles/breast_cancer_pregnancy/patients_wo_genetic_testing/02_interim'

###############################
#### Initial counsel notes #### 
###############################

colnames = [ # the alternative colnames for when focusing on note types etc. 
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date',
    'created_year', 'updated_year', 'signed_year', 'source_system_table', 'note_type', 'note_status', 'note_text'
]

## Read in SQL query results
struc = pd.read_csv(
    f'{din}/initialcounsel_data_patients_wo_genetic_testing.csv',
    header=None, names=colnames
)

## Store isolated sections in this dataframe 
separated_sections = pd.DataFrame(None, columns=['ir_id', 'title_type', 'created_year', 'gyn_history'])

## Loop through patient ir_id's and isolate sections 
for i in struc['ir_id'].unique():
    # select only notes that are for the given patient 
    subsetlines = list(struc['note_text'].iloc[[x==i for x in struc['ir_id']]].values)
    subsetyears = list(struc['created_year'].iloc[[x==i for x in struc['ir_id']]].values)
    sectionnotfound = []
    for j in range(len(subsetlines)):
        note = subsetlines[j]
        if 'Initial Consultation' in note:
            if 'GYNECOLOGICAL HISTORY:' in note:
                begin = note.split('GYNECOLOGICAL HISTORY:')[1]
                gyn_his = begin.split('REVIEW OF SYSTEMS:')[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Initial Consultation', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        elif 'Consultation' in note:
            # find the section of interest
            if 'Breast cancer risk factors:' in note:
                begin = re.split('Breast cancer risk factors:', note, 1)[1]
                gyn_his = re.split('FAMILY HISTORY|Family history:|Family History of Breast or Ovarian Cancer:|SOCIAL HISTORY:|Family Hx', begin)[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Consultation 1', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            elif 'Breast Cancer Risk Factors/Gynecologic History:' in note:
                begin = re.split('Breast Cancer Risk Factors/Gynecologic History:', note, 1)[1]
                gyn_his = begin.split('Family Cancer History:')[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'Consultation 2', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        elif 'New Patient History and Physical' in note:
            if 'BREAST RISK DATA' in note:
                begin = note.split('BREAST RISK DATA')[1]
                gyn_his = re.split('Other Risk exposure:|Other R isk exposure|SOCIAL HISTORY:|Social History|Social Hx:', begin)[0].strip()
                separated_sections = separated_sections.append(
                    {'ir_id' : i, 'title_type' : 'New Patient History and Physical', 'created_year' : int(subsetyears[j]), 'gyn_history' : gyn_his}, ignore_index=True
                )
                sectionnotfound.append(0)
            else:
                sectionnotfound.append(1)
        else:
            continue
    if np.all(sectionnotfound==0):
        print('Section of interest not found in patient %s' % i)


separated_sections.drop_duplicates(inplace=True, ignore_index=True)
separated_sections.to_csv(f'{dout}/isolated_gynecological_history_section_initial_counsel.csv', sep='!', header=True, index=False)

# init = separated_sections.iloc[[x=='Initial Consultation' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# con1 = separated_sections.iloc[[x=='Consultation 1' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# con2 = separated_sections.iloc[[x=='Consultation 2' for x in separated_sections['title_type'].values],:].reset_index(drop=True)
# newp = separated_sections.iloc[[x=='New Patient History and Physical' for x in separated_sections['title_type'].values],:].reset_index(drop=True)

#########################
#### Pathology notes ####
#########################

colnames = [ # the alternative colnames for when focusing on note types etc. 
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date',
    'department_key', 'department_name', 'department_external_name', 'department_specialty', 'location_name',
    'created_year', 'updated_year', 'signed_year', 'source_system_table', 'note_type', 'note_status', 'note_text'
]

## Read in SQL query results
data = pd.read_csv(
    f'{din}/pathol_data_patients_wo_genetic_testing.csv', 
    header=None, names=colnames
)

data['template'] = None

for i in range(data.shape[0]):
    if 'Breast Cancer Staging Summary' in data.iloc[i,-2]:
        data.iloc[i,-1] = 'Breast Cancer Staging Summary'
    elif len(re.findall(r'ANTIGEN RESULT/INTENSITY OF STAINING INTERPRETATION', re.sub(r'\s+', ' ', data.iloc[i,-2])))>0:
        data.iloc[i,-1] = 'ANTIGEN RESULT'
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
        'Specimen Dimensions'                       : find_and_slice(staging_sum, begin='Specimen Dimensions:', end=r'Tumor [Ss]ize:|Tumor Size \(greatest dimension of each\):|Lymphovascular Invasion:'),
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
        'Calcifications'                            : find_and_slice(staging_sum, begin=r'Calcifications:|Lobular neoplasia Focal ALH Calcifications:', end=r'Locations of|Margins of Excision:|Page'), # its 'Locations of Calcification:', but the "Calcification" part will be cut off with the 'begin' slice
        'Locations of Calcifications'               : find_and_slice(staging_sum, begin='Locations of Calcifications:', end=r'Margins of Excision:?|# Axillary Lymph Nodes Examined:'),
        'Margins of Excision'                       : find_and_slice(staging_sum, begin=r'Margins of Excision:|Negative Distance to Margin:', end=r'Invasive Cancer:|Invasive Cancer and DCIS:|DCIS:'),
        'Invasive Cancer'                           : find_and_slice(staging_sum, begin='Invasive Cancer:', end=r'Distance to Margin:|Axillary Lymph Nodes Number of Positive Versus Total:|DCIS'),
        'Distance to Margin'                        : find_and_slice(staging_sum, begin='Distance to Margin:', end=r'Margin Widely Free \(more than 0.5 cm\)|DCIS:|Axillary Lymph Nodes Number of Positive Versus Total:|# Axillary Lymph Nodes Examined:|Axillary Lymph Nodes:'),
        'Margin Widely Free (more than 0.5 cm)'     : find_and_slice(staging_sum, begin=r'Margin Widely Free \(more than 0.5 cm\)|Negative Distance to Margin:', end=r'DCIS:|Axillary Lymph Nodes|# Axillary Lymph Nodes Examined:'),
        'DCIS'                                      : find_and_slice(staging_sum, begin=r'DCIS:', end=r'Distance to Margin:|# Axillary Lymph Nodes Examined:'), # why is there another distance to margin here? ,
        '# Axillary Lymph Nodes Examined'           : find_and_slice(staging_sum, begin='# Axillary Lymph Nodes Examined:', end=r'Number of Positive Versus Total:|Breast Tumor Markers:'),
        'Number of Positive Versus Total'           : find_and_slice(staging_sum, begin=r'Number of Positive Versus Total:|Axillary Lymph Nodes Number of Positive Versus Total:|Axillary Lymph Nodes Number of Positive Versus Total:|Number of Positive:', end=r'Size of Largest Metastasis:|Breast Tumor Markers:|micromet|Breast Tumor Markers'),
        'Size of Largest Metastasis'                : find_and_slice(staging_sum, begin='Size of Largest Metastasis:', end='Extranodal Extension:'),
        'Extranodal Extension'                      : find_and_slice(staging_sum, begin='Extranodal Extension:', end=r'Breast Tumor Markers:|BREAST CANCER TUMOR MARKERS|TNM Staging:'),
        'Breast Tumor Markers'                      : find_and_slice(staging_sum, begin='Breast Tumor Markers:', end=r'ER:|TNM Staging:'),
        'ER'                                        : find_and_slice(staging_sum, begin='ER:', end=r'PR:|Alcohol use:'),
        'PR'                                        : find_and_slice(staging_sum, begin='PR:', end='HER2:'),
        'HER2'                                      : find_and_slice(staging_sum, begin='HER2:', end=r'Ki-67|Diagnosis'),
        'Ki-67'                                     : find_and_slice(staging_sum, begin='Ki-67', end=r'\.|p53:?|Labs:|TNM Staging:|CURRENT THERAPY:|CT|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|Assessment and Plan:'),
        'p53'                                       : find_and_slice(staging_sum, begin=r'p53:?', end=r'TNM Staging:|IMPRESSION AND RECOMMENDATION:|Final Diagnosis|ASSESSMENT/PLAN:|IMPRESSION:|RS of|Breast Biopsy Case Number:|Note:|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|Tumor Bank:|Impression and Plan|Genetic Testing:|Case Number:|Assessment and Plan|Mammogram|CT'),
        'TNM'                                       : find_and_slice(staging_sum, begin='TNM Staging:', end=r',|Grading of invasive carcinoma is based on|Labs:|Final Diagnosis|BREAST CANCER TUMOR MARKERS FOR INVASIVE CARCINOMA:|Past Medical History:|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|Oncotype|PMH/PSH:|SUBJECTIVE:|ASSESSMENT:|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|MX PE:|A/P:|Procedures'),
        'Oncotype RS'                               : find_and_slice(staging_sum, begin=r'Oncotype RS|Oncotype testing', end='Impression and Plan'),
        'Diagnosis'                                 : find_and_slice(staging_sum, begin=r'Joint Commission on Cancer \(AJCC, 7th Edition, 2010\). Diagnosis|Final Diagnosis', end=r'Impression and Plan|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|IMPRESSION:'),
        'Impression and Plan'                       : find_and_slice(staging_sum, begin=r'Impression and Plan|IMPRESSION AND RECOMMENDATION:|ASSESSMENT/PLAN:|IMPRESSION:', end='Seen and discussed with'),
        'Genetic Testing'                           : find_and_slice(notestring, begin=r'Genetic testing:', end=r'Ashkenazi Jewish heritage:|Time since diagnosis:|TREATMENT SIDE EFFECTS:|Labs:')
    }
    return tum_char_dict

tum_char = pd.DataFrame()

for i in range(bcs.shape[0]):
    tum_char_dict = get_tumor_char(bcs['note_text'].iloc[i])
    for varname in ['ir_id', 'first_name', 'last_name', 'birth_date']:
        tum_char_dict[varname] = bcs[varname].iloc[i]
    tum_char = tum_char.append(tum_char_dict, ignore_index=True)


## Save a subset of columns temporarily for exploratory purposes
columns_to_save = [
    'ir_id',
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
    'Genetic Testing'
]

tum_char = tum_char[columns_to_save].drop_duplicates(ignore_index=True)
tum_char_bcs = tum_char
bcs_irid = list(tum_char_bcs['ir_id'].unique())
# tum_char.to_csv(f'{dout}/isolated_sections_pathology.csv', index=False, header=True)

########################
#### ANTIGEN RESULT ####
########################

ant = data.iloc[data['template'].values=='ANTIGEN RESULT',:]
ant = ant.iloc[[x not in bcs_irid for x in ant['ir_id']],:]

def get_tumor_char(notestring):
    # notestring = bcs['note_text'].iloc[i]
    # first eliminate the part preceding the staging summary 
    notestring = re.sub(r'\s+', ' ', notestring).split('ANTIGEN RESULT/INTENSITY OF STAINING INTERPRETATION')[1]
    # now collect the entries
    tum_char_dict = {
        'ER'                                        : find_and_slice(notestring, begin='Estrogen Receptor', end=r'Progesterone Receptor'),
        'PR'                                        : find_and_slice(notestring, begin='Progesterone Receptor', end=r'HER-2/neu [Ss]core|HER-2/neu|C-erB-2'),
        'HER2'                                      : find_and_slice(notestring, begin=r'HER-2/neu [Ss]core|HER-2/neu|C-erB-2', end=r'Ki-67|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|CT|Assessment|Brain mri|Oncotype|Guardant:|REFERENCE RANGES:|A/P:|Results|Labs:|[0-9]/[0-9][0-9]|ASSESSMENT/PLAN:|Positive Diagnosis|Final Diagnosis'),
        'Ki-67'                                     : find_and_slice(notestring, begin='Ki-67', end=r'\.|p53|CT|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|Assessment and Plan:'),
        'p53'                                       : find_and_slice(notestring, begin=r'p53:?', end=r'Family History:|1\.|Impression/Plan:|[0-9][0-9]?/[0-9][0-9]?/[0-9][0-9]|Staging Scans|PMH/PSH:|TNM Staging:|Invasive Breast Cancer Staging Summary|Assessment and Plan|RS|Breast Biopsy Case Number:|REFERENCE RANGES:|Oncotype:|RECOMMENDATIONS/FOLLOW-UP:|HER2|Final Diagnosis|Diagnosis|Note for Specimen|BREAST CANCER TUMOR MARKERS|Impression:|\.|Part|ASSESSMENT/PLAN:|Assessment/Plan:|A/P:|Staging|Impression and Plan|Addendum Report|Note:|Labs:|Imaging|Mammogram|CT|Received|This report|CONCLUSIONS|Case|Genetic Testing:')
    }
    return tum_char_dict

tum_char = pd.DataFrame()

for i in range(ant.shape[0]):
    tum_char_dict = get_tumor_char(ant['note_text'].iloc[i])
    for varname in ['ir_id', 'first_name', 'last_name', 'birth_date']:
        tum_char_dict[varname] = ant[varname].iloc[i]
    tum_char = tum_char.append(tum_char_dict, ignore_index=True)

columns_to_save = [
    'ir_id',
    'first_name',
    'last_name',
    'birth_date',
    'ER',
    'PR',
    'HER2',
    'Ki-67',
    'p53'
]

tum_char = tum_char[columns_to_save].drop_duplicates(ignore_index=True)
tum_char = tum_char_bcs.append(tum_char, ignore_index=True)
tum_char.to_csv(f'{dout}/isolated_sections_pathology.csv', index=False, header=True)
