import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

############################################
#### Identifying staging summary tables ####
############################################

## select rows where path note includes a mastectomy or lumpectomy
# m_pt = m_pt.iloc[[(('astectomy' in x) | ('umpectomy' in x)) for x in m_pt['pathology note text']],:]

table_patterns = [
    ## generic "breast cancer" 
    '\nBREAST CANCER STAGING SUMMARY\n',
    '\nLEFT BREAST CANCER STAGING SUMMARY\n',
    '\nRIGHT BREAST CANCER STAGING SUMMARY\n',
    '\nBreast Carcinoma Checklist:',
    ## DCIS 
    '\nDUCTAL CARCINOMA IN SITU SUMMARY', # '\n' , in case it looks like 'DUCTAL CARCINOMA IN SITU SUMMARY (Combined with 01-S-11-33962)'
    '\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n',
    '\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary',
    '\n Ductal Carcinoma In Situ Staging Summary', # ' \n' ':' ','
    '\nIn Situ Breast Carcinoma Checklist:\n', # not necessarily DCIS but could be LCIS? 
    '\nIn-Situ Breast Carcinoma Checklist:\n',
    ## Invasive 
    '\n Invasive Breast Cancer Staging Summary', # ' ', '-', ':', ','
    '\nInvasive Breast Carcinoma Checklist', # '\n', ':'
    '\nINVASIVE BREAST CANCER SUMMARY', # '\n' , in case it looks like 'INVASIVE BREAST CANCER SUMMARY (LEFT BREAST)' 
    '\nINVASIVE RIGHT BREAST CANCER SUMMARY\n',
    '\nINVASIVE LEFT BREAST CANCER SUMMARY\n',
    '\n Microinvasive Breast Cancer Staging Summary',
    ## postneoadjuvant 
    '\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary', # '  \n'
    '\n POST NEOADJUUVANT Invasive Breast Cancer Staging Summary',
    '\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary',
    '\n Post-treatment Breast Cancer Staging Summary',
    '\n Post-Treatment Invasive Breast Cancer Staging Summary',
    '\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ', # ' \n'
    '\n Invasive Breast Cancer POST-TREATMENT Staging Summary',
    '\n Invasive Breast Cancer Post-Treatment Staging Summary',
    '\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary',
    # '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary - RIGHT BREAST',
    '\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ',
    '\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n',
]

##################################################
#### Record tumor characteristics from tables ####
##################################################

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
        string_slice = ''
    # clean up by removing large white spaces to single space and converting to lower case 
    string_slice = re.sub(r'\s+', ' ', string_slice).lower()
    return string_slice


def get_histology(note_text, table_title):
    if table_title=='\nBREAST CANCER STAGING SUMMARY\n':
        # get sections for invasive and in-situ histologies 
        invasive_section = find_and_slice(note_text, begin=r'INFILTRATING TUMOR TYPE:', end=r'GRADE:')
        insitu_section = find_and_slice(note_text, begin=r'DUCTAL CARCINOMA IN-SITU TYPE:', end=r'NUCLEAR GRADE:')
        # determine histology 
        if 'ductal' in invasive_section:
            hist='IDC'
        elif 'lobular' in invasive_section:
            hist='ILC'
        elif 'ductal' in insitu_section:
            hist='DCIS'
        elif 'lobular' in insitu_section:
            hist='LCIS'
        else:
            hist=None
        return hist
    elif (
        (table_title=='\nINVASIVE BREAST CANCER SUMMARY') | 
        (table_title=='\nINVASIVE LEFT BREAST CANCER SUMMARY\n') | 
        (table_title=='\nINVASIVE RIGHT BREAST CANCER SUMMARY\n')
    ):
        section = find_and_slice(note_text, begin=r'Histologic Type:', end=r'Grade:')
        if 'duct' in section:
            hist='IDC'
        elif 'lobular' in section:
            hist='ILC'
        elif 'mucinous' in section:
            hist='mucinous'
        elif 'tubular' in section:
            hist='tubular'
        else:
            hist='other'
        return hist
    elif (
        (table_title=='\nDUCTAL CARCINOMA IN SITU SUMMARY') |
        (table_title=='\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n') | 
        (table_title=='\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary') | 
        (table_title=='\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ Staging Summary') |
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary')
    ):
        return 'DCIS'
    elif (
        (table_title=='\n Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Post-Treatment Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer POST-TREATMENT Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer Post-Treatment Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ') | 
        # (table_title=='\n Microinvasive Breast Cancer Staging Summary') |
        (table_title=='\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n')
    ): # invasive by nature 
        section = find_and_slice(note_text, begin=r'Histologic Type:', end=r'Grade')
        if 'duct' in section:
            hist='IDC'
        elif 'lobular' in section:
            hist='ILC'
        elif 'mucinous' in section:
            hist='mucinous'
        elif 'tubular' in section:
            hist='tubular'
        else:
            hist='other'
        return hist
    elif (
        # (table_title=='\n Post-treatment Breast Cancer Staging Summary') | 
        (table_title=='\nIn Situ Breast Carcinoma Checklist:\n') |
        (table_title=='\nIn-Situ Breast Carcinoma Checklist:\n')
    ): # in-situ by nature 
        if 'DCIS' in note_text:
            hist='DCIS'
        elif 'LCIS' in note_text:
            hist='LCIS'
        elif 'duct' in note_text.lower():
            hist='DCIS'
        elif 'lobular' in note_text.lower():
            hist='LCIS'
        else:
            hist=None
        return hist
    elif (
        # (table_title=='\nBreast Carcinoma Checklist:') |
        (table_title=='\nInvasive Breast Carcinoma Checklist')
    ):
        invasive_section=find_and_slice(note_text, begin='Invasive carcinoma:\n Type:', end='Expected molecular subtype:')
        if 'duct' in invasive_section:
            hist='DCIS'
        elif 'lobular' in invasive_section:
            hist='LCIS'
        else:
            hist=None
        return hist


def get_grade(note_text, table_title):
    if table_title=='\nBREAST CANCER STAGING SUMMARY\n':
        section = find_and_slice(note_text, begin=r'GRADE:', end=r'SIZE:')
        section = re.sub(r'\s+', ' ', section)
        grades = re.findall(r'[0-9] ', section)
        grades = [int(x) for x in grades]
        if len(grades)>0:
            grades = [int(x) for x in grades]
            return np.max(grades)
        else:
            return None
    elif table_title=='\n Invasive Breast Cancer Staging Summary':
        section = find_and_slice(note_text, begin=r'Grade \(cumulative\):|Grade:', end=r'\(|of|Tubule Formation|Lymphovascular Invasion|Lymph-vascular Invasion:|Lymphatic Vascular Invasion:|Architectural Patterns:|TNM Staging:|LB:kl|Microinvasive carcinoma:')
        if section is not None:
            grades = re.findall(r' [0-9] ', f' {section} ') # pad with spaces to recognize grades in cases like '3 and 2' 
            if len(grades)>0:
                grades = [int(x) for x in grades]
                return np.max(grades)
        else:
            return None
    elif (
        (table_title=='\nINVASIVE BREAST CANCER SUMMARY') | 
        (table_title=='\nINVASIVE LEFT BREAST CANCER SUMMARY\n') | 
        (table_title=='\nINVASIVE RIGHT BREAST CANCER SUMMARY\n')
    ):
        section = find_and_slice(note_text, begin=r'Grade:', end=r'Lymphatic Vascular Invasion:|Lymphatic/vascular Invasion:')
        if section is not None:
            section = re.sub(r'III', '3', section); section = re.sub(r'II', '2', section); section = re.sub(r'I', '1', section)
            grades = re.findall(r' [0-9] ', f' {section} ') # pad with spaces to recognize grades in cases like '3 and 2' 
            if len(grades)>0:
                grades = [int(x) for x in grades]
                return np.max(grades)
        else:
            return None
    elif table_title=='\nDUCTAL CARCINOMA IN SITU SUMMARY':
        return None # DCIS doesn't have histologic grade 
    elif (
        (table_title=='\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n') | 
        (table_title=='\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary') | 
        (table_title=='\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ Staging Summary') |
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary')
    ):
        return None # DCIS doesn't have histologic grade 
    elif (
        (table_title=='\n Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Post-Treatment Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer POST-TREATMENT Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer Post-Treatment Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n') | 
        (table_title=='\n Microinvasive Breast Cancer Staging Summary')
    ):
        section = find_and_slice(note_text, begin=r'Grade (cumulative):|Grade:', end=r'Tubule Formation:|Architectural Patterns:')
        if section is not None:
            grades = re.findall(r' [0-9] ', f' {section} ') # pad with spaces to recognize grades in cases like '3 and 2' 
            if len(grades)>0:
                grades = [int(x) for x in grades]
                return np.max(grades)
        else:
            return None
    elif (
        (table_title=='\n Post-treatment Breast Cancer Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary') |
        (table_title=='\nIn Situ Breast Carcinoma Checklist:\n') |
        (table_title=='\nIn-Situ Breast Carcinoma Checklist:\n')
    ):
        return None # doesn't seem to have grade 
    elif (
        (table_title=='\nBreast Carcinoma Checklist:') |
        (table_title=='\nInvasive Breast Carcinoma Checklist')
    ):
        section = find_and_slice(note_text, begin=r'Histologic grade (scale of I to III):', end=r'Nottingham histologic score')
        if section is not None:
            section = re.sub(r'III', '3', section); section = re.sub(r'II', '2', section); section = re.sub(r'I', '1', section)
            grades = re.findall(r' [0-9] ', f' {section} ') # pad with spaces to recognize grades in cases like '3 and 2' 
            if len(grades)>0:
                grades = [int(x) for x in grades]
                return np.max(grades)
        else:
            return None


def get_size(note_text, table_title):
    if table_title=='\nBREAST CANCER STAGING SUMMARY\n':
        section = find_and_slice(note_text, begin=r'SIZE:\.*>?', end=r'DUCTAL CARCINOMA IN-SITU TYPE|IN-SITU/ INTRADUCTAL TUMOR TYPE|IN GREATEST DIMENSION')
        section = re.sub(r'\s+', ' ', section)
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\n Invasive Breast Cancer Staging Summary':
        section = find_and_slice(note_text, begin=r'Tumor [Ss]ize:|Tumor Size \(greatest dimension of each\):', end=r'Histologic Type:|\(Based on the|\(see|Tumor Type:|Final Diagnosis|Lymphovascular Invasion:|TNM Staging')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\nINVASIVE BREAST CANCER SUMMARY':
        section = find_and_slice(note_text, begin=r'Tumor Size :|Tumor [Ss]ize:|Tumor Size \([a-zA-Z ]{0,100}\):', end=r'Histologic Type')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\nDUCTAL CARCINOMA IN SITU SUMMARY':
        return None # this type of table doesn't seem to have size 
    elif (
        (table_title=='\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n') | 
        (table_title=='\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary') | 
        (table_title=='\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ Staging Summary')
    ):
        section = find_and_slice(note_text, begin=r'Overall estimate of the size of DCIS:|Overall estimate of the size:|Present on one slide \(largest diameter\):|Size/Extent of Lesion: Present on one slide \n Largest diameter:|Largest single diameter:', end=r'Calcifications')
        if section is not None:
            sizes = re.findall(r'[0-9]\.[0-9]', section)
            if len(sizes)>0:
                return np.max([float(x) for x in sizes])
            else:
                return None
    elif (
        (table_title=='\n Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Post-Treatment Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer POST-TREATMENT Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer Post-Treatment Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n')
    ):
        section = find_and_slice(note_text, begin=r'Tumor Size:|Overall estimate of the size of DCIS:', end=r'Histologic Type:|Calcifications:')
        if section is not None:
            if 'No residual tumor identified' in section:
                return 0.
            else:
                sizes = re.findall(r'[0-9]\.[0-9]', section)
                if len(sizes)>0:
                    return np.max([float(x) for x in sizes])
                else:
                    return None
    elif table_title=='\nINVASIVE BREAST CANCER SUMMARY':
        section = find_and_slice(note_text, begin=r'Tumor Size', end=r'Histologic Type:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\nINVASIVE LEFT BREAST CANCER SUMMARY\n':
        section = find_and_slice(note_text, begin=r'Tumor Size:', end=r'Histologic Type:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\nINVASIVE RIGHT BREAST CANCER SUMMARY\n':
        section = find_and_slice(note_text, begin=r'Tumor Size:', end=r'Histologic Type:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\n Microinvasive Breast Cancer Staging Summary':
        section = find_and_slice(note_text, begin=r'Tumor Size:', end=r'Histologic Type:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif table_title=='\n Post-treatment Breast Cancer Staging Summary':
        section = find_and_slice(note_text, begin=r'Tumor Size:', end=r'Lymphovascular Invasion:')
        if 'No residual invasive carcinoma' in section:
            return 0.
        else:
            sizes = re.findall(r'[0-9]\.[0-9]', section)
            if len(sizes)>0:
                return np.max([float(x) for x in sizes])
            else:
                return None
    elif table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary':
        section = find_and_slice(note_text, begin=r'Tumor Size:', end=r'Histologic Type:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None
    elif (
        (table_title=='\nBreast Carcinoma Checklist:') |
        (table_title=='\nIn Situ Breast Carcinoma Checklist:\n') |
        (table_title=='\nIn-Situ Breast Carcinoma Checklist:\n') | 
        (table_title=='\nInvasive Breast Carcinoma Checklist')
    ):
        section = find_and_slice(note_text, begin=r'Size:|Size (microscopic):', end=r'Grossly evident lesion:|Location within breast:')
        sizes = re.findall(r'[0-9]\.[0-9]', section)
        if len(sizes)>0:
            return np.max([float(x) for x in sizes])
        else:
            return None


def get_ki67(note_text, table_title):
    if table_title=='\nBREAST CANCER STAGING SUMMARY\n':
        return None # this table type doesn't have Ki-67
    elif table_title=='\n Invasive Breast Cancer Staging Summary':
        section = find_and_slice(note_text, begin=r'Ki-67', end=r'TNM Staging:|p53')
        pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
        if len(pc)>0: # if percentage available
            ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
        else: # if perfcentage not available 
            if (('see below' in section) | ('pending' in section) | ('see table below' in section)):
                if 'BREAST CANCER TUMOR MARKERS' in note_text:
                    markers=note_text.split('BREAST CANCER TUMOR MARKERS')[1]
                    section = find_and_slice(markers, begin=r'Ki-67', end=r'p53')
                    pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
                    if len(pc)>0:
                        ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
                    else: # no ki-67 in the addendum 
                        ki67=None
                else: # no tumor marker addendum 
                    ki67=None
            else: # if there's not "see below" etc. - find "low" "intermediate" or "high"
                if 'low' in section:
                    ki67='low'
                elif 'intermediate' in section:
                    ki67='intermediate'
                elif 'high' in section:
                    ki67='high'
                else:
                    ki67=None
        return ki67
    elif (
        (table_title=='\nINVASIVE BREAST CANCER SUMMARY') | 
        (table_title=='\nINVASIVE LEFT BREAST CANCER SUMMARY\n') | 
        (table_title=='\nINVASIVE RIGHT BREAST CANCER SUMMARY\n')
    ):
        section = find_and_slice(note_text, begin=r'Ki-67:', end=r'HER-2/neu:')
        pc = re.findall(r'[0-9]{1,3}%', section)
        if len(pc)>0: # if percentage available
            ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
        else: # if perfcentage not available 
            if 'BREAST CANCER TUMOR MARKERS' in note_text:
                markers=note_text.split('BREAST CANCER TUMOR MARKERS')[1]
                section = find_and_slice(markers, begin=r'Ki-67', end=r'p53')
                pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
                if len(pc)>0:
                    ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
                else: # no ki-67 in the addendum 
                    ki67=None
            else: # no tumor marker addendum 
                if 'low' in section:
                    ki67='low'
                elif 'intermediate' in section:
                    ki67='intermediate'
                elif 'high' in section:
                    ki67='high'
                else:
                    ki67=None
        return ki67
    elif table_title=='\nDUCTAL CARCINOMA IN SITU SUMMARY':
        return None # this type doesn't have Ki-67
    elif (
        (table_title=='\n Ductal Carcinoma In Situ with Microinvasion Staging Summary  \n') | 
        (table_title=='\n Ductal Carcinoma In Situ with DCIS and Microinvasion Staging Summary') | 
        (table_title=='\n Post-Neoadjuvant Ductal Carcinoma In Situ Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ Staging Summary') |
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary')
    ):
        return None # DCIS doesn't seems to have Ki-67
    elif (
        (table_title=='\n Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n POST NEO-ADJUUVANT Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Post-Treatment Invasive Breast Cancer Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT THERAPY Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer POST-TREATMENT Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer Post-Treatment Staging Summary') | 
        (table_title=='\n Invasive Breast Cancer POST-NEOADJUVANT TREATMENT Staging Summary ') | 
        (table_title=='\n Invasive Breast Cancer Staging Summary - POST-TREATMENT  \n') | 
        (table_title=='\n Microinvasive Breast Cancer Staging Summary')
    ):
        section = find_and_slice(note_text, begin=r'Ki-67', end=r'p53:|TNM Staging:')
        pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
        if len(pc)>0: # if percentage available
            ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
        else: # if perfcentage not available 
            addendum=re.findall(r'BREAST CANCER POST-TREATMENT TUMOR MARKERS|BREAST CANCER TUMOR MARKERS|BREAST CANCER TUMOR POST-TREATMENT MARKERS|BREAST CANCER POST-TREATMEN TUMOR MARKERS', note_text)
            if len(addendum)>0:
                markers=note_text.split(addendum[0])[1]
                section = find_and_slice(markers, begin=r'Ki-67', end=r'p53')
                pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
                if len(pc)>0:
                    ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
                else: # no ki-67 in the addendum 
                    ki67=None
            else:
                ki67=None
        return ki67
    elif (
        (table_title=='\n Post-treatment Breast Cancer Staging Summary') | 
        (table_title=='\n Ductal Carcinoma In Situ POST-TREATMENT Staging Summary') |
        (table_title=='\nIn Situ Breast Carcinoma Checklist:\n') |
        (table_title=='\nIn-Situ Breast Carcinoma Checklist:\n')
    ):
        return None # doesn't seem to have ki-67
    elif (
        # (table_title=='\nBreast Carcinoma Checklist:') | # this one doesn't have Ki-67
        (table_title=='\nInvasive Breast Carcinoma Checklist')
    ):
        section = find_and_slice(note_text, begin=r'Ki-67 Manual Quantitative Immunohistochemistry', end=r'C-erbB2')
        pc = re.findall(r'[0-9]{1,3}%', section) # record percentage point if available
        if len(pc)>0: # if percentage available
            ki67 = np.max([int(x[:-1]) for x in pc]) # if there are multiple percentage points, take the larger one 
        else:
            ki67=None

