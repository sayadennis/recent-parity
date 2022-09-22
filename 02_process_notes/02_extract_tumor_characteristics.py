import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din='/share/fsmresfiles/breast_cancer_pregnancy/data/notes/pathol'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary/exploratory'

## Get column names for pathology notes table
with open('/share/fsmresfiles/breast_cancer_pregnancy/data/notes/column_names_note_type_detail.txt', 'r') as f:
    lines = f.readlines()

notes_colnames=[]
for item in lines:
    notes_colnames.append(item.strip())

## Get column names for the table obtained from oncology tables 
onc_colnames = [
    'ir_id', 'accessioned_date', # 'patient_key'
    'pathologist', 'physician', 'resident', 'case_collect_date', 
    'age_at_case_collect', 'case_type', 'report_type', 'path_report_date', 
    'frozen_section_diagnosis', 'clinical_info', 'gross_description', 'diagnosis'
]

## Read in data
pathol = pd.read_csv(f'{din}/notes_pathology.csv', header=None, names=notes_colnames)
onc_pathol = pd.read_csv(f'{din}/pathology_from_onc_tables.csv', index_col=None, header=None, names=onc_colnames)

## Tumor characteristics: size
# "Largest single diameter:	1.4 cm"
# "Largest single diameter:	_"
# "Largest diameter:	0.08 cm (excision); 0.12 cm (core)"
# "MEASURING 1.1 CM IN DIAMETER." "LARGEST INDIVIDUAL DCIS DIAMETER 1.2 CM." "Largest single diameter: 1.2 cm"
# "Present on one slide (largest diameter):	N/A" "Largest single diameter:	2.5 cm on slide"
# "Largest single diameter:	1.0 cm" "Overall estimate of the size of DCIS:	2.6 cm"
# "Size/Extent of Lesion:		Present on one slide (largest diameter):	0.3 cm in greatest dimension (see note for part D)" "Largest single diameter:	N/A"
# "Largest single diameter:	2.2 cm"
# "Present	 Size/Extent of Lesion:		 	Present on one slide (largest diameter):	no" "Largest single diameter:	1.1 cm"
# "Largest single diameter: 6.0 CM (GROSSLY)"
# "DCIS focus measures 0.3 cm in greatest dimension" "Largest single diameter:	2.3 cm" "Overall estimate of the size of DCIS:	2.4 cm"
# "Overall estimate of the size of DCIS:	0.8 cm"
# "Largest single diameter:	0.8 cm in a single slide	 	Overall estimate of the size of DCIS:	Approximately 1.7 cm"
# "Largest single diameter:  1.0 cm  	Overall estimate of the size of DCIS (based on 	the thickness of the slices and the location of the involved blocks):  largest focus measuring 1.0 cm"
# "Largest single diameter:	0.4	 	Overall estimate of the size of DCIS:	At least 0.4 cm in greatest dimensio"

def get_tumorsize(diag):
    diag = re.sub('\s+', ' ', diag) # replace all whitespaces with single space
    size0 = re.findall(r'largest single diameter: ?[0-9]+.[0-9]+ cm', diag, re.IGNORECASE) # ,?\?? ?
    if len(size0)>0:
        size = re.findall(r'[\.0-9]+', size0[0])[0]
    else:
        size = 0.
    return size

size_df = pd.DataFrame(index=None, columns=['ir_id', 'date', 'size (cm)'])
for ir_id in onc_pathol['ir_id'].unique():
    patient_diags = list(onc_pathol['diagnosis'].iloc[onc_pathol['ir_id'].values==ir_id])
    diag_dates = list(onc_pathol['path_report_date'].iloc[onc_pathol['ir_id'].values==ir_id])
    for (diag, date) in zip(patient_diags, diag_dates):
        if type(diag)==str:
            size = get_tumorsize(diag)
            size_df = size_df.append({'ir_id' : ir_id, 'date' : date, 'size (cm)' : size}, ignore_index=True)
        else:
            continue

# how many pathology notes does a single person have in average? 
nums = []
for ir_id in onc_pathol['ir_id'].unique():
    patient_diags = list(onc_pathol['diagnosis'].iloc[onc_pathol['ir_id'].values==ir_id])
    nums.append(np.sum([type(x)==str for x in patient_diags]))

nums = []
for ir_id in onc_pathol['ir_id'].unique():
    patient_diags = list(onc_pathol['clinical_info'].iloc[onc_pathol['ir_id'].values==ir_id])
    nums.append(np.sum([type(x)==str for x in patient_diags]))

print('clinical info')
print('missing rate: ', np.sum(np.array(nums)==0)/len(nums))
print('')

nums = []
for ir_id in onc_pathol['ir_id'].unique():
    patient_diags = list(onc_pathol['gross_description'].iloc[onc_pathol['ir_id'].values==ir_id])
    nums.append(np.sum([type(x)==str for x in patient_diags]))

print('gross description')
print('missing rate: ', np.sum(np.array(nums)==0)/len(nums))


stsm=[]
for ir_id in pathol['ir_id'].unique():
    patient_notes = list(pathol['note_text'].iloc[pathol['ir_id'].values==ir_id])
    if np.any(['Staging Summary' in x for x in patient_notes]):
        stsm.append(ir_id)

onc=[]
for ir_id in onc_pathol['ir_id'].unique():
    patient_diags = list(onc_pathol['diagnosis'].iloc[onc_pathol['ir_id'].values==ir_id])
    if np.any([type(x)==str for x in patient_diags]):
        onc.append(ir_id)


ol=[]
for ir_id in stsm:
    if ir_id in onc:
        ol.append(ir_id)

len(ol)

from matplotlib_venn import venn2

venn2(subsets = (len(stsm)-len(ol), len(onc)-len(ol), len(ol)), set_labels = ('From general notes', 'From cancer registry'))
plt.savefig(f'{dout}/venn_tumor_characteristics_presence_notes_vs_onc_pathol.png')

## Tumor characteristics: tumor grade
## Tumor characteristics: histology

## Receptor status (ER, PR, HER2) 
## Ki-67 (low, intermediate, high)
## p53 mutation presence
## Lymph node status 
## Presence of lymphovascular invasion 
