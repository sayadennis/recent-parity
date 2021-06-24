import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

dn = '/share/fsmresfiles/breast_cancer_pregnancy/notes'

with open(os.path.join(dn, 'column_names_note_type_detail.txt'), 'r') as f:
    lines = f.readlines()

colnames = []
for line in lines:
    colnames.append(line.strip())


#### Genetic counseling notes #### 

struc = pd.read_csv(
    os.path.join(dn, 'genetic_counsel/notes_genetic_counsel_notetypedetail_missingtext.csv'), header=None
)
struc.columns = colnames[:12]

# textlist = []
# with open(os.path.join(dn, 'genetic_counsel/notes_genetic_counsel_notetypedetail_onlytext.txt'), 'r') as f:
#     lines = f.readlines()

# for line in lines:
#     textlist.append(line.strip())

counts = pd.DataFrame(
    0, dtype=int, 
    index=[
        'RTF & plain (NMFF Clarity)', 'RTF only (NMFF Clarity)', 'plain only (NMFF Clarity)', 
        'RTF & plain (Clarity Foundation)', 'RTF only (Clarity Foundation)', 'plain only (Clarity Foundation)'
    ], 
    columns=['total', 'pre_project1', 'post_project1']
)

for i in struc['source_system_id'].unique():
    substruc = struc.iloc[struc['source_system_id'].values==i,:]
    if len(substruc)==1: # exists in only 1 form
        # investicate this one case 
        ix = substruc.index
        rtf = (substruc.loc[ix,'source_system_table']=='y_rtf_note_text').values[0] # boolean
        system_nmff = (substruc.loc[ix,'source_system_key']==2).values[0] # boolean
        pre_proj1 = datetime.strptime(substruc.loc[ix,'created_datetime'].values[0], '%Y-%m-%d %H:%M:%S.%f') < datetime(2018,3,1)
        if (rtf & system_nmff):
            counts.loc['RTF only (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF only (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['RTF only (NMFF Clarity)', 'post_project1'] += 1
        elif (~rtf & system_nmff):
            counts.loc['plain only (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['plain only (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['plain only (NMFF Clarity)', 'post_project1'] += 1
        elif (~rtf & ~system_nmff):
            counts.loc['plain only (Clarity Foundation)', 'total'] += 1
            if pre_proj1:
                counts.loc['plain only (Clarity Foundation)', 'pre_project1'] += 1
            else:
                counts.loc['plain only (Clarity Foundation)', 'post_project1'] += 1
        else:
            print('found RTF in Clarity Foundation:\n')
            print('source system table is {} and source system key is {}'.format(struc.loc[i,'source_system_table'], struc.loc[i,'source_system_key']))
    else: # exists in 2+ forms
        ix_rtf = substruc.iloc[substruc['source_system_table'].values=='y_rtf_note_text',:].index
        ix_plain = substruc.iloc[substruc['source_system_table'].values=='hno_note_text',:].index
        system_nmff = np.all(substruc.loc[ix_rtf, 'source_system_key']==2) #.values[0] # boolean
        pre_proj1 = datetime.strptime(substruc.loc[ix_rtf, 'created_datetime'].values[0], '%Y-%m-%d %H:%M:%S.%f') < datetime(2018,3,1)
        if system_nmff:
            counts.loc['RTF & plain (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF & plain (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['RTF & plain (NMFF Clarity)', 'post_project1'] += 1
        else:
            counts.loc['RTF & plain (Clarity Foundation)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF & plain (Clarity Foundation)', 'pre_project1'] += 1
            else:
                counts.loc['RTF & plain (Clarity Foundation)', 'post_project1'] += 1


counts.to_csv(os.path.join(dn, 'genetic_counsel/text_format_counts_summary.csv'), header=True, index=True)



#### Initial counsel notes #### 

struc = pd.read_csv(
    os.path.join(dn, 'initial_counsel/notes_initial_counsel_notetypedetail_missingtext.csv'), header=None
)
struc.columns = colnames[:12]
struc = struc.iloc[[x in ['hno_note_text', 'y_rtf_note_text'] for x in struc['source_system_table'].values],:]

# textlist = []
# with open(os.path.join(dn, 'initial_counsel/notes_initial_counsel_notetypedetail_onlytext.txt'), 'r') as f:
#     lines = f.readlines()

# for line in lines:
#     textlist.append(line.strip())


counts = pd.DataFrame(
    0, dtype=int, 
    index=[
        'RTF & plain (NMFF Clarity)', 'RTF only (NMFF Clarity)', 'plain only (NMFF Clarity)', 
        'RTF & plain (Clarity Foundation)', 'RTF only (Clarity Foundation)', 'plain only (Clarity Foundation)'
    ], 
    columns=['total', 'pre_project1', 'post_project1']
)

for i in struc['source_system_id'].unique():
    substruc = struc.iloc[struc['source_system_id'].values==i,:]
    if len(substruc)==1: # exists in only 1 form
        # investicate this one case 
        ix = substruc.index
        rtf = (substruc.loc[ix,'source_system_table']=='y_rtf_note_text').values[0] # boolean
        system_nmff = (substruc.loc[ix,'source_system_key']==2).values[0] # boolean
        pre_proj1 = datetime.strptime(substruc.loc[ix,'created_datetime'].values[0], '%Y-%m-%d %H:%M:%S.%f') < datetime(2018,3,1)
        if (rtf & system_nmff):
            counts.loc['RTF only (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF only (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['RTF only (NMFF Clarity)', 'post_project1'] += 1
        elif (~rtf & system_nmff):
            counts.loc['plain only (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['plain only (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['plain only (NMFF Clarity)', 'post_project1'] += 1
        elif (~rtf & ~system_nmff):
            counts.loc['plain only (Clarity Foundation)', 'total'] += 1
            if pre_proj1:
                counts.loc['plain only (Clarity Foundation)', 'pre_project1'] += 1
            else:
                counts.loc['plain only (Clarity Foundation)', 'post_project1'] += 1
        else:
            print('found RTF in Clarity Foundation:\n')
            print('source system table is {} and source system key is {}'.format(struc.loc[i,'source_system_table'], struc.loc[i,'source_system_key']))
    else: # exists in 2+ forms
        ix_rtf = substruc.iloc[substruc['source_system_table'].values=='y_rtf_note_text',:].index
        ix_plain = substruc.iloc[substruc['source_system_table'].values=='hno_note_text',:].index
        system_nmff = np.all(substruc.loc[ix_rtf, 'source_system_key']==2) #.values[0] # boolean
        pre_proj1 = datetime.strptime(substruc.loc[ix_rtf, 'created_datetime'].values[0], '%Y-%m-%d %H:%M:%S.%f') < datetime(2018,3,1)
        if system_nmff:
            counts.loc['RTF & plain (NMFF Clarity)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF & plain (NMFF Clarity)', 'pre_project1'] += 1
            else:
                counts.loc['RTF & plain (NMFF Clarity)', 'post_project1'] += 1
        else:
            counts.loc['RTF & plain (Clarity Foundation)', 'total'] += 1
            if pre_proj1:
                counts.loc['RTF & plain (Clarity Foundation)', 'pre_project1'] += 1
            else:
                counts.loc['RTF & plain (Clarity Foundation)', 'post_project1'] += 1


counts.to_csv(os.path.join(dn, 'initial_counsel/text_format_counts_summary.csv'), header=True, index=True)
