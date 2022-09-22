import os
import numpy as np
import pandas as pd
import re
from datetime import datetime

din='/share/fsmresfiles/breast_cancer_pregnancy/data'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

colnames = [ # the alternative colnames for when focusing on note types etc. 
    'ir_id', 'EPIC_mrn', 'Powerchart_mrn', 'first_name', 'last_name', 'birth_date',
    'department_key', 'department_name', 'department_external_name', 'department_specialty', 'location_name',
    'created_year', 'updated_year', 'signed_year', 'source_system_table', 'note_type', 'note_status', 'note_text'
]

gennotes=pd.read_csv(f'{din}/genetic/notes_genetic_counsel.csv', header=None, names=colnames, index_col=None)
patholnotes=pd.read_csv(f'{din}/pathology/notes_pathology.csv', header=None, names=colnames, index_col=None)

# find out what note types exist 
gennotes['note_type'].value_counts().sort_values(ascending=False)
patholnotes['note_type'].value_counts().sort_values(ascending=False)

# find out what department specialty exist 
gennotes['department_specialty'].value_counts().sort_values(ascending=False)
patholnotes['department_specialty'].value_counts().sort_values(ascending=False)

# find out the distribution of note types and department specialty
grouped = pd.DataFrame(gennotes.groupby(['department_specialty', 'note_type'])['note_type'].count())
grouped.columns = ['note_counts']
print(grouped)
grouped.to_csv(f'{dout}/counts_note_type_department_specialty_genetic.csv', index=True, header=True)

grouped = pd.DataFrame(patholnotes.groupby(['department_specialty', 'note_type'])['note_type'].count())
grouped.columns = ['note_counts']
print(grouped)
grouped.to_csv(f'{dout}/counts_note_type_department_specialty_pathology.csv', index=True, header=True)
