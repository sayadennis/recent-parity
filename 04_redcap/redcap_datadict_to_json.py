import numpy as np
import pandas as pd
import json

dn = '/share/fsmresfiles/breast_cancer_pregnancy/'
datadict_fn = 'FrequencyAndResultsOfGeneticTe_DataDictionary_2022-09-26.csv'

datadict_df = pd.read_csv(f'{dn}/{datadict_fn}')
datadict_df = datadict_df.set_index('Variable / Field Name')
datadict_df.index.name = None

dd = datadict_df.to_dict(orient='index')

for key in dd.keys(): # key is the column nane 
    if isinstance(dd[key]['Choices, Calculations, OR Slider Labels'], str): # if multiple choice
        choices = dd[key]['Choices, Calculations, OR Slider Labels']
        choices = choices.split(' | ') # element like "1, Pathogenic"
        choices = [item.split(', ') for item in choices] # element like ["1", "Pathogenic"]
        choices = {int(items[0]): items[1] for items in choices} # element like {1: "Pathogenic"} 
        dd[key]['Choices, Calculations, OR Slider Labels'] = choices

####################################
#### Tables for Andrea's report ####
####################################

## Items to include (note)
# Demographics tables - like you mentioned, age 
# Parity distribution - independent variable, so separating into those 4 categories 
# NAC regimen distribution
# Common mutations
# Biomarker stages - three groups (ER/PR+ HER2-; triple negative (no hormone receptor no HER2); and HER2+) distribution



