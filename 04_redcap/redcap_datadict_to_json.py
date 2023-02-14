import numpy as np
import pandas as pd
import pickle

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data/06_exported_from_redcap'

datadict_fn = 'FrequencyAndResultsOfGeneticTe_DataDictionary_2022-12-09.csv'

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

with open(f'{dn}/data_dictionary.p', 'wb') as f:
    pickle.dump(dd, f)
