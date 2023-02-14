import pickle
import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy'

###################
#### Load data ####
###################

datadir = 'data/06_exported_from_redcap'
fn = 'FrequencyAndResultsO_DATA_2022-10-20_1504.csv'
redcap = pd.read_csv(f'{dn}/{datadir}/{fn}')

with open(f'{dn}/{datadir}/data_dictionary.p', 'rb') as f:
    dd = pickle.load(f)
