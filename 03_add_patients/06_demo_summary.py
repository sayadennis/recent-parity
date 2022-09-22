import os
import sys
import numpy as np
import pandas as pd
from sklearn.utils import resample
from sklearn.linear_model import LogisticRegression

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

###################
#### Read data ####
###################

data_dict = pd.read_csv(f'{dn}/recent_parity_data_dictionary.csv')

data = pd.read_csv(f'{dn}/tkhr_pat_wo_gentest_added/recent_parity_combined_data_reviewed_tkhr_added.csv')

## Calculate years since last pregnancy


######################
#### Data summary ####
######################

print('###################################################')
print('#### Recent parity project demographic summary ####')
print('###################################################\n')

birth_years = []
for birth_date in data['DOB']:
    if '/' in birth_date:
        yr_string = birth_date[-2:]
        if yr_string.startswith('0'):
            yr = int(yr_string[1])
        else:
            yr = int(yr_string)
        #
        if yr > 3:
            yr = int('19' + yr_string)
        else:
            yr = int('20' + yr_string)
        # 
        birth_years.append(yr)
    elif '-' in birth_date:
        yr = int(birth_date.split('-')[0])
        birth_years.append(yr)
    else:
        print('birth date string does not contain / or -: ', birth_date)
#

print('## Birth date summary statistics ##')
print(f'Birth year range: {min(birth_years)}-{max(birth_years)}')
print(f'Mean {np.mean(birth_years):.2f}')
print(f'Median {np.median(birth_years):.2f}')
print(f'Standard deviation {np.std(birth_years):.2f}\n')

for item in ['Gravida', 'Para', 'TimeBCDx', 'AgeBCDx', 'AgeFirstPreg', 'AgeLastPreg']: # 'Menarche', 
    print(f'## {item} summary statistics ##')
    values = list(data[item])
    print(f'{item} range: {min(values)}-{max(values)}')
    print(f'Mean {np.nanmean(values):.2f}')
    print(f'Median {np.nanmedian(values):.2f}')
    print(f'Standard deviation {np.nanstd(values):.2f}\n')

