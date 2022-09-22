import os
import numpy as np
import pandas as pd

dn = '/share/fsmresfiles/breast_cancer_pregnancy/data'

data = pd.read_csv(f'{dn}/recent_parity_combined_data_reviewed.csv')
reviewed = data.iloc[data['Review complete'].values=='x',:]

# genes = ['general', 'brca', 'brca1', 'brca2', 'palb2', 'tp53', 'pten', 'cdh1', 'stk11', 'chek2', 'atm']

# fix the value of age of first pregnancy where it's zero 
reviewed.loc[reviewed.index[reviewed['agefirstpreg'].values==0],'agefirstpreg'] = np.nan

# fix age of menarche values that are off 
reviewed.loc[reviewed.index[reviewed['menarche'].values>20],'menarche'] = np.nan

print('###################################################')
print('#### Recent parity project demographic summary ####')
print('###################################################\n')

birth_years = []
for birth_date in reviewed['birth_date']:
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
#
print('## Birth date summary statistics ##')
print(f'Birth year range: {min(birth_years)}-{max(birth_years)}')
print(f'Mean {np.mean(birth_years):.2f}')
print(f'Median {np.median(birth_years):.2f}')
print(f'Standard deviation {np.std(birth_years):.2f}\n')

for item in ['gravida', 'para', 'menarche', 'diagnosis_year', 'age_at_diagnosis', 'agefirstpreg']:
    print(f'## {item} summary statistics ##')
    values = list(reviewed[item])
    print(f'{item} range: {min(values)}-{max(values)}')
    print(f'Mean {np.nanmean(values):.2f}')
    print(f'Median {np.nanmedian(values):.2f}')
    print(f'Standard deviation {np.nanstd(values):.2f}\n')

