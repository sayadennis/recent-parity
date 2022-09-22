import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din='/share/fsmresfiles/breast_cancer_pregnancy/explore_cohort'
dout='/share/fsmresfiles/breast_cancer_pregnancy/data_summary'

colnames=['ir_id', 'mrn1', 'mrn2', 'first_name', 'last_name', 'dob', 'location']

cohort_loc = pd.read_csv(f'{din}/location_genetic_counsel_notes_cohort.csv', header=None, names=colnames)

query_locs = list((cohort_loc['location'].unique())[~pd.isnull(cohort_loc['location'].unique())]) # unique location names that are not NULL 

for location in query_locs:
    nmh = (('NORTHWESTERN MEMORIAL HOSPITAL' in location) | ('NMH' in location))
    nmp = ('NMP' in location)
    nmg = ('NMG' in location)
    rmg = ('RMG' in location)
    dupage = (('DUPAGE' in location) | ('CDH' in location))
    lakeforest = (('LAKE FOREST' in location) | ('LFH' in location) | ('LFG' in location))
    valleywest = ('VALLEY WEST' in location)
    delnor = ('DELNOR' in location)
    kishwaukee = ('KISHWAUKEE' in location)
    grayslake = ('GRAYSLAKE' in location)
    zzcdpg = ('ZZCDPG' in location)
    if not (nmh | nmp | nmg | rmg | dupage | lakeforest | valleywest | delnor | kishwaukee | grayslake | zzcdpg):
        print(location)

loclist = ['NMH', 'NMP', 'NMG', 'RMG', 'CDH', 'LFH', 'VW', 'DELN', 'KISH', 'ZZCDPG', 'Other']

def get_loc(location_string):
    if (('NORTHWESTERN MEMORIAL HOSPITAL' in location) | ('NMH' in location)):
        loc = 'NMH'
    elif ('NMP' in location):
        loc = 'NMP'
    elif ('NMG' in location):
        loc = 'NMG'
    elif ('RMG' in location):
        loc = 'RMG'
    elif (('DUPAGE' in location) | ('CDH' in location)):
        loc = 'CDH'
    elif (('LAKE FOREST' in location) | ('LFH' in location) | ('LFG' in location)):
        loc = 'LFH'
    elif ('VALLEY WEST' in location):
        loc = 'VW'
    elif ('DELNOR' in location):
        loc = 'DELN'
    elif ('KISHWAUKEE' in location):
        loc = 'KISH'
    elif ('ZZCDPG' in location):
        loc = 'ZZCDPG'
    else:
        loc = 'Other'
    return loc

locations = pd.DataFrame(0, index=cohort_loc['ir_id'].unique(), columns=loclist)

for ir_id in cohort_loc['ir_id'].unique():
    patient_locs = list(cohort_loc.iloc[cohort_loc['ir_id'].values==ir_id, [x=='location' for x in cohort_loc.columns]].dropna().values.T[0])
    for location in patient_locs:
        loc = get_loc(location)
        locations.loc[ir_id, loc] = 1

fig, ax = plt.subplots()
ax.pie(locations.sum().sort_values(), autopct='%1.1f%%', shadow=True, startangle=90)
ax.legend(locations.sum().sort_values().index, loc='center right')
ax.axis('equal')

fig.savefig(f'{dout}/cohort_genetic_notes_piechart.png')
