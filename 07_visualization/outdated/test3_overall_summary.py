import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
fin = "20230509_super_summary_table.csv"
fout = "20230822_test3_summary_plot.png"

data = pd.read_csv(f"{din}/{fin}")

#########################
#### Visualization 2 ####
#########################

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10,7))

pastel_colors = list(plt.cm.Pastel1.colors[:4])

## Demographic - age 

varname = 'Mean age (±SD)'
means = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[0]) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]
sds = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[1].strip('()±'))
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

ax[0, 0].bar(
    [1,3,4,5,6],
    means,
    yerr=sds, capsize=0,
    color=['gray'] + pastel_colors,
    edgecolor='black',
    width=1.0,
)

ax[0, 0].set_ylim(0, None)
ax[0, 0].set_xticks([1,3,4,5,6])
ax[0, 0].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[0, 0].set_ylabel('Mean Age')
ax[0, 0].set_title('Age at diagnosis')

## Demographic - gravida 

varname = 'Mean Gravida (±SD)'
means = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[0]) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]
sds = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[1].strip('()±'))
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

ax[0, 1].bar(
    [1,3,4,5,6],
    means,
    yerr=sds, capsize=0,
    color=['gray'] + pastel_colors,
    edgecolor='black',
    width=1.0,
)

ax[0, 1].set_xticks([1,3,4,5,6])
ax[0, 1].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[0, 1].set_ylim(0, None)
ax[0, 1].set_ylabel('Mean Gravida')
ax[0, 1].set_title('Gravida')

## Demographic - para 

varname = 'Mean Para (±SD)'
means = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[0]) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]
sds = [
    float(data.iloc[data['main_category'].values==varname,:][category].values[0].split()[1].strip('()±'))
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

lines = ax[0, 2].bar(
    [1,3,4,5,6],
    means,
    yerr=sds, capsize=0,
    color=['gray'] + pastel_colors,
    edgecolor='black',
    width=1.0,
)
for item, label in zip(lines, ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']):
    item.set_label(label)

ax[0, 2].set_xticks([1,3,4,5,6])
ax[0, 2].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[0, 2].set_ylim(0, None)
ax[0, 2].set_ylabel('Mean Para')
ax[0, 2].set_title('Para')

## Biomarker

set_colors = list(plt.cm.Set2.colors[:3])

cts = {}
for bm in ['ER/PR+ HER2-', 'HER2+', 'Triple Negative']:
    cts[bm] = [
        float(data.iloc[data['sub_category'].values==bm,:][category].values[0].split()[0])
        for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ]

ax[1,0].bar(
    [1,3,4,5,6],
    cts['ER/PR+ HER2-'],
    label='ER/PR+ HER2-',
    color=set_colors[0],
    edgecolor='black',
)
ax[1,0].bar(
    [1,3,4,5,6],
    cts['HER2+'],
    bottom=cts['ER/PR+ HER2-'],
    label='HER2+',
    color=set_colors[1],
    edgecolor='black',
)
ax[1,0].bar(
    [1,3,4,5,6],
    cts['Triple Negative'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['ER/PR+ HER2-'], cts['HER2+'])],
    label='Triple Negative',
    color=set_colors[2],
    edgecolor='black',
)
ax[1,0].set_ylabel('Number of patients')
ax[1,0].set_title('Biomarker category')
ax[1,0].set_xticks([1,3,4,5,6])
ax[1,0].set_xticklabels(['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'], rotation=45, ha='right')
ax[1,0].legend(loc='upper right')

## Histology - invasive 

set_colors = list(plt.cm.Set2.colors[3:6])

cts = {}

varnames_dict = {
    'Invasive' : ['IDC', 'ILC', 'Invasive (non-specific)'],
    'Non-invasive' : ['DCIS'],
    'Other' : ['Other']
}

for key in varnames_dict.keys():
    cts[key] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames_dict[key] for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histology']).ravel(),:
        ][category].values])
        for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ]

ax[1,1].bar(
    [1,3,4,5,6],
    cts['Invasive'],
    label='Invasive',
    color=set_colors[0],
    edgecolor='black',
)
ax[1,1].bar(
    [1,3,4,5,6],
    cts['Non-invasive'],
    bottom=cts['Invasive'],
    label='Non-invasive',
    color=set_colors[1],
    edgecolor='black',
)
ax[1,1].bar(
    [1,3,4,5,6],
    cts['Other'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['Invasive'], cts['Non-invasive'])],
    label='Other',
    color=set_colors[2],
    edgecolor='black',
)
ax[1,1].set_ylabel('Number of patients')
ax[1,1].set_title('Histology')
ax[1,1].set_xticks([1,3,4,5,6])
ax[1,1].set_xticklabels(['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'], rotation=45, ha='right')
ax[1,1].legend()

## Histologic Grade - 1

set_colors = list(plt.cm.Set3.colors[-3:])

cts = {}
varnames_dict = {
    'Grade I' : ['1'],
    'Grade II' : ['2'],
    'Grade III' : ['3']
}

for key in varnames_dict.keys():
    cts[key] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames_dict[key] for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histologic_grade']).ravel(),:
        ][category].values])
        for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ]

ax[1,2].bar(
    [1,3,4,5,6],
    cts['Grade I'],
    label='Grade I',
    color=set_colors[0],
    edgecolor='black',
)
ax[1,2].bar(
    [1,3,4,5,6],
    cts['Grade II'],
    bottom=cts['Grade I'],
    label='Grade II',
    color=set_colors[1],
    edgecolor='black',
)
ax[1,2].bar(
    [1,3,4,5,6],
    cts['Grade III'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['Grade I'], cts['Grade II'])],
    label='Grade III',
    color=set_colors[2],
    edgecolor='black',
)
ax[1,2].set_ylabel('Number of patients')
ax[1,2].set_title('Histologic Grade')
ax[1,2].set_xticks([1,3,4,5,6])
ax[1,2].set_xticklabels(['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'], rotation=45, ha='right')
ax[1,2].legend()

## Remove the top and right spines for all subplots
for row in ax:
    for axis in row:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"/share/fsmresfiles/breast_cancer_pregnancy/plots/{fout}")
plt.close()
