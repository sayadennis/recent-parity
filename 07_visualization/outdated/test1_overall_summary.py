import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
fin = "20230509_super_summary_table.csv"
fout = "20230822_test1_summary_plot.png"

data = pd.read_csv(f"{din}/{fin}")

##################################################
#### Visualization 1: all parallel bar graphs ####
##################################################

fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(8,10))

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
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
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
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
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
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
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

for i, varname in enumerate(['ER/PR+ HER2-', 'HER2+', 'Triple Negative']): # , 'Missing'
    cts = [
        float(data.iloc[data['sub_category'].values==varname,:][category].values[0].split()[0]) 
        for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
    ]

    ax[1, i].bar(
        [1,3,4,5,6],
        cts,
        color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
        edgecolor='black',
        width=1.0,
    )

    ax[1, i].set_ylim(0, None)
    ax[1, i].set_xticks([1,3,4,5,6])
    ax[1, i].set_xticklabels(
        ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
        rotation=40, ha='right'
    )
    ax[1, i].set_ylabel('Patient counts')
    ax[1, i].set_title(varname)

## Histology - invasive 

cts = {}
for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histology']).ravel(),:
        ][category].values])
        for varnames in [['IDC', 'ILC', 'Invasive (non-specific)'],['DCIS'],['Other']]
    ]


ax[2, 0].bar(
    [1,3,4,5,6],
    [cts[category][0] for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']],
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[2, 0].set_ylim(0, None)
ax[2, 0].set_xticks([1,3,4,5,6])
ax[2, 0].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[2, 0].set_ylabel('Patient counts')
ax[2, 0].set_title('Invasive')

## Histology - non-invasive 

ax[2, 1].bar(
    [1,3,4,5,6],
    [cts[category][1] for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']],
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[2, 1].set_ylim(0, None)
ax[2, 1].set_xticks([1,3,4,5,6])
ax[2, 1].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[2, 1].set_ylabel('Patient counts')
ax[2, 1].set_title('Non-invasive')

## Histology - other

ax[2, 2].bar(
    [1,3,4,5,6],
    [cts[category][2] for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']],
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[2, 2].set_ylim(0, None)
ax[2, 2].set_xticks([1,3,4,5,6])
ax[2, 2].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[2, 2].set_ylabel('Patient counts')
ax[2, 2].set_title('Other')

## Histologic Grade - 1

varname = '1'

cts = [
    float(
        data.iloc[
                np.array(data['sub_category'].values==varname) & 
                np.array(data['main_category'].values=='histologic_grade'),:
            ][category].values[0].split()[0]
    ) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

ax[3, 0].bar(
    [1,3,4,5,6],
    cts,
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[3, 0].set_ylim(0, None)
ax[3, 0].set_xticks([1,3,4,5,6])
ax[3, 0].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[3, 0].set_ylabel('Patient counts')
ax[3, 0].set_title('Grade I')

## Histologic Grade - 2

varname = '2'

cts = [
    float(
        data.iloc[
                np.array(data['sub_category'].values==varname) & 
                np.array(data['main_category'].values=='histologic_grade'),:
            ][category].values[0].split()[0]
    ) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

ax[3, 1].bar(
    [1,3,4,5,6],
    cts,
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[3, 1].set_ylim(0, None)
ax[3, 1].set_xticks([1,3,4,5,6])
ax[3, 1].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[3, 1].set_ylabel('Counts')
ax[3, 1].set_title('Grade II')

## Histologic Grade - 3

varname = '3'

cts = [
    float(
        data.iloc[
                np.array(data['sub_category'].values==varname) & 
                np.array(data['main_category'].values=='histologic_grade'),:
            ][category].values[0].split()[0]
    ) 
    for category in ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years']
]

ax[3, 2].bar(
    [1,3,4,5,6],
    cts,
    color=['gray'] + list(plt.cm.Pastel1.colors[:4]),
    edgecolor='black',
    width=1.0,
)

ax[3, 2].set_ylim(0, None)
ax[3, 2].set_xticks([1,3,4,5,6])
ax[3, 2].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 years', '5-10 years', '>=10 years'],
    rotation=40, ha='right'
)
ax[3, 2].set_ylabel('Counts')
ax[3, 2].set_title('Grade III')


## Remove the top and right spines for all subplots
for row in ax:
    for axis in row:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"/share/fsmresfiles/breast_cancer_pregnancy/plots/{fout}")
plt.close()
