import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
fin = "20230509_super_summary_table.csv"
fout = "20230822_test2_summary_plot.png"

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
    ['Overall', 'Nulliparous', '<5 yrs', '5-10 yrs', '>=10 yrs'],
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
    ['Overall', 'Nulliparous', '<5 yrs', '5-10 yrs', '>=10 yrs'],
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
for item, label in zip(lines, ['Overall', 'Nulliparous', '<5 yrs', '5-10 yrs', '>=10 yrs']):
    item.set_label(label)

ax[0, 2].set_xticks([1,3,4,5,6])
ax[0, 2].set_xticklabels(
    ['Overall', 'Nulliparous', '<5 yrs', '5-10 yrs', '>=10 yrs'],
    rotation=40, ha='right'
)
ax[0, 2].set_ylim(0, None)
ax[0, 2].set_ylabel('Mean Para')
ax[0, 2].set_title('Para')

## Biomarker

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        float(data.iloc[data['sub_category'].values==bm,:][category].values[0].split()[0])
        for bm in ['ER/PR+ HER2-', 'HER2+', 'Triple Negative']
    ]

ax[1,0].bar(
    [1,2,3],
    cts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
)
ax[1,0].bar(
    [1,2,3],
    cts['<5 years'],
    bottom=cts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
)
ax[1,0].bar(
    [1,2,3],
    cts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['Nulliparous'], cts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
)
ax[1,0].bar(
    [1,2,3],
    cts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(cts['Nulliparous'], cts['<5 years'], cts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
)
ax[1,0].set_ylabel('Number of patients')
ax[1,0].set_title('Biomarker category')
ax[1,0].set_xticks([1,2,3])
ax[1,0].set_xticklabels(['ER/PR+ HER2-', 'HER2+', 'Triple Negative'], rotation=45, ha='right')
# ax[1,0].legend()

## Histology - invasive 

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histology']).ravel(),:
        ][category].values])
        for varnames in [['IDC', 'ILC', 'Invasive (non-specific)'],['DCIS'],['Other']]
    ]

ax[1,1].bar(
    [1,2,3],
    cts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
)
ax[1,1].bar(
    [1,2,3],
    cts['<5 years'],
    bottom=cts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
)
ax[1,1].bar(
    [1,2,3],
    cts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['Nulliparous'], cts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
)
ax[1,1].bar(
    [1,2,3],
    cts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(cts['Nulliparous'], cts['<5 years'], cts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
)
ax[1,1].set_ylabel('Number of patients')
ax[1,1].set_title('Histology')
ax[1,1].set_xticks([1,2,3])
ax[1,1].set_xticklabels(['Invasive', 'Non-invasive', 'Other'], rotation=45, ha='right')
# ax[1,1].legend()

## Histologic Grade - 1

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histologic_grade']).ravel(),:
        ][category].values])
        for varnames in [['1'],['2'],['3']]
    ]


ax[1,2].bar(
    [1,2,3],
    cts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
)
ax[1,2].bar(
    [1,2,3],
    cts['<5 years'],
    bottom=cts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
)
ax[1,2].bar(
    [1,2,3],
    cts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(cts['Nulliparous'], cts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
)
ax[1,2].bar(
    [1,2,3],
    cts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(cts['Nulliparous'], cts['<5 years'], cts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
)
ax[1,2].set_ylabel('Number of patients')
ax[1,2].set_title('Histologic Grade')
ax[1,2].set_xticks([1,2,3])
ax[1,2].set_xticklabels(['Grade I', 'Grade II', 'Grade III'], rotation=45, ha='right')

## Remove the top and right spines for all subplots
for row in ax:
    for axis in row:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"/share/fsmresfiles/breast_cancer_pregnancy/plots/{fout}")
plt.close()
