import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
fin = "20230509_super_summary_table.csv"
fout = "figure1_demo_gyn_tumor.png"

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
    linewidth=.5,
    error_kw=dict(ecolor='black', elinewidth=.5)
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
    linewidth=.5,
    error_kw=dict(ecolor='black', elinewidth=.5)
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
    linewidth=.5,
    error_kw=dict(ecolor='black', elinewidth=.5)
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

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        float(data.iloc[data['sub_category'].values==bm,:][category].values[0].split()[0])
        for bm in ['ER/PR+ HER2-', 'HER2+', 'Triple Negative']
    ]

# Comment the below chunk if we want a stacked bar graph (of counts) as opposed to a stacked percentage chart
pcts = {}
for key in cts.keys():
    pcts[key] = []
    for i in range(len(cts['Nulliparous'])):
        pcts[key].append(100 * cts[key][i]/np.sum([cts[x][i] for x in cts.keys()]))

ax[1,0].bar(
    [1,2,3],
    pcts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
    linewidth=.5,
)
ax[1,0].bar(
    [1,2,3],
    pcts['<5 years'],
    bottom=pcts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
    linewidth=.5,
)
ax[1,0].bar(
    [1,2,3],
    pcts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(pcts['Nulliparous'], pcts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
    linewidth=.5,
)
ax[1,0].bar(
    [1,2,3],
    pcts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(pcts['Nulliparous'], pcts['<5 years'], pcts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
    linewidth=.5,
)
ax[1,0].set_ylabel('Percentage of patients')
ax[1,0].set_title('Biomarker category')
ax[1,0].set_xticks([1,2,3])
ax[1,0].set_ylim([0.,100.])
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

# Comment the below chunk if we want a stacked bar graph (of counts) as opposed to a stacked percentage chart
pcts = {}
for key in cts.keys():
    pcts[key] = []
    for i in range(len(cts['Nulliparous'])):
        pcts[key].append(100 * cts[key][i]/np.sum([cts[x][i] for x in cts.keys()]))

ax[1,1].bar(
    [1,2,3],
    pcts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
    linewidth=.5,
)
ax[1,1].bar(
    [1,2,3],
    pcts['<5 years'],
    bottom=pcts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
    linewidth=.5,
)
ax[1,1].bar(
    [1,2,3],
    pcts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(pcts['Nulliparous'], pcts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
    linewidth=.5,
)
ax[1,1].bar(
    [1,2,3],
    pcts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(pcts['Nulliparous'], pcts['<5 years'], pcts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
    linewidth=.5,
)
ax[1,1].set_ylabel('Number of patients')
ax[1,1].set_title('Histology')
ax[1,1].set_xticks([1,2,3])
ax[1,1].set_ylim([0.,100.])
ax[1,1].set_xticklabels(['Invasive', 'Non-invasive', 'Other'], rotation=45, ha='right')
# ax[1,1].legend()

## Histologic Grade

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [
        np.sum([int(x.split()[0]) for x in  data.iloc[
            np.array([x in varnames for x in data['sub_category'].values]) & 
            np.array([data['main_category'].values=='histologic_grade']).ravel(),:
        ][category].values])
        for varnames in [['1'],['2'],['3']]
    ]

# Comment the below chunk if we want a stacked bar graph (of counts) as opposed to a stacked percentage chart
pcts = {}
for key in cts.keys():
    pcts[key] = []
    for i in range(len(cts['Nulliparous'])):
        pcts[key].append(100 * cts[key][i]/np.sum([cts[x][i] for x in cts.keys()]))

ax[1,2].bar(
    [1,2,3],
    pcts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
    linewidth=.5,
)
ax[1,2].bar(
    [1,2,3],
    pcts['<5 years'],
    bottom=pcts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
    linewidth=.5,
)
ax[1,2].bar(
    [1,2,3],
    pcts['5-10 years'],
    bottom=[np.sum([x,y]) for x,y in zip(pcts['Nulliparous'], pcts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
    linewidth=.5,
)
ax[1,2].bar(
    [1,2,3],
    pcts['>=10 years'],
    bottom=[np.sum([x,y,z]) for x,y,z in zip(pcts['Nulliparous'], pcts['<5 years'], pcts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
    linewidth=.5,
)
ax[1,2].set_ylabel('Number of patients')
ax[1,2].set_title('Histologic Grade')
ax[1,2].set_xticks([1,2,3])
ax[1,2].set_ylim([0.,100.])
ax[1,2].set_xticklabels(['Grade I', 'Grade II', 'Grade III'], rotation=45, ha='right')

## Remove the top and right spines for all subplots
for row in ax:
    for axis in row:
        axis.spines['right'].set_visible(False)

## Remove the top spines only for the 3 plots in the top row
for axis in ax[0]:
    axis.spines['top'].set_visible(False)

plt.tight_layout()
fig.savefig(f"/share/fsmresfiles/breast_cancer_pregnancy/plots/{fout}")
plt.close()
