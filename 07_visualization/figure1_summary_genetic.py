import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
dout = "/share/fsmresfiles/breast_cancer_pregnancy/plots"
fin = "20230509_super_summary_table.csv"
fout = "figure1_genetic.png"

data = pd.read_csv(f"{din}/{fin}")

genes = [
    'BRCA1', 
    'BRCA2', 
    'PALB2', 
    'CHEK2', 
    'ATM',
    'TP53', 
    'PTEN', 
    'STK11', 
    'CDH1', 
]
genes.reverse() # reverse so that the horizontal plot starts from BRCA1

cts = {}
for category in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']:
    cts[category] = [ 
        float(data.iloc[data['main_category'].values==gene,:][category].values[0].split()[0])
        for gene in genes
    ]   

fig, ax = plt.subplots(figsize=(10,3.5))
pastel_colors = list(plt.cm.Pastel1.colors[:4])

ax.barh(
    range(len(genes)),
    cts['Nulliparous'],
    label='Nulliparous',
    color=pastel_colors[0],
    edgecolor='black',
    linewidth=.5,
)
ax.barh(
    range(len(genes)),
    cts['<5 years'],
    left=cts['Nulliparous'],
    label='<5 years',
    color=pastel_colors[1],
    edgecolor='black',
    linewidth=.5,
)
ax.barh(
    range(len(genes)),
    cts['5-10 years'],
    left=[np.sum([x,y]) for x,y in zip(cts['Nulliparous'], cts['<5 years'])],
    label='5-10 years',
    color=pastel_colors[2],
    edgecolor='black',
    linewidth=.5,
)
ax.barh(
    range(len(genes)),
    cts['>=10 years'],
    left=[np.sum([x,y,z]) for x,y,z in zip(cts['Nulliparous'], cts['<5 years'], cts['5-10 years'])],
    label='>=10 years',
    color=pastel_colors[3],
    edgecolor='black',
    linewidth=.5,
)
ax.set_xlabel('Number of patients')
ax.set_title('Genetic Testing Results Summary')
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes)
ax.legend(loc='lower right')

# Add numbers at the top of each bar
totals = [sum([cts[key][i] for key in ['Nulliparous', '<5 years', '5-10 years', '>=10 years']]) for i in range(len(genes))]
for index, value in enumerate(totals):
    ax.text(value + 0.5, index, str(int(value)), va='center')

## Remove the top and right spines for all subplots
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(f"{dout}/{fout}")
plt.close()

