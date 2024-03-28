# pylint: disable=C0114
# pylint: disable=W0621
# pylint: disable=R0914
# pylint: disable=duplicate-code

import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

dn = "/share/fsmresfiles/breast_cancer_pregnancy"
dout = "/share/fsmresfiles/breast_cancer_pregnancy/plots"

###################
#### Read data ####
###################

datadir = "data/06_exported_from_redcap"
fn = "FrequencyAndResultsO_DATA_2023-03-17_0949.csv"

data = pd.read_csv(f"{dn}/{datadir}/{fn}")

# remove exluded patients
data = data.iloc[(data.exclude_demo.values != 1) & (data.exclude_tum.values != 1), :]

with open(f"{dn}/{datadir}/data_dictionary.p", "rb") as f:
    dd = pickle.load(f)

# mark family history as none for patients who are missing it
data.fam_hx = data.fam_hx.fillna(0.0)
data.parous = data.parous.fillna(0.0)

#########################################
#### Fill in useful categorical info ####
#########################################

genes = [
    "brca1",
    "brca2",
    "palb2",
    "chek2",
    "atm",
    # "tp53",
    # "pten",
    # "stk11",
    # "cdh1",
]

# Change mutation values to binary encoding
for gene in genes:
    data[gene] = (
        data[gene].map(dd[gene]["Choices, Calculations, OR Slider Labels"])
        == "Pathogenic"
    ).astype(int)

data_parous = data.iloc[data.parous.values == 1, :]

############################
#### Plot distributions ####
############################

plot_features = {
    "Age at last parity": "age_at_most_recent_pregnancy",
    "Age at breast cancer Dx": "age_at_diagnosis",
    "Years between parity/cancer": "years_since_pregnancy",
}

colors = {
    "All patients": "tab:gray",
    "brca1": "c",
    "brca2": "m",
    "palb2": "darkorange",
    "chek2": "mediumseagreen",
    "atm": "mediumpurple",
}

fig, axs = plt.subplots(
    nrows=len(genes) + 1,  # overall + genes
    ncols=len(plot_features.items()),
    figsize=(10, 8),
)

for j, (feature_name, varname) in enumerate(plot_features.items()):
    for i, gene in enumerate(["All patients"] + genes):
        axs[i, j].spines["top"].set_visible(False)
        axs[i, j].spines["right"].set_visible(False)
        if feature_name == "Years between parity/cancer":
            axs[i, j].set_xlim(0, 45)
            bins = np.arange(0, 36, 5)
        else:
            axs[i, j].set_xlim(15, 51)
            bins = np.arange(15, 51, 5)
        if i == 0:  # if overall distribution
            axs[i, j].hist(data_parous[varname], bins=bins, color=colors[gene])
            # axs[i,j].set_title(feature_name, fontsize=14)
        else:
            axs[i, j].hist(
                data_parous.iloc[data_parous[gene].values == 1, :][varname],
                bins=bins,
                color=colors[gene],
            )
            if data_parous[gene].sum() < 6:
                axs[i, j].set_ylim(0, 3.1)
                axs[i, j].set_yticks([0, 3])
                axs[i, j].set_yticklabels([0, 3])
            elif data_parous[gene].sum() == 0:
                axs[i, j].set_ylim(0, 1.1)
                axs[i, j].set_yticks([0, 1])
                axs[i, j].set_yticklabels([0, 1])

# Add row labels with genes
fig.text(0.03, 0.90, "All patients", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.74, "BRCA1", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.58, "BRCA2", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.43, "PALB2", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.27, "CHEK2", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.10, "ATM", rotation=90, va="center", ha="center", fontsize=14)

# Add column labels with plot titles
fig.text(0.23, 0.98, "Age at last pregnancy", va="center", ha="center", fontsize=14)
fig.text(0.55, 0.98, "Age at breast cancer Dx", va="center", ha="center", fontsize=14)
fig.text(
    0.85, 0.98, "Years between parity/cancer", va="center", ha="center", fontsize=14
)

plt.tight_layout()
plt.subplots_adjust(left=0.1, top=0.95)
fig.savefig(f"{dout}/age_distributions_by_mutation.png")
plt.close()
