# pylint: disable=missing-module-docstring
# pylint: disable=not-an-iterable
# pylint: disable=bare-except
# pylint: disable=not-an-iterable
# pylint: disable=C0330
# pylint: disable=simplifiable-if-statement


import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, ttest_ind

dn = "/share/fsmresfiles/breast_cancer_pregnancy"

## Load data & data dictionary
datadir = "data/06_exported_from_redcap"
fn = "FrequencyAndResultsO_DATA_2023-03-17_0949.csv"
data = pd.read_csv(f"{dn}/{datadir}/{fn}")

with open(f"{dn}/{datadir}/data_dictionary.p", "rb") as f:
    dd = pickle.load(f)

# remove exluded patients
data = data.iloc[(data.exclude_demo.values != 1) & (data.exclude_tum.values != 1), :]

## Fill in biomarker subtypes
data["biomarker_subtypes"] = None
for i in data.index:
    try:
        er = (
            dd["er_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "er_status"]
            ]
            == "Positive"
        )
    except:
        er = np.nan
    try:
        pr = (
            dd["pr_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "pr_status"]
            ]
            == "Positive"
        )
    except:
        pr = np.nan
    if not pd.isnull(data.loc[i, "her2_status"]):
        her2 = (
            dd["her2_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "her2_status"]
            ]
            == "Positive"
        )
    else:
        her2 = np.nan
    #
    if (er is True) | (pr is True):
        erpr = True
    else:
        erpr = False
    #
    if np.all(np.isnan([er, pr])) | np.isnan(her2):
        data.loc[i, "biomarker_subtypes"] = None
    elif (erpr) & (not her2):
        data.loc[i, "biomarker_subtypes"] = "ER/PR+ HER2-"
    elif (not erpr) & (not her2):
        data.loc[i, "biomarker_subtypes"] = "Triple Negative"
    elif her2:
        data.loc[i, "biomarker_subtypes"] = "HER2+"
    else:
        print("Unknown pattern at index:", i)

data["er_pr"] = None
for i in data.index:
    try:
        er = (
            dd["er_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "er_status"]
            ]
            == "Positive"
        )
    except:
        er = np.nan
    try:
        pr = (
            dd["pr_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "pr_status"]
            ]
            == "Positive"
        )
    except:
        pr = np.nan
    if (er is True) | (pr is True):
        data.loc[i, "er_pr"] = "ER/PR+"
    elif np.all(np.isnan([er, pr])):
        data.loc[i, "er_pr"] = None
    else:
        data.loc[i, "er_pr"] = "ER/PR-"

data["nodal_involvement_nac"] = None  # create new column for NAC summary purposes
nodalinv_map = {
    "N0": "N0",
    "N1": "N1-3",
    "N2": "N1-3",
    "N3": "N1-3",
    "NX": "NX",
    "pN0": "N0",
    "pN1": "N1-3",
    "pN1mi": "N1-3",
    "pN2": "N1-3",
    "pN3": "N1-3",
    "pNX": "NX",
    "ypN0": "N0",
    "ypN1": "N1-3",
    "ypN1mi": "N1-3",
    "ypN2": "N1-3",
    "ypN3": "N1-3",
    "ypNX": "NX",
}
for i in data.index:
    if pd.isnull(data.loc[i, "clin_node_stag_cat"]):
        if not pd.isnull(data.loc[i, "node_staging_category"]):
            raw_record = data.node_staging_category.map(
                dd["node_staging_category"]["Choices, Calculations, OR Slider Labels"]
            ).loc[i]
            data.loc[i, "nodal_involvement_nac"] = nodalinv_map[raw_record]
    else:
        raw_record = data.clin_node_stag_cat.map(
            dd["clin_node_stag_cat"]["Choices, Calculations, OR Slider Labels"]
        ).loc[i]
        data.loc[i, "nodal_involvement_nac"] = nodalinv_map[raw_record]

#####################
#### NAC summary ####
#####################

############################
## Biologic subtype table ##
############################

nac_biomarker_crosstab = pd.crosstab(data.biomarker_subtypes, data.nat)
nac_biomarker_crosstab.columns = nac_biomarker_crosstab.columns.map(
    {0: "No NAC", 1: "Had NAC"}
)
nac_biomarker_crosstab.columns.name = None
nac_biomarker_crosstab.index.name = None
nac_biomarker_crosstab_pct = nac_biomarker_crosstab.copy()
for ix in nac_biomarker_crosstab.index:
    for colname in nac_biomarker_crosstab.columns:
        ct = nac_biomarker_crosstab.loc[ix, colname]
        pct = (
            100
            * nac_biomarker_crosstab.loc[ix, colname]
            / nac_biomarker_crosstab.loc[ix, :].sum()
        )
        nac_biomarker_crosstab_pct.loc[ix, colname] = f"{ct} ({pct:.1f}%)"

nac_biomarker_crosstab_pct.to_csv(
    f"{dn}/summary_tables/nac_by_biomarker.csv", index=True, header=True
)

chi2, p, dof, expected = chi2_contingency(nac_biomarker_crosstab)

#############################
## Tumor size violin plots ##
#############################

fig, ax = plt.subplots(figsize=(8, 5))

nonac_patho = ax.violinplot(
    data.iloc[data.nat.values != 1, :]["tumor_size"].dropna().values,
    positions=[1],
    showmedians=True,
)
nac_clin = ax.violinplot(
    data.iloc[data.nat.values == 1, :]["clin_tumor_size"].dropna().values,
    positions=[2],
    showmedians=True,
)
nac_patho = ax.violinplot(
    data.iloc[data.nat.values == 1, :]["tumor_size"].dropna().values,
    positions=[3],
    showmedians=True,
)

ax.set_xticks([1, 2, 3])  # Position of the columns on the x-axis
ax.set_xticklabels(["No NAC\n(pathological)", "NAC\n(clinical)", "NAC\n(pathological)"])
ax.set_ylabel("Tumor size (cm)")
ax.set_title("Tumor size distributions by NAC category")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{dn}/plots/tumor_size_by_nac.png")
plt.close()

t, p = ttest_ind(
    data.iloc[data.nat.values != 1, :]["tumor_size"].dropna().values,
    data.iloc[data.nat.values == 1, :]["clin_tumor_size"].dropna().values,
)
print(
    (
        f"Tumor size comparison between no NAC "
        "(pathological size) and NAC (clinical size): t={t:.2f} (p={p:.2e})"
    )
)

t, p = ttest_ind(
    data.iloc[data.nat.values == 1, :]["clin_tumor_size"].dropna().values,
    data.iloc[data.nat.values == 1, :]["tumor_size"].dropna().values,
)
print(
    (
        f"Tumor size comparison between NAC "
        f"(clinical size) and NAC (pathological size): t={t:.2f} (p={p:.2e})"
    )
)

#############################
## Nodal involvement table ##
#############################

# Distribution among NAC patients
nac_clinical = (
    data.iloc[data.nat.values == 1, :]
    .clin_node_stag_cat.map(
        dd["clin_node_stag_cat"]["Choices, Calculations, OR Slider Labels"]
    )
    .value_counts()
)
nac_patho = (
    data.iloc[data.nat.values == 1, :]
    .node_staging_category.map(
        dd["node_staging_category"]["Choices, Calculations, OR Slider Labels"]
    )
    .value_counts()
)
# Distribution among all patients
nonac_patho = (
    data.iloc[data.nat.values != 1, :]
    .node_staging_category.map(
        dd["node_staging_category"]["Choices, Calculations, OR Slider Labels"]
    )
    .value_counts()
)

nodal_involvement = pd.DataFrame(
    index=["N0", "N1-3", "NX"],
    columns=["No NAC (patho)", "Had NAC (clinical)", "Had NAC (patho)"],
)

for label, df in zip(
    ["No NAC (patho)", "Had NAC (clinical)", "Had NAC (patho)"],
    [nonac_patho, nac_clinical, nac_patho],
):
    nodal_involvement.loc["N0", label] = df.loc[["N0" in x for x in df.index]].sum()
    nodal_involvement.loc["N1-3", label] = df.iloc[
        [("N1" in x) | ("N2" in x) | ("N3" in x) for x in df.index]
    ].sum()
    nodal_involvement.loc["NX", label] = df.iloc[["NX" in x for x in df.index]].sum()

nodal_involvement_pct = nodal_involvement.copy()
for i in nodal_involvement.index:
    for j in nodal_involvement.columns:
        ct = nodal_involvement.loc[i, j]
        pct = 100 * ct / nodal_involvement.loc[:, j].sum(axis=0)
        nodal_involvement_pct.loc[i, j] = f"{ct} ({pct:.1f}%)"

for label in ["No NAC (patho)", "Had NAC (clinical)", "Had NAC (patho)"]:
    total_ct = nodal_involvement.loc[:, label].sum(axis=0)
    nodal_involvement_pct.loc["Total", label] = f"{total_ct} (100.0%)"

nodal_involvement_pct.to_csv(f"{dn}/summary_tables/nodal_involvement_nac_vs_not.csv")

#####################################################################################
## Nodal involvement separated by NAC status and biologic subtype proportion chart ##
#####################################################################################

nodal_groups = data.groupby(
    ["biomarker_subtypes", "nat", "nodal_involvement_nac"]
).size()

total_by_subtype = {
    subtype: nodal_groups.loc[subtype].values.sum()
    for subtype in ["ER/PR+ HER2-", "HER2+", "Triple Negative"]
}

cts = {
    nac_status: [
        nodal_groups.loc[subtype, i, node_status]
        for subtype in ["ER/PR+ HER2-", "HER2+", "Triple Negative"]
        for node_status in ["N0", "N1-3", "NX"]
    ]
    for i, nac_status in enumerate(["No NAC", "Had NAC"])
}

pcts = {}
for key in cts.keys():
    pcts[key] = []
    for i in range(len(cts[key])):
        pcts[key].append(100 * cts[key][i] / np.sum([cts[x][i] for x in cts.keys()]))

fig, ax = plt.subplots(figsize=(8, 4))
pos = [1, 2, 3, 5, 6, 7, 9, 10, 11]

ax.bar(
    pos,
    pcts["No NAC"],
    label="No NAC",
    color="royalblue",
    edgecolor="black",
    linewidth=0.5,
)
ax.bar(
    pos,
    pcts["Had NAC"],
    bottom=pcts["No NAC"],
    label="Had NAC",
    color="gold",
    edgecolor="black",
    linewidth=0.5,
)

ax.legend(loc="lower left")
ax.set_ylim(0, 100)
ax.set_ylabel("Percentages")

fig.subplots_adjust(bottom=0.20)

ax.set_xticks(pos)
fig.text(
    0.25,
    0.03,
    "ER/PR+ HER2-\n(n={})".format(total_by_subtype["ER/PR+ HER2-"]),
    ha="center",
    color="black",
    fontsize=12,
)
fig.text(
    0.51,
    0.03,
    "HER2+\n(n={})".format(total_by_subtype["HER2+"]),
    ha="center",
    color="black",
    fontsize=12,
)
fig.text(
    0.76,
    0.03,
    "Triple Negative\n(n={})".format(total_by_subtype["Triple Negative"]),
    ha="center",
    color="black",
    fontsize=12,
)

ax.set_xticklabels(["N0", "N1-3", "NX", "N0", "N1-3", "NX", "N0", "N1-3", "NX"])
ax.set_title("Proportion of patients who received NAC by nodal involvement")
fig.savefig(f"{dn}/plots/nac_proportion_by_node_biomarker.png")
plt.close()


###############################
## Diagnosis age violin plot ##
###############################

fig, ax = plt.subplots(figsize=(6, 4))

nonac_age = ax.violinplot(
    data.iloc[data.nat.values != 1, :]["age_at_diagnosis"].dropna().values,
    positions=[1],
    showmedians=True,
)
nac_age = ax.violinplot(
    data.iloc[data.nat.values == 1, :]["age_at_diagnosis"].dropna().values,
    positions=[2],
    showmedians=True,
)

# Customize the violin plot
ax.set_xticks([1, 2])  # Position of the columns on the x-axis
ax.set_xticklabels(["No NAC", "NAC"])
ax.set_ylabel("Age at diagnosis")
ax.set_title("Age distributions by NAC category")

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{dn}/plots/age_by_nac.png")
plt.close()

t, p = ttest_ind(
    data.iloc[data.nat.values != 1, :]["age_at_diagnosis"].dropna().values,
    data.iloc[data.nat.values == 1, :]["age_at_diagnosis"].dropna().values,
)
print(f"Age comparison: t={t:.2f} (p={p:.2e})")


###################################################################
## Age separated by NAC status and biologic subtype violin plots ##
###################################################################

fig, ax = plt.subplots(figsize=(6, 4))

pos = [1, 2, 4, 5, 7, 8]
violins = ax.violinplot(
    [
        data.iloc[
            (data.nat.values != 1) & (data.biomarker_subtypes.values == "ER/PR+ HER2-"),
            :,
        ]["age_at_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1) & (data.biomarker_subtypes.values == "ER/PR+ HER2-"),
            :,
        ]["age_at_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values != 1) & (data.biomarker_subtypes.values == "HER2+"), :
        ]["age_at_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1) & (data.biomarker_subtypes.values == "HER2+"), :
        ]["age_at_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values != 1)
            & (data.biomarker_subtypes.values == "Triple Negative"),
            :,
        ]["age_at_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1)
            & (data.biomarker_subtypes.values == "Triple Negative"),
            :,
        ]["age_at_diagnosis"]
        .dropna()
        .values,
    ],
    pos,
    widths=0.7,
    showmedians=True,
    showextrema=True,
)

for i in [0, 2, 4]:
    violins["bodies"][i].set_facecolor("blue")

for i in [1, 3, 5]:
    violins["bodies"][i].set_facecolor("orange")

ax.set_xticks([1.5, 4.5, 7.5])
ax.set_xticklabels(
    ["ER/PR+ HER2-", "HER2+", "Triple Negative"], ha="right", rotation=30
)
ax.set_ylabel("Age at diagnosis")
ax.set_title("Diagnosis age by NAC status and biologic subtype")
ax.legend(["No NAC", "Had NAC"], loc="lower right")

plt.tight_layout()
fig.savefig(f"{dn}/plots/age_by_nac_and_biomarker.png")
plt.close()

#######################################################################
## Year Dx separated by NAC status and biologic subtype violin plots ##
#######################################################################

fig, ax = plt.subplots(figsize=(6, 4))

data = data.iloc[data.year_of_diagnosis.values >= 2010, :]

pos = [1, 2, 4, 5, 7, 8]
violins = ax.violinplot(
    [
        data.iloc[
            (data.nat.values != 1) & (data.biomarker_subtypes.values == "ER/PR+ HER2-"),
            :,
        ]["year_of_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1) & (data.biomarker_subtypes.values == "ER/PR+ HER2-"),
            :,
        ]["year_of_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values != 1) & (data.biomarker_subtypes.values == "HER2+"), :
        ]["year_of_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1) & (data.biomarker_subtypes.values == "HER2+"), :
        ]["year_of_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values != 1)
            & (data.biomarker_subtypes.values == "Triple Negative"),
            :,
        ]["year_of_diagnosis"]
        .dropna()
        .values,
        data.iloc[
            (data.nat.values == 1)
            & (data.biomarker_subtypes.values == "Triple Negative"),
            :,
        ]["year_of_diagnosis"]
        .dropna()
        .values,
    ],
    pos,
    widths=0.7,
    showmedians=True,
    showextrema=True,
)

for i in [0, 2, 4]:
    violins["bodies"][i].set_facecolor("blue")

for i in [1, 3, 5]:
    violins["bodies"][i].set_facecolor("orange")

ax.set_xticks([1.5, 4.5, 7.5])
ax.set_xticklabels(
    ["ER/PR+ HER2-", "HER2+", "Triple Negative"], ha="right", rotation=30
)
ax.set_ylabel("Year of diagnosis")
ax.set_title("Diagnosis year by NAC status and biologic subtype")
ax.legend(["No NAC", "Had NAC"], loc="lower right")

plt.tight_layout()
fig.savefig(f"{dn}/plots/year_dx_by_nac_and_biomarker.png")
plt.close()


#####################################
## Year Dx as a population pyramid ##
#####################################

cts = {
    "ER/PR+ HER2-": {"Had NAC": [], "No NAC": []},
    "HER2+": {"Had NAC": [], "No NAC": []},
    "Triple Negative": {"Had NAC": [], "No NAC": []},
}

years = np.arange(2010, 2021)
for dx_year in years:
    for biologic_type in list(cts):
        no_nac = data.iloc[
            (data.nat.values != 1)
            & (data.biomarker_subtypes.values == biologic_type)
            & (data.year_of_diagnosis.values == dx_year),
            :,
        ].shape[0]
        had_nac = data.iloc[
            (data.nat.values == 1)
            & (data.biomarker_subtypes.values == biologic_type)
            & (data.year_of_diagnosis.values == dx_year),
            :,
        ].shape[0]
        cts[biologic_type]["Had NAC"].append(had_nac)
        cts[biologic_type]["No NAC"].append(no_nac)

fig, ax = plt.subplots(1, 3, figsize=(6, 4))

for i, biologic_type in enumerate(list(cts)):
    max_num = np.max(cts[biologic_type]["Had NAC"] + cts[biologic_type]["No NAC"])
    ax[i].barh(
        np.arange(len(years)), cts[biologic_type]["Had NAC"], color="gold", height=1.0
    )
    ax[i].barh(
        np.arange(len(years)),
        -1 * np.array(cts[biologic_type]["No NAC"]),
        color="royalblue",
        height=1.0,
    )
    ax[i].set_xlim(-1.05 * max_num, 1.05 * max_num)
    xtick_num = int(1.05 * max_num - ((1.05 * max_num) % 10))
    ax[i].set_xticks([-1 * xtick_num, 0, xtick_num])
    ax[i].set_xticklabels([xtick_num, 0, xtick_num])
    ax[i].spines["top"].set_visible(False)
    ax[i].spines["right"].set_visible(False)
    ax[i].set_xlabel(biologic_type)
    if i == 0:
        ax[i].set_yticks([0, 5, 10])
        ax[i].set_yticklabels([2010, 2015, 2020])
        ax[i].legend(["Had NAC", "No NAC"], loc="lower left")
    else:
        ax[i].spines["left"].set_visible(False)
        ax[i].set_yticks([])

fig.suptitle("Number of patients who received NAC by diagnosis year")

plt.tight_layout()
fig.savefig(f"{dn}/plots/year_dx_by_nac_and_biomarker_horizontalbars.png")
plt.close()
