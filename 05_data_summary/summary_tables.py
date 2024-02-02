# pylint: disable=missing-module-docstring
# pylint: disable=bare-except
# pylint: disable=duplicate-code

import pickle
from datetime import datetime

import numpy as np
import pandas as pd

dn = "/share/fsmresfiles/breast_cancer_pregnancy"

## Load data & data dictionary
datadir = "data/06_exported_from_redcap"
fn = "FrequencyAndResultsO_DATA_2023-03-17_0949.csv"
data = pd.read_csv(f"{dn}/{datadir}/{fn}")

with open(f"{dn}/{datadir}/data_dictionary.p", "rb") as f:
    dd = pickle.load(f)

# remove exluded patients
data = data.iloc[(data.exclude_demo.values != 1) & (data.exclude_tum.values != 1), :]

# Select only patients with neoadjuvant therapy if necessary
sumtabdir = "summary_tables"

#########################################
#### Fill in useful categorical info ####
#########################################

## Fill in parity category
data["parity_category"] = None
for i in data.index:
    if data.loc[i, "parous"] == 0:
        data.loc[i, "parity_category"] = "Nulliparous"
    elif data.loc[i, "years_since_pregnancy"] < 5:
        data.loc[i, "parity_category"] = "<5 years"
    elif data.loc[i, "years_since_pregnancy"] < 10:
        data.loc[i, "parity_category"] = "5-10 years"
    elif data.loc[i, "years_since_pregnancy"] >= 10:
        data.loc[i, "parity_category"] = ">=10 years"
    else:
        print(f"Years since pregnancy not available for index: {i}")

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
    if (er is True) | (pr is True):  # pylint: disable=simplifiable-if-statement
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

super_table = pd.DataFrame(
    columns=[
        "super_category",
        "main_category",
        "sub_category",
        "Nulliparous",
        "<5 years",
        "5-10 years",
        ">=10 years",
        "Overall",
    ]
)

#####################
#### Demographic ####
#####################

#### Age ####
main_table = pd.DataFrame(index=["Mean age (±SD)"], columns=super_table.columns)
mean = data.age_at_diagnosis.mean()
std = data.age_at_diagnosis.std()
main_table.loc["Mean age (±SD)", "Overall"] = f"{mean:.2f} (±{std:.2f})"

for par_cat in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
    mean = data.iloc[data.parity_category.values == par_cat, :].age_at_diagnosis.mean()
    std = data.iloc[data.parity_category.values == par_cat, :].age_at_diagnosis.std()
    main_table.loc["Mean age (±SD)", par_cat] = f"{mean:.2f} (±{std:.2f})"

main_table["super_category"] = "demo"
main_table["main_category"] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((super_table, main_table))

###############################
#### Gynecological history ####
###############################

main_table = pd.DataFrame(
    index=["Mean Gravida (±SD)", "Mean Para (±SD)"],
    columns=["Nulliparous", "<5 years", "5-10 years", ">=10 years"],
)

for varname, tablevarname in zip(
    ["number_pregnancies", "number_births"], ["Mean Gravida (±SD)", "Mean Para (±SD)"]
):
    mean = data[varname].mean()
    std = data[varname].std()
    main_table.loc[tablevarname, "Overall"] = f"{mean:.2f} (±{std:.2f})"


for par_cat in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
    # Gravida
    mean = data.iloc[
        data.parity_category.values == par_cat, :
    ].number_pregnancies.mean()
    std = data.iloc[data.parity_category.values == par_cat, :].number_pregnancies.std()
    main_table.loc["Mean Gravida (±SD)", par_cat] = f"{mean:.2f} (±{std:.2f})"
    # Para
    mean = data.iloc[data.parity_category.values == par_cat, :].number_births.mean()
    std = data.iloc[data.parity_category.values == par_cat, :].number_births.std()
    main_table.loc["Mean Para (±SD)", par_cat] = f"{mean:.2f} (±{std:.2f})"

main_table["super_category"] = "gynecological"
main_table["main_category"] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((super_table, main_table))

##############################
#### Genetic test results ####
##############################

# Common mutations
main_table = pd.DataFrame(
    columns=["Nulliparous", "<5 years", "5-10 years", ">=10 years"]
)

for key in dd.keys():
    if dd[key]["Form Name"] == "geneticsfam_hx":
        if dd[key]["Field Type"] == "radio":
            if key in ["any_patho_mutation", "any_vus_mutation"]:
                cts = {}
                for par_cat in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
                    cts[par_cat] = np.sum(
                        data.iloc[data.parity_category.values == par_cat, :][key].map(
                            dd[key]["Choices, Calculations, OR Slider Labels"]
                        )
                        == "Present"
                    )
            else:
                cts = {}
                for par_cat in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
                    cts[par_cat] = np.sum(
                        data.iloc[data.parity_category.values == par_cat, :][key].map(
                            dd[key]["Choices, Calculations, OR Slider Labels"]
                        )
                        == "Pathogenic"
                    )
            main_table = pd.concat((main_table, pd.DataFrame(cts, index=[key.upper()])))

main_table["Overall"] = (
    main_table[["Nulliparous", "<5 years", "5-10 years", ">=10 years"]]
    .sum(axis=1)
    .astype(int)
)

for j in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
    ol = (
        main_table.iloc[np.invert([x.startswith("ANY") for x in main_table.index]), :][
            j
        ].sum()
        - main_table.loc["ANY_PATHO_MUTATION", j]
    )
    print(f"{j}: {ol} patients belong in more than one category")

for i in main_table.index:
    for j in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
        percentage = 100 * main_table.loc[i, j] / (data.parity_category == j).sum()
        main_table.loc[i, j] = f"{main_table.loc[i,j]} ({percentage:.1f}%)"

main_table["super_category"] = "genetic"
main_table["main_category"] = main_table.index
main_table.reset_index(drop=True, inplace=True)

super_table = pd.concat((super_table, main_table))

###############################
#### Tumor characteristics ####
###############################

main_table = pd.DataFrame(
    columns=[
        "main_category",
        "sub_category",
        "Nulliparous",
        "<5 years",
        "5-10 years",
        ">=10 years",
        "Overall",
    ]
)

for colname in [
    "biomarker_subtypes",
    "histology",
    "histologic_grade",
    "tumor_staging_category",
    "node_staging_category",
]:  # 'clin_tumor_stag_cat', 'clin_node_stag_cat'
    if colname == "histologic_grade":
        grade_converted_value = None  # initialize variable
        for i in data.index:
            is_dcis = (
                dd["histology"]["Choices, Calculations, OR Slider Labels"][
                    int(data.loc[i, "histology"])
                ]
                == "DCIS"
                if not pd.isnull(data.loc[i, "histology"])
                else False
            )
            is_missing = pd.isnull(data.loc[i, "histologic_grade"])
            if not is_missing:
                grade_converted_value = dd[colname][
                    "Choices, Calculations, OR Slider Labels"
                ][int(data.loc[i, "histologic_grade"])]
            else:
                grade_converted_value = None
            #
            if is_dcis:
                data.loc[i, "histologic_grade"] = "No histologic grade (DCIS)"
            elif is_missing | (grade_converted_value == "No histologic grade (DCIS)"):
                data.loc[i, "histologic_grade"] = "missing"
            else:
                data.loc[i, "histologic_grade"] = grade_converted_value
        #
        subtable = pd.crosstab(data[colname], data["parity_category"]).reset_index()
    elif colname != "biomarker_subtypes":
        subtable = (
            pd.crosstab(data[colname], data["parity_category"])
            .rename(dd[colname]["Choices, Calculations, OR Slider Labels"], axis=0)
            .reset_index()
        )
        subtable = pd.concat(
            (
                subtable,
                pd.DataFrame(
                    {colname: "missing"}
                    | pd.crosstab(
                        pd.isnull(data[colname]).astype(int), data["parity_category"]
                    )
                    .loc[1, :]
                    .to_dict(),
                    index=[0],
                ),
            )
        ).reset_index(drop=True)
    else:
        subtable = pd.crosstab(data[colname], data["parity_category"]).reset_index()
        subtable = pd.concat(
            (
                subtable,
                pd.DataFrame(
                    {colname: "missing"}
                    | pd.crosstab(
                        pd.isnull(data[colname]).astype(int), data["parity_category"]
                    )
                    .loc[1, :]
                    .to_dict(),
                    index=[0],
                ),
            )
        ).reset_index(drop=True)
    #
    subtable = subtable.rename({colname: "sub_category"}, axis=1)
    subtable.columns.name = None
    subtable["main_category"] = colname
    main_table = pd.concat((main_table, subtable))

main_table["Overall"] = (
    main_table[["Nulliparous", "<5 years", "5-10 years", ">=10 years"]]
    .sum(axis=1)
    .astype(int)
)

main_table.reset_index(drop=True, inplace=True)

for tumvarcat in main_table.main_category.unique():
    subtable = main_table.iloc[main_table.main_category.values == tumvarcat, :]
    totals = subtable.sum(axis=0)
    cts = subtable[["Nulliparous", "<5 years", "5-10 years", ">=10 years", "Overall"]]
    cts.index = subtable.sub_category
    pcts = (
        subtable[["Nulliparous", "<5 years", "5-10 years", ">=10 years", "Overall"]]
        / totals[["Nulliparous", "<5 years", "5-10 years", ">=10 years", "Overall"]]
    )
    pcts.index = subtable.sub_category
    for subcat in subtable.sub_category.values:
        for par_cat in [
            "Nulliparous",
            "<5 years",
            "5-10 years",
            ">=10 years",
            "Overall",
        ]:
            ct = cts.loc[subcat, par_cat]
            pct = 100 * pcts.loc[subcat, par_cat]
            main_table.iloc[
                (
                    (main_table.main_category.values == tumvarcat)
                    & (main_table.sub_category.values == subcat)
                ),
                main_table.columns == par_cat,
            ] = f"{ct} ({pct:.1f}%)"

main_table = pd.concat(
    (main_table, pd.DataFrame({"main_category": "tumor_size"}, index=["tumor_size"]))
)
mean = data.tumor_size.mean()
std = data.tumor_size.std()
main_table.loc["tumor_size", "Overall"] = f"{mean:.2f} (±{std:.2f})"
for par_cat in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
    mean = data.iloc[data.parity_category.values == par_cat, :].tumor_size.mean()
    std = data.iloc[data.parity_category.values == par_cat, :].tumor_size.std()
    main_table.loc["tumor_size", par_cat] = f"{mean:.2f} (±{std:.2f})"

main_table["super_category"] = "tumor_char"

super_table = pd.concat((super_table, main_table))

###################
#### Treatment ####
###################

main_table = pd.DataFrame(
    columns=[
        "main_category",
        "sub_category",
        "Nulliparous",
        "<5 years",
        "5-10 years",
        ">=10 years",
    ]
)

for colname in ["nat", "nat_reg", "her2_therapy", "rcb", "rcb_category"]:
    if colname == "nat":
        sub_table = pd.crosstab(data[colname], data["parity_category"]).rename(
            {1: "Yes", 0: "No"}, axis=0
        )
    elif colname == "rcb":
        continue
        ## The below needs to NOT be part of this loop. Fix ##
        # sub_table = pd.DataFrame(
        #     index=['Mean RCB (±SD)'],
        #     columns=['Nulliparous', '<5 years', '5-10 years', '>=10 years']
        # )
        # for par_cat in sub_table.columns:
        #     mean = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.mean()
        #     std = data.iloc[data.parity_category.values==par_cat,:].number_pregnancies.std()
        # sub_table.loc['Mean RCB (±SD)',par_cat] = f'{mean} (±{std})'
    else:
        sub_table = pd.crosstab(data[colname], data["parity_category"]).rename(
            dd[colname]["Choices, Calculations, OR Slider Labels"], axis=0
        )
    sub_table = sub_table.reset_index()
    sub_table.columns.name = None
    sub_table = sub_table.rename({colname: "sub_category"}, axis=1)
    sub_table["main_category"] = colname
    main_table = pd.concat((main_table, sub_table))


main_table["Overall"] = (
    main_table[["Nulliparous", "<5 years", "5-10 years", ">=10 years"]]
    .sum(axis=1)
    .astype(int)
)

main_table["super_category"] = "treatment"
main_table.reset_index(drop=True, inplace=True)

for natvarcat in main_table.main_category.unique():
    subtable = main_table.iloc[main_table.main_category.values == natvarcat, :]
    totals = subtable.sum(axis=0)
    cts = subtable[["Nulliparous", "<5 years", "5-10 years", ">=10 years", "Overall"]]
    cts.index = subtable.sub_category
    pcts = cts.div(cts["Overall"], axis=0)
    pcts.index = subtable.sub_category
    for subcat in subtable.sub_category.values:
        for par_cat in [
            "Nulliparous",
            "<5 years",
            "5-10 years",
            ">=10 years",
            "Overall",
        ]:
            ct = cts.loc[subcat, par_cat]
            pct = 100 * pcts.loc[subcat, par_cat]
            main_table.iloc[
                (
                    (main_table.main_category.values == natvarcat)
                    & (main_table.sub_category.values == subcat)
                ),
                main_table.columns == par_cat,
            ] = f"{ct} ({pct:.1f}%)"


super_table = pd.concat((super_table, main_table))

datestring = datetime.now().strftime("%Y%m%d")

super_table.to_csv(
    f"{dn}/{sumtabdir}/{datestring}_super_summary_table.csv", index=False
)
