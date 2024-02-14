# pylint: disable=C0114
# pylint: disable=W0621
# pylint: disable=R0914
# pylint: disable=bare-except
# pylint: disable=duplicate-code

import pickle
from datetime import datetime

import numpy as np
import pandas as pd
import recent_parity_stat
from scipy.stats import ttest_ind

dn = "/share/fsmresfiles/breast_cancer_pregnancy"
dout = "/share/fsmresfiles/breast_cancer_pregnancy/stat_results"

###################
#### Read data ####
###################

datadir = "data/06_exported_from_redcap"
fn = "FrequencyAndResultsO_DATA_2023-03-17_0949.csv"

data = pd.read_csv(f"{dn}/{datadir}/{fn}")

# remove exluded patients
data = data.iloc[(data.exclude_demo.values != 1) & (data.exclude_tum.values != 1), :]

data = data.iloc[data.nat.values == 1,]

with open(f"{dn}/{datadir}/data_dictionary.p", "rb") as f:
    dd = pickle.load(f)

## Fill in parity category
data = recent_parity_stat.fill_in_parity_categorical(data)

## Fill in biomarker subtypes
data["biomarker_subtypes"] = None
for i in data.index:
    # ER/PR status
    try:
        er = (
            dd["er_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "er_status"]
            ]
            == "Positive"
        )
        pr = (
            dd["pr_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "pr_status"]
            ]
            == "Positive"
        )
    except:
        er = False
        pr = False
        print("Missing ER/PR at index:", i)
    # HER2 status
    try:
        her2 = (
            dd["her2_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "her2_status"]
            ]
            == "Positive"
        )
    except:
        her2 = False
        print("Missing HER2 at index:", i)
    # fill in subtypes
    if (er | pr) & (not her2):
        data.loc[i, "biomarker_subtypes"] = "ER/PR+ HER2-"
    elif (not er) & (not pr) & (not her2):
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
        pr = (
            dd["pr_status"]["Choices, Calculations, OR Slider Labels"][
                data.loc[i, "pr_status"]
            ]
            == "Positive"
        )
    except:
        er = False
        pr = False
        print("Missing ER/PR at index:", i)
    if er | pr:
        data.loc[i, "er_pr"] = "ER/PR+"
    else:
        data.loc[i, "er_pr"] = "ER/PR-"

data["pCR"] = (
    data.response_nat.map(dd["response_nat"]["Choices, Calculations, OR Slider Labels"])
    == "No residual disease"
).astype(int)

data["positive_response"] = (
    (
        (
            data.rcb_category.map(
                dd["rcb_category"]["Choices, Calculations, OR Slider Labels"]
            )
            == "0"
        ).values
        | (
            data.rcb_category.map(
                dd["rcb_category"]["Choices, Calculations, OR Slider Labels"]
            )
            == "1"
        ).values
    )
).astype(int)

data["poor_response"] = (
    (
        (
            data.rcb_category.map(
                dd["rcb_category"]["Choices, Calculations, OR Slider Labels"]
            )
            == "2"
        ).values
        | (
            data.rcb_category.map(
                dd["rcb_category"]["Choices, Calculations, OR Slider Labels"]
            )
            == "3"
        ).values
    )
).astype(int)

##########################
#### Perform Analysis ####
##########################

datestring = datetime.now().strftime("%Y%m%d")

parity_comparisons = {
    # 'Parous vs. Nulliparous' : {'ref' : 'nulliparous', 'comp' : 'parous'},
    # '<5 vs. >=5 years' : {'ref' : 'parity >=5 years', 'comp' : 'parity <5 years'},
    # '<10 vs. >=10 years' : {'ref' : 'parity >=10 years', 'comp' : 'parity <10 years'},
    # '<5 vs. >=10 years' : {'ref' : 'parity >=10 years', 'comp' : 'parity <5 years'},
    "<5 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity <5 years"},
    "5-10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity 5-10 years"},
    "<10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity <10 years"},
    ">=5 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity >=5 years"},
    ">=10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity >=10 years"},
}

stratifications = data.biomarker_subtypes.unique()
feature_names = ["poor_response", "positive_response", "pCR"]

for parity_comparison, _ in parity_comparisons.items():
    print(f"#### {parity_comparison} ####")
    ref = parity_comparisons[parity_comparison]["ref"]
    comp = parity_comparisons[parity_comparison]["comp"]
    results = pd.DataFrame(
        columns=[
            "Variable of interest",
            "Stratification",
            "Category",
            comp,
            ref,
            "OR (95% CI)",
            "p-value",
        ]
    )
    for feature_name in feature_names:
        for stratification in stratifications:
            subdata = data.iloc[data.biomarker_subtypes.values == stratification, :]
            X, y = recent_parity_stat.generate_lrdata(
                subdata, parity_ref=ref, parity_comp=comp, feature_name=feature_name
            )
            if (np.sum((X[comp] == 0) & (y == 1)) == 0) | (
                np.sum((X[comp] == 1) & (y == 1)) == 0
            ):
                print(
                    f"# Cannot perform comparison for {feature_name} due to lack of data\n"
                )
            else:
                crosstab = pd.crosstab(X[comp], y)
                totals = crosstab.sum(axis=0)
                for i in crosstab.index:
                    for j in crosstab.columns:
                        cts = crosstab.loc[i, j]
                        pcts = 100 * cts / totals[j]
                        crosstab.loc[i, j] = f"{cts} ({pcts:.1f}%)"
                oddsratios, cis, pvals = recent_parity_stat.get_oddsratio_ci(X, y)
                or_formatted = f"{oddsratios[0]:.2f} ({cis[0][0]:.2f}-{cis[0][1]:.2f})"
                pval_formatted = f"{pvals[0]:.4f}"
                results = pd.concat(
                    (
                        results,
                        pd.DataFrame(
                            {
                                "Variable of interest": [feature_name, None],
                                "Stratification": [stratification, None],
                                "Category": ["Yes", "No"],
                                parity_comparisons[parity_comparison]["comp"]: [
                                    crosstab.loc[1, 1],
                                    crosstab.loc[1, 0],
                                ],
                                parity_comparisons[parity_comparison]["ref"]: [
                                    crosstab.loc[0, 1],
                                    crosstab.loc[0, 0],
                                ],
                                "OR (95% CI)": [or_formatted, None],
                                "p-value": [pval_formatted, None],
                            }
                        ),
                    )
                )
    fout = datestring + "_nat_" + parity_comparison.replace(" ", "_") + ".csv"
    results.to_csv(f"{dout}/{fout}", index=False)

#############################
#### t-test of RCB score ####
#############################

print("\n## Parous vs. Nulliparous Comparison ##")
# Non-stratified
t, p = ttest_ind(
    data.rcb.iloc[data.parous.values == 1], data.rcb.iloc[data.parous.values == 0]
)
print(f"Overall: t={t:.2f} (p={p:.4f})")

# Stratified
for bm_cat in ["ER/PR+ HER2-", "Triple Negative", "HER2+"]:
    subdata = data.iloc[data.biomarker_subtypes.values == bm_cat, :]
    t, p = ttest_ind(
        subdata.rcb.iloc[subdata.parous.values == 1],
        subdata.rcb.iloc[subdata.parous.values == 0],
    )
    print(f"{bm_cat}: t={t:.2f} (p={p:.4f})")

print("\n\n## Recent Parity vs. Non-Recent Parity Comparison ##")
for thres in [5, 10]:
    print(f"\n# {thres}-year recency threshold")
    t, p = ttest_ind(
        data.rcb.iloc[data.years_since_pregnancy.values < thres],
        data.rcb.iloc[data.years_since_pregnancy.values >= thres],
    )
    print(f"Overall: t={t:.2f} (p={p:.4f})")
    # Stratified
    for bm_cat in ["ER/PR+ HER2-", "Triple Negative", "HER2+"]:
        subdata = data.iloc[data.biomarker_subtypes.values == bm_cat, :]
        t, p = ttest_ind(
            subdata.rcb.iloc[subdata.years_since_pregnancy.values < thres],
            subdata.rcb.iloc[subdata.years_since_pregnancy.values >= thres],
        )
        print(f"{bm_cat}: t={t:.2f} (p={p:.4f})")
