# pylint: disable=C0114
# pylint: disable=W0621
# pylint: disable=R0914

import pickle
from datetime import datetime

import numpy as np
import pandas as pd
import recent_parity_stat

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

with open(f"{dn}/{datadir}/data_dictionary.p", "rb") as f:
    dd = pickle.load(f)

# mark family history as none for patients who are missing it
data.fam_hx = data.fam_hx.fillna(0.0)

#########################################
#### Fill in useful categorical info ####
#########################################

data = recent_parity_stat.fill_in_parity_categorical(data)

genes = [
    "any_patho_mutation",
    "brca1",
    "brca2",
    "palb2",
    "chek2",
    "atm",
    "tp53",
    "pten",
    "stk11",
    "cdh1",
]

data["any_patho_mutation"] = (
    data["any_patho_mutation"].map(
        dd["any_patho_mutation"]["Choices, Calculations, OR Slider Labels"]
    )
    == "Pathogenic"
).astype(int)

for feature_name in genes[1:]:
    data[feature_name] = (
        data[feature_name].map(
            dd[feature_name]["Choices, Calculations, OR Slider Labels"]
        )
        == "Pathogenic"
    ).astype(int)

genes = [item for item in genes if data[item].sum(axis=0) >= 2]

##########################
#### Perform Analysis ####
##########################

datestring = datetime.now().strftime("%Y%m%d")

parity_comparisons = {
    "Parous vs. Nulliparous": {"ref": "nulliparous", "comp": "parous"},
    "<5 vs. >=5 years": {"ref": "parity >=5 years", "comp": "parity <5 years"},
    "<10 vs. >=10 years": {"ref": "parity >=10 years", "comp": "parity <10 years"},
    "<5 vs. >=10 years": {"ref": "parity >=10 years", "comp": "parity <5 years"},
    "<5 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity <5 years"},
    "5-10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity 5-10 years"},
    "<10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity <10 years"},
    ">=5 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity >=5 years"},
    ">=10 years vs. Nulliparous": {"ref": "nulliparous", "comp": "parity >=10 years"},
}

for parity_comparison, _ in parity_comparisons.items():
    print(f"#### {parity_comparison} ####")
    ref = parity_comparisons[parity_comparison]["ref"]
    comp = parity_comparisons[parity_comparison]["comp"]
    results = pd.DataFrame(
        columns=[
            "Variable of interest",
            "Category",
            comp,
            ref,
            "OR (95% CI)",
            "p-value",
        ]
    )
    for feature_name in genes:
        X, y = recent_parity_stat.generate_lrdata(
            data, parity_ref=ref, parity_comp=comp, feature_name=feature_name
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
    fout = datestring + "_genetic_" + parity_comparison.replace(" ", "_") + ".csv"
    results.to_csv(f"{dout}/{fout}", index=False)
