# pylint: disable=C0114
# pylint: disable=W0621
# pylint: disable=R0914

import pickle

import pandas as pd

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
data.parous = data.parous.fillna(0.0)

#############################################
#### Select necessary columns and export ####
#############################################

# Define gene names
genes = [
    # "any_patho_mutation",
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

# Map to binary encoding
for feature_name in genes:
    data[feature_name] = (
        data[feature_name].map(
            dd[feature_name]["Choices, Calculations, OR Slider Labels"]
        )
        == "Pathogenic"
    ).astype(int)

# select necessary columns
subdata_parous = data.iloc[data.parous.values == 1.0, :][
    [
        "record_id",
        "age_at_diagnosis",
        "year_of_diagnosis",
        "age_at_most_recent_pregnancy",
        "years_since_pregnancy",
        "fam_hx",
    ]
    + genes
]

subdata_all = data[
    [
        "record_id",
        "parous",
        "age_at_diagnosis",
        "year_of_diagnosis",
        "age_at_most_recent_pregnancy",
        "years_since_pregnancy",
        "fam_hx",
    ]
    + genes
]
subdata_parous.to_csv("/home/srd6051/recent_parity_data_parous.csv", index=False)
subdata_all.to_csv("/home/srd6051/recent_parity_data_all.csv", index=False)
