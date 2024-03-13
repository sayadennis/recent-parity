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

#############################################
#### Select necessary columns and export ####
#############################################

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

subdata = data.iloc[data.parous.values == 1.0, :][
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

subdata.to_csv("/home/srd6051/recent_parity_quest.csv", index=False)
