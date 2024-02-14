"""
Functions to facilitate the statistical analysis.
"""

# pylint: disable=too-many-locals

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.utils import resample

PARITY_CATEGORIES = [
    "nulliparous",
    "parity <5 years",
    "parity 5-10 years",
    "parity >=10 years",
]


def fill_in_parity_categorical(data: pd.DataFrame) -> pd.DataFrame:
    """
    Add new columns related to parity categories to facilitate statistical analysis.

    Parameters:
        - data: REDCap table

    Returns:
        - The updated data frame
    """
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

    data["nulliparous"] = np.invert(data["parous"].astype(bool)).astype(float)

    data["parity <5 years"] = np.where(
        (data["parous"] == 0) | np.isnan(data["years_since_pregnancy"]),
        np.nan,
        data["years_since_pregnancy"] < 5,
    )

    data["parity >=5 years"] = np.where(
        (data["parous"] == 0) | np.isnan(data["years_since_pregnancy"]),
        np.nan,
        data["years_since_pregnancy"] >= 5,
    )

    data["parity <10 years"] = np.where(
        (data["parous"] == 0) | np.isnan(data["years_since_pregnancy"]),
        np.nan,
        data["years_since_pregnancy"] < 10,
    )

    data["parity >=10 years"] = np.where(
        (data["parous"] == 0) | np.isnan(data["years_since_pregnancy"]),
        np.nan,
        data["years_since_pregnancy"] >= 10,
    )

    data["parity 5-10 years"] = np.where(
        (data["parous"] == 0) | np.isnan(data["years_since_pregnancy"]),
        np.nan,
        (data["years_since_pregnancy"] >= 5) & (data["years_since_pregnancy"] < 10),
    )

    return data


def generate_lrdata(df, parity_ref, parity_comp, feature_name):
    """
    This function takes a dataframe and
    generates the input and target of the logistic regression model
    based on the given parity comparison and feature name.
    """
    lrdata = pd.DataFrame(
        None, index=df.index, columns=[parity_comp, feature_name, "age", "fam_hx"]
    )
    lrdata[feature_name] = df[feature_name].astype(float)
    lrdata["age"] = df["age_at_diagnosis"].astype(float)
    lrdata["fam_hx"] = df["fam_hx"].astype(float)
    for i in df.index:
        lrdata.loc[i, parity_comp] = (
            0
            if df.loc[i, parity_ref] == 1
            else 1
            if df.loc[i, parity_comp] == 1
            else np.nan
        )
    # drop any rows with NaN
    lrdata.dropna(inplace=True, axis=0)
    # separate X and y
    X = lrdata[[parity_comp, "age", "fam_hx"]]
    y = lrdata[feature_name]
    # standard scale age
    scaler = StandardScaler()
    X["age"] = scaler.fit_transform(X["age"].values.reshape(-1, 1))
    return X, y


def get_oddsratio_ci(X, y, alpha=0.95, rep=5000):
    """
    This function takes the input and target of the logistic regression
    and returns the odds ratio, confidence intervals, and the p-value.
    """
    or1, or2, or3 = [], [], []
    i = 0
    for i in range(rep):
        X_bs, y_bs = resample(
            X, y, random_state=i, stratify=y
        )  # create bootstrap (bs) sample
        if ~np.all(y_bs == 0):
            lrm = LogisticRegression(penalty="l2", solver="lbfgs")
            lrm.fit(X_bs, y_bs)
            or1.append(np.exp(lrm.coef_[0][0]))
            or2.append(np.exp(lrm.coef_[0][1]))
            or3.append(np.exp(lrm.coef_[0][2]))
        else:
            continue
    oddsratios = [np.mean(or1), np.mean(or2), np.mean(or3)]
    # first get ci1
    ci_lower = ((1.0 - alpha) / 2.0) * 100
    ci_higher = (alpha + ((1.0 - alpha) / 2.0)) * 100
    ci, pvals = [], []
    for bs_sample in [or1, or2, or3]:
        lower = max(0.0, np.percentile(bs_sample, ci_lower))
        upper = np.percentile(bs_sample, ci_higher)
        ci.append((lower, upper))
        pvals.append(
            np.min([(np.array(bs_sample) < 1).mean(), (np.array(bs_sample) > 1).mean()])
            * 2
        )
    return oddsratios, ci, pvals
