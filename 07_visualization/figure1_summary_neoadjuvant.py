# pylint: disable=missing-module-docstring

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

din = "/share/fsmresfiles/breast_cancer_pregnancy/summary_tables"
dout = "/share/fsmresfiles/breast_cancer_pregnancy/plots"
fin = "20240202_super_summary_table.csv"
fout = "figure1_rcb_category.png"

data = pd.read_csv(f"{din}/{fin}")

fig, ax = plt.subplots(figsize=(3, 3))
pastel_colors = list(plt.cm.Pastel1.colors[:4])  # pylint: disable=no-member

rcb_dict = {"0": "RCB-0/I", "1": "RCB-0/I", "2": "RCB-II/III", "3": "RCB-II/III"}
# rcb_dict = {'0' : 'RCB-0', '1' : 'RCB-I', '2' : 'RCB-II', '3' : 'RCB-III'}

cts = {}
for category in ["Nulliparous", "<5 years", "5-10 years", ">=10 years"]:
    cts[category] = [
        np.sum(
            [
                int(x.split()[0])
                for x in data.iloc[
                    np.array([x in varnames for x in data["sub_category"].values])
                    & np.array(
                        [data["main_category"].values == "rcb_category"]
                    ).ravel(),
                    :,
                ][category].values
            ]
        )
        for varnames in [["0", "1"], ["2", "3"]]
    ]

ax.bar(
    range(len(cts["Nulliparous"])),
    cts["Nulliparous"],
    label="Nulliparous",
    color=pastel_colors[0],
    edgecolor="black",
    linewidth=0.5,
)
ax.bar(
    range(len(cts["Nulliparous"])),
    cts["<5 years"],
    bottom=cts["Nulliparous"],
    label="<5 years",
    color=pastel_colors[1],
    edgecolor="black",
    linewidth=0.5,
)
ax.bar(
    range(len(cts["Nulliparous"])),
    cts["5-10 years"],
    bottom=[np.sum([x, y]) for x, y in zip(cts["Nulliparous"], cts["<5 years"])],
    label="5-10 years",
    color=pastel_colors[2],
    edgecolor="black",
    linewidth=0.5,
)
ax.bar(
    range(len(cts["Nulliparous"])),
    cts[">=10 years"],
    bottom=[
        np.sum([x, y, z])
        for x, y, z in zip(cts["Nulliparous"], cts["<5 years"], cts["5-10 years"])
    ],
    label=">=10 years",
    color=pastel_colors[3],
    edgecolor="black",
    linewidth=0.5,
)
ax.set_ylabel("Number of patients")
ax.set_title("RCB Category Post-NAC")
ax.set_xticks(range(len(cts["Nulliparous"])))
ax.set_xticklabels(["RCB-0/I", "RCB-II/III"], rotation=45, ha="right")
# ax.set_xticklabels(['RCB-0', 'RCB-I', 'RCB-II', 'RCB-III'], rotation=45, ha='right')

## Remove the top and right spines for all subplots
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{dout}/{fout}")
plt.close()
