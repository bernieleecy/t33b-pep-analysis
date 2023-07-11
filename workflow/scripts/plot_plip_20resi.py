# Plot PLIP data from a single excel file
# Not excluding any data from the PLIP analysis (analysed 505 frames)
# Adapted for use with Snakemake

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def make_df(file, start_res=1, end_res=20):
    """
    Reads an excel file containing the number of hydrophilic and hydrophobic
    residues (split by residues)
    Reads the first sheet

    Processes the data for the specified residues.
    - Supply the residue numbers as is (do not correct for python's
      zero-indexing, the function deals with it)
    """

    df = pd.read_excel(file, sheet_name=0, index_col=0)
    df.fillna(value=0, inplace=True)
    df.reset_index(inplace=True)
    df["Peptide_number"].astype("int")

    # minus one from start_res to handle the zero indexing
    df = df.iloc[start_res - 1 : end_res]
    melt_df = pd.melt(
        df,
        id_vars=["Peptide_number"],
        value_vars=["Hydrophobic", "Hydrophilic"],
        var_name="Contact type",
    )
    return melt_df


def plot_data(df, ax):
    total_contacts = df.groupby(["Peptide_number"], as_index=False).sum()
    hydrophobic = df[df["Contact type"] == "Hydrophobic"]
    sns.barplot(
        x="Peptide_number",
        y="value",
        data=total_contacts,
        ax=ax,
        palette=["lightsteelblue"],
        label="Hydrophilic",
    )
    sns.barplot(
        x="Peptide_number",
        y="value",
        data=hydrophobic,
        ax=ax,
        palette=["thistle"],
        label="Hydrophobic",
    )


# plot all data
df = make_df(snakemake.input[0], start_res=1, end_res=20)

if snakemake.params.k18_type == "K18":
    labels = [
        "A1",
        "R2",
        "T3",
        "K4",
        "Q5",
        "T6",
        "A7",
        "R8",
        "K9",
        "S10",
        "T11",
        "G12",
        "G13",
        "K14",
        "A15",
        "P16",
        "R17",
        "K18",
        "Q19",
        "L20",
    ]
elif snakemake.params.k18_type == "K18Ac":
    labels = [
        "A1",
        "R2",
        "T3",
        "K4",
        "Q5",
        "T6",
        "A7",
        "R8",
        "K9Me$_\mathregular{3}$",
        "S10",
        "T11",
        "G12",
        "G13",
        "K14",
        "A15",
        "P16",
        "R17",
        "K18Ac",
        "Q19",
        "L20",
    ]

fig, ax = plt.subplots(figsize=(10, 4.5), constrained_layout=True)

plot_data(df, ax=ax)

ymin = getattr(snakemake.params, "ymin", 0)
ymax = getattr(snakemake.params, "ymax", 1600)

ax.set(xlabel="", ylabel="", ylim=(ymin, ymax))
ax.set_ylabel("Count", fontweight="bold")
ax.set_xlabel("Peptide residue", fontweight="bold")
ax.set_xticklabels(labels, fontsize=11)

leg = ax.legend(frameon=False, loc="upper left")
leg._legend_box.align = "left"

fig.savefig(snakemake.output[0], dpi=600)
