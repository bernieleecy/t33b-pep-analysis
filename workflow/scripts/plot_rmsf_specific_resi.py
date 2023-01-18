import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

files = sorted(snakemake.input)


def process_rmsf(file_list):
    """
    For processing a list of xvg files to extract RMSF

    As I want to use df.describe(), don't transpose the data

    Returns the summary statistics of the dataframe (Average is row 1, stdev is
    row 2)
    """
    resi = []
    rmsf = []

    for i, file in enumerate(file_list):
        with open(file) as infile:
            val = []

            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    data = line.strip().split("  ")  # two spaces between each entry
                    rmsf_val = float(data[1].strip())
                    val.append(rmsf_val * 10)  # convert to A
                    if i == 0:
                        resi.append(int(data[0].strip()))

            rmsf.append(val)

    df = pd.DataFrame(rmsf, columns=resi)

    return df


# all data, exclude N and C terminal caps
resi_of_interest = snakemake.params.resi_of_interest
rmsf_data = process_rmsf(files).filter(items=resi_of_interest)

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)

sns.boxplot(data=rmsf_data)

ymin = getattr(snakemake.params, "ymin", 0)
ymax = getattr(snakemake.params, "ymax", 5)

ax.set(
    xlabel="Residue number",
    ylabel="RMSF (Å)",
    xticklabels=resi_of_interest,
    ylim=(0, 3),
)

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
