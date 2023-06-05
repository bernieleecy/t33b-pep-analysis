"""
Script to plot T33B ligand time series data

Adapted for use in a SnakeMake workflow
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

def read_xvg_file(file, dist=False):
    """
    For processing a single xvg file to extract data for histogram

    """
    time_data = []

    with open(file) as infile:
        for line in infile:
            if not (line.startswith("#") or line.startswith("@")):
                data = line.strip().split()
                data_val = float(data[1].strip())
                if dist:
                    time_data.append(data_val * 10)  # convert to A
                else:
                    time_data.append(data_val)

    return time_data

dist_on = bool(getattr(snakemake.params, "dist_on", False))
time_file = snakemake.input[0]
data = read_xvg_file(time_file, dist=dist_on)

# plot data
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)

# need np.append as amber time data starts from frame 1
bins = int(getattr(snakemake.params, "bins", 20))
plt.hist(data, bins=bins)

# check for the xmax and xlabel params, otherwise supply a default
xmin = float(getattr(snakemake.params, "xmin", 0))
xmax = float(getattr(snakemake.params, "xmax", 10))
ymin = float(getattr(snakemake.params, "ymin", 0))
ymax = float(getattr(snakemake.params, "ymax", 30000))
xlabel = str(getattr(snakemake.params, "xlabel", "Dihedral (°)"))

# setting ylim so I can compare different datasets
# use a param to set it, since I use this script for both protein and lig rmsd
ax.set(xlabel=xlabel, ylabel="Frequency", xlim=(xmin, xmax), ylim=(ymin, ymax))

sns.despine()

fig.savefig(snakemake.output[0], dpi=600)
