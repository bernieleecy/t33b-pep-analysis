# For cMD runs, for distances only, now modified to take variable length simulations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
plt.style.use(snakemake.config["mpl_style"])
sns.set_palette("colorblind")


def gen_df(file, n_runs=3, dist=True):
    """Makes the df for the heatmap from a xvg file

    Args:
        file (str): xvg file with time series data
        n_runs (int): Number of runs in the xvg file (default 3)

    Returns:
        df (pandas DataFrame): DataFrame for heatmap
    """

    gmx_data = {}

    with open(file, "r") as f:
        for line in f:
            if not line.startswith(("@", "#")):
                data = line.strip().split()
                time = float(data[0]) / 1000
                if dist:
                    gmx_val = float(data[1]) * 10  # convert to angstrom
                else:
                    gmx_val = float(data[1])

                # values need to be initialised as a list ([dist])
                if time not in gmx_data:
                    gmx_data[time] = [gmx_val]
                # subsequent values are then appended
                else:
                    gmx_data[time].append(gmx_val)

    df = pd.DataFrame.from_dict(gmx_data)
    df["Run_ID"] = [f"Run {i+1}" for i in range(n_runs)]
    df.set_index("Run_ID", inplace=True)
    print(df.T.describe())

    return df


dist_on = bool(getattr(snakemake.params, "dist_on", True))
data_df = gen_df(snakemake.input[0], n_runs=snakemake.params.n_runs, dist=dist_on)

fig, axes = plt.subplots(1, 2, figsize=(8, 5), gridspec_kw=dict(width_ratios=[1, 0.05]))

vmin = int(getattr(snakemake.params, "vmin", 2))
vmax = int(getattr(snakemake.params, "vmax", 15))
sns.heatmap(data_df, vmin=vmin, vmax=vmax, cmap="rocket", ax=axes[0], cbar=False)

label = getattr(snakemake.params, "label", "Distance (Å)")
fig.colorbar(axes[0].collections[0], cax=axes[1], label=label)

n_frames = data_df.shape[1]
sim_time = (n_frames * snakemake.config["xtc_step"]) / 1000  # in ns

if sim_time < 500:
    xticks = np.arange(0, n_frames, (n_frames // 5))
    xticklabels = np.arange(0, sim_time, (sim_time // 5), dtype=int)
else:
    xticks = np.arange(0, n_frames, (n_frames // 10))
    xticklabels = np.arange(0, sim_time, (sim_time // 10), dtype=int)

axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
axes[0].set_title(snakemake.params.title)

axes[0].set(xlabel="Time (ns)", ylabel="", xticks=xticks)
axes[0].set_xticklabels(xticklabels, rotation=0)

for _, spine in axes[0].spines.items():
    spine.set_visible(True)
    spine.set_linewidth(1.5)

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=600)
