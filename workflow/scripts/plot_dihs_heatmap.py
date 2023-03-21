# For cMD runs, for dihedrals only

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
plt.style.use(snakemake.config["mpl_style"])
sns.set_palette("colorblind")


def gen_df(file, n_runs=3):
    """Makes the df for the heatmap from a xvg file

    Args:
        file (str): xvg file with time series data
        n_runs (int): Number of runs in the xvg file (default 3)

    Returns:
        df (pandas DataFrame): DataFrame for heatmap
    """

    dih_data = {}

    with open(file, "r") as f:
        for line in f:
            if not line.startswith(("@", "#")):
                data = line.strip().split()
                time = float(data[0]) / 1000
                dih = float(data[1]) # already in deg

                # values need to be initialised as a list ([dih])
                if time not in dih_data:
                    dih_data[time] = [dih]
                # subsequent values are then appended
                else:
                    dih_data[time].append(dih)

    df = pd.DataFrame.from_dict(dih_data)
    df["Run_ID"] = [f"Run {i+1}" for i in range(n_runs)]
    df.set_index("Run_ID", inplace=True)
    print(df.T.describe())
    print(df.shape)

    return df


data_df = gen_df(snakemake.input[0], n_runs=snakemake.params.n_runs)

fig, axes = plt.subplots(1, 2, figsize=(8, 5), gridspec_kw=dict(width_ratios=[1, 0.05]))

vmin = int(getattr(snakemake.params, "vmin", -180))
vmax = int(getattr(snakemake.params, "vmax", 180))
sns.heatmap(data_df, vmin=vmin, vmax=vmax, cmap="rocket", ax=axes[0], cbar=False)

fig.colorbar(axes[0].collections[0], cax=axes[1])

n_frames = data_df.shape[1]
sim_time = (n_frames * snakemake.config["xtc_step"]) / 1000 # in ns
xticks = np.arange(0, n_frames, (n_frames//5))
xticklabels = np.arange(0, sim_time, (sim_time//5), dtype=int)

axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
axes[0].set_title(snakemake.params.title)

axes[0].set(xlabel="Time (ns)", ylabel="", xticks=xticks)
axes[0].set_xticklabels(xticklabels, rotation=0)

for _, spine in axes[0].spines.items():
    spine.set_visible(True)
    spine.set_linewidth(1.5)

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=600)
