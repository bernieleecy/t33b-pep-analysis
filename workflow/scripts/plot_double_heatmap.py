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

    dist_data = {}

    with open(file, "r") as f:
        for line in f:
            if not line.startswith(("@", "#")):
                data = line.strip().split()
                time = float(data[0]) / 1000
                dist = float(data[1]) * 10  # convert to angstrom

                # values need to be initialised as a list ([dist])
                if time not in dist_data:
                    dist_data[time] = [dist]
                # subsequent values are then appended
                else:
                    dist_data[time].append(dist)

    df = pd.DataFrame.from_dict(dist_data)
    df["Run_ID"] = [f"Run {i+1}" for i in range(n_runs)]
    df.set_index("Run_ID", inplace=True)
    print(df.T.describe())

    return df


# Plot for benzylic amine distances
# Here, it makes more sense to have 2 side-by-side plots
# May need to have n_runs as a parameter in the future

df_1 = gen_df(snakemake.input.f1, n_runs=snakemake.params.n_runs)
df_2 = gen_df(snakemake.input.f2, n_runs=snakemake.params.n_runs)

fig, axes = plt.subplots(
    1, 3, figsize=(15, 5), gridspec_kw=dict(width_ratios=[1, 1, 0.05])
)

vmin = int(getattr(snakemake.params, "vmin", 2))
vmax = int(getattr(snakemake.params, "vmax", 15))

sns.heatmap(df_1, vmin=vmin, vmax=vmax, cmap="rocket", ax=axes[0], cbar=False)
sns.heatmap(df_2, vmin=vmin, vmax=vmax, cmap="rocket", ax=axes[1], cbar=False)

fig.colorbar(axes[1].collections[0], cax=axes[2])

xticks = np.arange(0, 5001, 1000)
xticklabels = np.arange(0, 101, 20)

axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
axes[0].set_title(f"{snakemake.params.pro_1}\u2013{snakemake.params.lig_1}")

axes[1].set_yticklabels("")
axes[1].set_title(f"{snakemake.params.pro_2}\u2013{snakemake.params.lig_2}")

for ax in [axes[0], axes[1]]:
    ax.set(xlabel="Time (ns)", ylabel="", xticks=xticks)
    ax.set_xticklabels(xticklabels, rotation=0)

    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1.5)

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=600)
