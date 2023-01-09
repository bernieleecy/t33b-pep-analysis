"""
Time series plot (combining multiple runs)
Adpated for dealing with gromacs xvg files
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import itertools

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

# Snakemake input is a list of files here
files = sorted(snakemake.input)
n_runs = len(files)


def read_xvg_files(file_list, tu="ns"):
    """
    For processing multiple xvg files to extract distance data
    - Converts from nm to Å
    Also need the times, so I can plot ligand heavy atom RMSD against time
    """
    timesteps = []
    time_data = []

    for i, file in enumerate(file_list):
        single_data = []

        with open(file) as infile:
            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    data = line.strip().split("    ")  # four spaces between each entry
                    if tu == "ns":
                        time = float(data[0].strip())
                    elif tu == "ps":
                        time = float(data[0].strip()) / 1000
                    data_val = float(data[1].strip())
                    single_data.append(data_val * 10)  # convert to A
                    if i == 0:
                        timesteps.append(time)

        time_data.append(single_data)

    return timesteps, time_data


# get times and data
time_unit = str(getattr(snakemake.params, "time_unit", "ns"))
times, data = read_xvg_files(files)

# plot data at a specified interval (from config files)
# 100 ns trajectories here, so plot at 500 ps intervals
timestep = int(snakemake.config["xtc_step"])
plot_step = int(snakemake.config["plot_step"])
# must be an integer to use for slicing
interval = int(plot_step / timestep)

# plot data, (7,4) due to legend being outside the plot
fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)

for i, y in enumerate(data):
    ax.plot(times[::interval], y[::interval], label=f"Run {i+1}")

# 100 ns intervals so change to 5 ticks
xticks_intervals = int(times[-1] / 5)
xticks = np.arange(0, int(times[-1]) + 1, xticks_intervals)

# check for the ylabel param, otherwise supply a default
# can set ymin and ymax if desired
ylabel = str(getattr(snakemake.params, "ylabel", "Distance (Å)"))
ymin = int(getattr(snakemake.params, "ymin", 0))
ymax = int(getattr(snakemake.params, "ymax", 15))

ax.set(
    xlabel=f"Time ({time_unit})",
    ylabel=ylabel,
    ylim=(ymin, ymax),
    xticks=xticks,
)

ax.legend(frameon=False, bbox_to_anchor=(1, 0.5), loc="center left")
sns.despine()

fig.savefig(snakemake.output[0], dpi=600)
