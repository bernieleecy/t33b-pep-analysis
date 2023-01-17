"""
Time series plot (individual runs)
Adapted for dealing with gromacs xvg files
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])


def read_xvg_file(file, tu="ns"):
    """
    For processing a single xvg file to extract distance data
    - Converts from nm to Å
    Also need the times, so I can plot ligand heavy atom RMSD against time
    """
    timesteps = []
    time_data = []

    with open(file) as infile:
        for line in infile:
            if not (line.startswith("#") or line.startswith("@")):
                data = line.strip().split("    ")  # four spaces between each entry
                if tu == "ns":
                    time = float(data[0].strip())
                elif tu == "ps":
                    time = float(data[0].strip()) / 1000
                data_val = float(data[1].strip())
                timesteps.append(time)
                time_data.append(data_val * 10)  # convert to A

    return timesteps, time_data


# get times and data
time_unit = str(getattr(snakemake.params, "time_unit", "ns"))
time_file = snakemake.input[0]
times, data = read_xvg_file(time_file, tu=time_unit)

# plot data at a specified interval (from config files)
# 100 ns trajectories in this workflow, so plot at 500 ps intervals
timestep = int(snakemake.config["xtc_step"])
plot_step = int(snakemake.config["plot_step"])
# must be an integer to use for slicing
interval = int(plot_step / timestep)

# plot data
fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax.plot(times[::interval], data[::interval], color="black")

# 100 ns trajectories so change to 5 ticks
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

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
