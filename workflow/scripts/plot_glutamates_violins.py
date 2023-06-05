import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

# y_data attribute is of interest here
# data in nm so multiply by 10 for angstrom
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

glu_1_dist = read_xvg_file(snakemake.input.glu_1_r17, dist=dist_on)
glu_2_dist = read_xvg_file(snakemake.input.glu_2_r17, dist=dist_on)
glu_3_dist = read_xvg_file(snakemake.input.glu_3_r17, dist=dist_on)
labels = [
    f"{snakemake.params.glu_1}\u2013T3A-R17",
    f"{snakemake.params.glu_2}\u2013T3A-R17",
    f"{snakemake.params.glu_3}\u2013T3A-R17",
]

data_list = [glu_1_dist, glu_2_dist, glu_3_dist]

df = pd.DataFrame(data_list).T
print(df.describe())
df.columns = labels

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)

sns.violinplot(data=df, scale="area")
ax.set(
    xlabel=None,
    ylabel="Distance (Å)",
    ylim=(2, 15),
)

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
