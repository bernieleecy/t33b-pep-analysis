import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cygnus as cyg
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config["mpl_style"])

# y_data attribute is of interest here
# data in nm so multiply by 10 for angstrom

glu_1_dist = cyg.XvgFile(snakemake.input.glu_1_r17).process_data().y_data * 10
glu_2_dist = cyg.XvgFile(snakemake.input.glu_2_r17).process_data().y_data * 10
glu_3_dist = cyg.XvgFile(snakemake.input.glu_3_r17).process_data().y_data * 10
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
