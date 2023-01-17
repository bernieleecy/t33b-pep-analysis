# plotting a single kde plot

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
# GROMACS data in nm so multiply by 10 for angstrom
data = cyg.XvgFile(snakemake.input[0]).process_data().y_data * 10

fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)

sns.kdeplot(x=data)

xmin = float(getattr(snakemake.params, "xmin", 0))
xmax = float(getattr(snakemake.params, "xmax", 20))
xlabel = str(getattr(snakemake.params, "xlabel", "Distance (Å)"))

ax.set(xlabel=xlabel, xlim=(xmin, xmax))

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
