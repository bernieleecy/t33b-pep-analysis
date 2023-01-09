import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

sns.set_style("ticks")
sns.set_palette("colorblind")
plt.style.use(snakemake.config['mpl_style'])

files = sorted(snakemake.input)

def process_rmsf(file_list):
    '''
    For processing a list of xvg files to extract RMSF

    As I want to use df.describe(), don't transpose the data

    Returns the summary statistics of the dataframe (Average is row 1, stdev is
    row 2)
    '''
    rmsf = []

    for file in file_list:
        with open(file) as infile:
            val = []

            for line in infile:
                if not (line.startswith("#") or line.startswith("@")):
                    data = line.strip().split('  ') # two spaces between each entry
                    rmsf_val = float(data[1].strip())
                    val.append(rmsf_val * 10) # convert to A

            rmsf.append(val)

    df = pd.DataFrame(rmsf).describe()

    return df

# all data, exclude C terminal cap
rmsf_data = process_rmsf(files)
protein_rmsf_avg = rmsf_data.iloc[1,0:-1]
protein_rmsf_stdev = rmsf_data.iloc[2,0:-1]

# peptide residues hard coded in
resi = np.arange(1,21)
print(len(resi))

fig, ax = plt.subplots(figsize=(8,4), constrained_layout=True)

ax.errorbar(resi, protein_rmsf_avg, yerr=protein_rmsf_stdev,
            marker='o', capsize=3,
            color='black')

ax.set(xlabel='Residue number',
       ylabel='RMSF (Å)',
       xticks=resi,
       xticklabels=resi,
       ylim=(0,10))

sns.despine()
fig.savefig(snakemake.output[0], dpi=600)
