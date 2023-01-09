# Snakemake workflow: t33b-pep-analysis

A Snakemake workflow for analysing T33B—peptide MD simulations, initialised on 9 January 2023 (and is thus under active development).
* Currently using snakemake >=7.8, as this reruns code automatically if the parameters, code or input file set has changed – this saves a lot of typing

On my M1 mac, the appropriate conda environment is `snake_env`
* I have a configuration file in the `TRIM3` folder that activates the appropriate conda env and sets `$GMXLIB`, so that should be run before using this workflow
* As of 11/08/2022, I'm using a new version of `snake_env` that also includes AmberTools22 (uses intel packages)
* This version of `snake_env` has the `osx-64` subdir set in the `.condarc` file, see install instructions below

```
CONDA_SUBDIR=osx-64 mamba env create -f t33b_md.yml
conda activate snake_env
python -c "import platform; print(platform.machine())" # this should be x86_64
conda config --env --set subdir osx-64
```

## Organisation

Folder structure on my computer, with run and results folders shown as an example. 

```
.
├── config
├── results
│   ├── apo
│   ├── mod_pep
│   └── unmod_pep
├── runs
│   ├── apo
│   ├── mod_pep
│   └── unmod_pep
└── workflow
    ├── rules
    └── scripts
```
# Usage

## Files required

Files needed in the `runs/run_type` folder:

These files do not need to be generated separately: 
- `md_{id}.xtc`: Raw trajectory files from GROMACS 
- `complex.gro`
- `complex_ions.gro`
- `em.tpr`
- `md_1.tpr` (rather than `md_1_extend.tpr`)

Trying to remove dependency on doing basic trajectory processing, so working on incorporating the PBC correction and trajectory fitting into the workflow. 

These files are generated when I do basic trajectory processing after a run (before moving the runs into the workflow folder): 
- `{id}-whole_fit.xtc`: PBC-corrected, solvent and ions removed, aligned to `md_{id}.tpr` protein backbone
- `index.ndx`: Contains the 'Protein_T3A_ZN' and 'T3A&!H*' index groups 

## Filling in the samples.csv files

Details for the samples.csv file:
* folder: contains all the runs and the required tpr file, etc.
* peptide: whether the runs contain a peptide (yes for most, no for apo runs)
* k9_type: modification on k9 (K9Me3, K9, NA)
* k18_type: modification on k18 (K18Ac, K18, na)
* mutations: whether the protein has a mutation (WT for no mutations, mutations of interest are N1039A and E981A)

## Running rules

To run the basic analysis, with parallelisation:
```
snakemake -np # dry run first
snakemake -c4 # run with 4 cores
```

To rerun some analysis with a reduced dataset
```
snakemake clip_10ns_analysis -c4
```

The PLIP rule is separated from the others, and can be run with:
```
snakemake get_plip -c4
```

# Other notes 

To configure the environment on a new machine, run
```
conda env create -f t33a_md.yml
```

In addition, GROMACS is required to use the `gmx` tools.
PLIP is also needed (I use Docker to install it). 
