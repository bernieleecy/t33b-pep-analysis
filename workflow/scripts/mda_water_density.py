# Using MDAnalysis to obtain water density relative to TIP3P water

import MDAnalysis as mda
from MDAnalysis.analysis.density import DensityAnalysis

tpr_file = snakemake.input.tpr
xtc_file = snakemake.input.xtc

u = mda.Universe(tpr_file, xtc_file)

ow = u.select_atoms("name OW")

D = DensityAnalysis(ow, delta=1)
D.run()
D.results.density.convert_density("TIP3P")

# type='double' ensures the dx file can be read by PyMOL
D.results.density.export(snakemake.output[0], type="double")
