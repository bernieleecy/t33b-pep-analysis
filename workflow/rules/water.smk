# Getting water density


rule get_water_density:
    """
    Gets water density
    """
    input:
        tpr="runs/{folder}/em.tpr",
        xtc="runs/{folder}/combined_wat.xtc",
    output:
        "results/{folder}/water_density/water.dx",
    script:
        "../scripts/mda_water_density.py"


rule plot_water_density:
    """
    Plots water density in PyMOL
    Important to use complex_ions.pdb here
    - Waters removed, so the only waters are from complex.gro (the
      crystallographic waters)
    - remove complex_ions after aligning, otherwise I have many atoms
    """
    input:
        complex="runs/{folder}/complex.gro",
        complex_ions="runs/{folder}/complex_ions.pdb",
        water=rules.get_water_density.output,
    output:
        pymol_script="results/{folder}/water_density/water.pml",
        pse="results/{folder}/water_density/water_density.pse",
        png_20="results/{folder}/water_density/water_xtal_2.0.png",
        png_15="results/{folder}/water_density/water_xtal_1.5.png",
    shell:
        """
        echo -e "load {input.complex} \
                 \nload {input.complex_ions} \
                 \nload {input.water} \
                 \nalign complex, complex_ions \
                 \ndelete complex_ions \
                 \nextract bd_wat, resname HOH \
                 \nhide sticks lines \
                 \nshow lines, resi 18 and not name H* \
                 \nshow spheres, bd_wat and name OW \
                 \nset sphere_scale, 0.2, bd_wat \
                 \n \
                 \nisomesh mesh_2.0, water, 2.0, bd_wat \
                 \nisomesh mesh_1.5, water, 1.5, bd_wat \
                 \nisomesh mesh_all_2.0, water, 2.0 \
                 \nisomesh mesh_all_1.5, water, 1.5 \
                 \nhide mesh, mesh_all* \
                 \nset mesh_width, 1.5 \
                 \nset mesh_color, gray50 \
                 \nset_view (\\ \
                     \n0.634894371,    0.499791741,   -0.589161217, \\ \
                    \n-0.771749020,    0.374542743,   -0.513925195, \\ \
                    \n-0.036189321,    0.780976713,    0.623511136, \\ \
                    \n-0.000149934,   -0.000237733,  -59.691345215, \\ \
                    \n66.286834717,   61.318622589,   41.697376251, \\ \
                  \n-2813.448486328, 2932.944091797,   20.000000000 ) \
                 \n \
                 \nsel PHD, resi 885-934 \
                 \nsel linker, resi 935-958 \
                 \nsel bromodomain, resi 959-1071 \
                 \nsel peptide, resi 1-21 \
                 \n \
                 \nutil.cba('palegreen', 'PHD') \
                 \nutil.cba('helium', 'linker') \
                 \nutil.cba('lightpink', 'bromodomain') \
                 \nutil.cba('36', 'peptide') \
                 \n \
                 \nsave {output.pse} \
                 \nhide mesh, mesh_1.5 \
                 \npng {output.png_20}, ray=1, dpi=600, height=1000, width=1000 \
                 \nshow mesh, mesh_1.5 \
                 \nhide mesh, mesh_2.0 \
                 \npng {output.png_15}, ray=1, dpi=600, height=1000, width=1000 \
                 \ndeselect" > {output.pymol_script} 
        pymol -c {output.pymol_script} 
        """
