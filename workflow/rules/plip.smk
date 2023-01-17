import os


checkpoint make_plip_pdbs:
    input:
        xtc="runs/{folder}/plip_traj.xtc",
        ndx="runs/{folder}/index.ndx",
    output:
        directory("results/{folder}/plip/pdbs"),
    params:
        prefix="runs/{folder}",
    shell:
        """
        mkdir -p {output}
        echo 'Protein_ZN' |
        gmx trjconv -f {input.xtc} -s {params.prefix}/{config[em_tpr]} \
                    -n {input.ndx} \
                    -o {output}/extracted.pdb -sep -nzero 3
        sed -i '' '/TITLE/d;/MODEL/d;/ENDMDL/d' {output}/extracted*.pdb
        sed -i '' 's/CYZ/CYS/g;s/HEZ/HIS/g;s/HIE/HIS/g;s/HID/HIS/g' {output}/extracted*.pdb
        """


rule run_plip:
    """
    Run PLIP on the extracted PDB files
    Wildcard constraint required to prevent weird recursive things
    """
    input:
        "results/{folder}/plip/pdbs/{i}.pdb",
    output:
        outdir=directory("results/{folder}/plip/pdbs/{i}"),
    wildcard_constraints:
        i="[A-Za-z0-9]+",
    shell:
        """
        docker run --rm -v $PWD:/results -w /results \
                   pharmai/plip:latest \
                   -f {input} -o {output.outdir} \
                   -x --nofixfile --hbond_dist_max 3.5 \
                   --hbond_don_angle_min 150 \
                   --peptides B
        """


def agg_plip(wildcards):
    """
    Forcing alphanumeric names (so I only search the desired directory for
    folders)
    Otherwise, if I have run PLIP before, the wildcard goes weird, e.g. it
    becomes "extracted000/extracted000_protonated" rather than "extracted000"
    """
    checkpoint_output = checkpoints.make_plip_pdbs.get(**wildcards).output[0]
    return expand(
        "results/{folder}/plip/pdbs/{i}",
        folder=wildcards.folder,
        i=glob_wildcards(os.path.join(checkpoint_output, "{i,[A-Za-z0-9]+}.pdb")).i,
    )


rule process_plip:
    """
    For aggregating PLIP results from report.xml files
    This will not rerun unless make_plip_pdbs reruns (due to an xtc file
    update), or if the script is updated
    As far as I can tell, this is intended snakemake behaviour!
    """
    input:
        agg_plip,
    output:
        "results/{folder}/plip/plip_data.xlsx",
    script:
        "../scripts/process_plip_data_peptides.py"


rule plot_plip_all:
    """
    Plot hydrophobic and hydrophilic contacts across all residues
    """
    input:
        rules.process_plip.output,
    output:
        "results/{folder}/plip/all_resi_hydrophobic_hydrophilic.png",
    params:
        k18_type=lambda wildcards: samples.loc[wildcards.folder, "k18_type"],
        ymax=2000,
    script:
        "../scripts/plot_plip_20resi.py"
