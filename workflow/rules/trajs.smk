# For trajectory handling
# 1-whole_fit.xtc etc. are aligned to md_1_extend.tpr etc., not em.tpr
# need to redo alignment for clustering
import os


rule make_ndx:
    """
    Params depend on the protein
    """
    input:
        "runs/{folder}/md_1.tpr",
    output:
        "runs/{folder}/index.ndx",
    params:
        grp_1="1 | 13",
        grp_2="ri 1-188",
        grp_3="ri 1-188 & 4",
        grp_4="ri 189-209",
        grp_5="ri 189-209 & 4",
        unmod_4="ri 191-211",
        unmod_5="ri 191-211 & 4",
    run:
        if (
            samples.loc[wildcards.folder, "peptide"] == "yes"
            and samples.loc[wildcards.folder, "k18_type"] == "K18Ac"
        ):
            shell(
                'echo -e "{params.grp_1}\n{params.grp_2}\n{params.grp_3}\n{params.grp_4}\n{params.grp_5}\n q "'
                "| gmx make_ndx -f {input} -o {output}"
            )
        elif (
            samples.loc[wildcards.folder, "peptide"] == "yes"
            and samples.loc[wildcards.folder, "k18_type"] == "K18"
        ):
            shell(
                'echo -e "{params.grp_1}\n{params.grp_2}\n{params.grp_3}\n{params.unmod_4}\n{params.unmod_5}\n \
                name 23 r_189-209\n name 24 r_189-209_&_Backbone\n q "'
                "| gmx make_ndx -f {input} -o {output}"
            )
        elif samples.loc[wildcards.folder, "peptide"] == "no":
            shell(
                'echo -e "{params.grp_1}\n{params.grp_2}\n{params.grp_3}\n q "'
                "| gmx make_ndx -f {input} -o {output}"
            )


rule convert_gro:
    """
    For water density plotting
    """
    input:
        gro="runs/{folder}/complex_ions.gro",
        ndx="runs/{folder}/index.ndx",
    output:
        "runs/{folder}/complex_ions.pdb",
    params:
        prefix="runs/{folder}",
    shell:
        """
        echo "Protein_ZN" |
        gmx editconf -f {input.gro} -n {input.ndx} \
                     -o {output}
        """


rule make_protein_noh2o:
    """
    Just need one file to help with VMD visualisation
    """
    input:
        gro="runs/{folder}/md_1.gro",
        ndx="runs/{folder}/index.ndx",
    output:
        "runs/{folder}/protein_noh2o.gro",
    shell:
        "echo Protein_ZN | gmx editconf -f {input.gro} -n {input.ndx} -o {output}"


rule make_fitted_trajs:
    """
    Fitted to protein backbone (rather than protein and peptide backbone)
    """
    input:
        traj="runs/{folder}/md_{i}.xtc",
        ndx="runs/{folder}/index.ndx",
    output:
        whole=temp("runs/{folder}/{i}-whole.xtc"),
        whole_fit="runs/{folder}/{i}-whole_fit.xtc",
    params:
        prefix="runs/{folder}",
    shell:
        """
        echo -e "Protein_ZN \nProtein_ZN \nProtein_ZN" |
        gmx trjconv -f {input.traj} -s {params.prefix}/{config[md_tpr]} \
                -n {input.ndx} -o {output.whole} -pbc cluster -center
        echo -e "r_1-188_&_Backbone \nProtein_ZN" |
        gmx trjconv -f {output.whole} -s {params.prefix}/{config[md_tpr]} \
                -n {input.ndx} -o {output.whole_fit} -fit rot+trans
        """


def get_whole_trajs(wildcards):
    return [
        os.path.join("runs", wildcards.folder, run + "-whole_fit.xtc") for run in IDS
    ]


rule make_combined_traj:
    input:
        get_whole_trajs,
    output:
        "runs/{folder}/combined.xtc",
    shell:
        "gmx trjcat -f {input} -o {output} -cat"


rule make_combined_traj_clip_10ns:
    input:
        get_whole_trajs,
    output:
        "runs/{folder}/combined_clip_10ns.xtc",
    shell:
        "gmx trjcat -f {input} -o {output} -cat -b 10000"


rule make_plip_traj:
    input:
        get_whole_trajs,
    output:
        "runs/{folder}/plip_traj.xtc",
    shell:
        "gmx trjcat -f {input} -o {output} -cat -dt 1000"


rule align_wat_trajs:
    """
    Previously used Backbone for final fitting

    Here, use the protein backbone only for final fitting as this is more appropriate
    """
    input:
        xtc="runs/{folder}/md_{i}.xtc",
        ndx="runs/{folder}/index.ndx",
    output:
        align="runs/{folder}/{i}-align.xtc",
        traj_a=temporary("runs/{folder}/a_{i}.xtc"),
        traj_b=temporary("runs/{folder}/b_{i}.xtc"),
        traj_c=temporary("runs/{folder}/c_{i}.xtc"),
    params:
        prefix="runs/{folder}",
    shell:
        """
        echo 'Protein_ZN' 0 |
        gmx trjconv -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                    -n {input.ndx} -o {output.traj_a} \
                    -pbc cluster -dt 200

        echo 0 |
        gmx trjconv -f {output.traj_a} -s {params.prefix}/{config[md_tpr]} \
                    -o {output.traj_b} -pbc whole -ur compact

        echo 'Protein_ZN' 0 |
        gmx trjconv -f {output.traj_b} -s {params.prefix}/{config[md_tpr]} \
                    -n {input.ndx} -o {output.traj_c} \
                    -pbc mol -center -ur compact

       echo 'r_1-188_&_Backbone' 0 |
        gmx trjconv -f {output.traj_c} -s {params.prefix}/{config[em_tpr]} \
                    -n {input.ndx} -o {output.align} \
                    -fit rot+trans
        """


def get_aligned_wat_trajs(wildcards):
    return [os.path.join("runs", wildcards.folder, run + "-align.xtc") for run in IDS]


rule make_concat_wat_trajs:
    """
    Input removed after combining files to save space
    """
    input:
        get_aligned_wat_trajs,
    output:
        wat_full="runs/{folder}/combined_wat.xtc",
    shell:
        """
        gmx trjcat -f {input} -o {output.wat_full} -cat
        rm {input}
        """
