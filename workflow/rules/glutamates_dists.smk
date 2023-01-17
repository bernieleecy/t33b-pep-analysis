# For looking at R17 to E981, E984 and E1044 distances
# Usually using rules.name_of_rule.output[0] to refer to files for plotting


rule make_t33b_glutamates_ndx:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        "runs/{folder}/glutamate_dist.ndx",
    shell:
        """
        mkdir -p results/{wildcards.folder}/glutamate_dists/data
        echo -e "r 981 & a CD \n \
                 r 984 & a CD \n \
                 r 1044 & a CD \n \
                 r 17 & a CZ \n q" |
                 gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_e981_xvgs:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/e981_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_981_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule make_t33b_e984_xvgs:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/e984_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_984_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule make_t33b_e1044_xvgs:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/e1044_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1044_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule plot_glu_violins:
    """
    Output of previous rules needs to be supplied as a file not a list
    (this applies to all the plotting rules in this file)
    """
    input:
        glu_1_r17=rules.make_t33b_e981_xvgs.output[0],
        glu_2_r17=rules.make_t33b_e984_xvgs.output[0],
        glu_3_r17=rules.make_t33b_e1044_xvgs.output[0],
    output:
        "results/{folder}/glutamate_dists/glutamate_violins.png",
    params:
        glu_1="E981",
        glu_2="E984",
        glu_3="E1044",
    script:
        "../scripts/plot_glutamates_violins.py"


rule plot_t33b_e981_e1044_r17_heatmaps:
    """
    Make 2 heatmaps side-by-side
    Here, lig_1 and lig_2 atoms are the same (R17)
    """
    input:
        f1=rules.make_t33b_e981_xvgs.output[0],
        f2=rules.make_t33b_e1044_xvgs.output[0],
    output:
        "results/{folder}/glutamate_dists/E981_E1044_R17_heatmap.png",
    params:
        pro_1="E981",
        pro_2="E1044",
        lig_1="R17",
        lig_2="R17",
        vmin=2,
        vmax=15,
        n_runs=N_RUNS,
    script:
        "../scripts/plot_double_heatmap.py"


rule plot_t33b_e981_e984_r17_heatmaps:
    """
    Make 2 heatmaps side-by-side
    Here, lig_1 and lig_2 atoms are the same (R17)
    """
    input:
        f1=rules.make_t33b_e981_xvgs.output[0],
        f2=rules.make_t33b_e984_xvgs.output[0],
    output:
        "results/{folder}/glutamate_dists/E981_E984_R17_heatmap.png",
    params:
        pro_1="E981",
        pro_2="E984",
        lig_1="R17",
        lig_2="R17",
        vmin=2,
        vmax=15,
        n_runs=N_RUNS,
    script:
        "../scripts/plot_double_heatmap.py"


rule make_t33b_e981_xvgs_clip_10ns:
    input:
        xtc="runs/{folder}/combined_clip_10ns.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/clip_e981_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_981_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule make_t33b_e984_xvgs_clip_10ns:
    input:
        xtc="runs/{folder}/combined_clip_10ns.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/clip_e984_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_984_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule make_t33b_e1044_xvgs_clip_10ns:
    input:
        xtc="runs/{folder}/combined_clip_10ns.xtc",
        ndx="runs/{folder}/glutamate_dist.ndx",
    output:
        "results/{folder}/glutamate_dists/data/clip_e1044_r17_cz.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1044_&_CD" plus group "r_17_&_CZ"' \
                     -oall {output} -tu ns
        """


rule plot_glu_violins_clip_10ns:
    input:
        glu_1_r17=rules.make_t33b_e981_xvgs_clip_10ns.output[0],
        glu_2_r17=rules.make_t33b_e984_xvgs_clip_10ns.output[0],
        glu_3_r17=rules.make_t33b_e1044_xvgs_clip_10ns.output[0],
    output:
        "results/{folder}/glutamate_dists/glutamate_violins_clip_10ns.png",
    params:
        glu_1="E981",
        glu_2="E984",
        glu_3="E1044",
    script:
        "../scripts/plot_glutamates_violins.py"
