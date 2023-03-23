# For dihedral angles, e.g. n1039 chi1
# Optional input, since chi1 angles are not relevant to the N1039A mutant
# Data already aggregated, no input functions required


rule make_t33b_chi1:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        n1039_ndx="runs/{folder}/n1039_dihs.ndx",
    params:
        grp_1="r 1039 & a N CA CB CG",
    shell:
        """
        mkdir -p results/{wildcards.folder}/n1039_dihs/data
        echo -e "{params.grp_1} \nq" |
        gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_n1039_dihs_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/n1039_dihs.ndx",
    output:
        dih_hist = "results/{folder}/n1039_dihs/data/n1039_chi1_hist.xvg",
        dih_v_time = "results/{folder}/n1039_dihs/data/n1039_chi1.xvg",
    params:
        ndx_group="r_1039_&_N_CA_CB_CG"
    shell:
        """
        echo '{params.ndx_group}' |
        gmx angle -f {input.xtc} -n {input.ndx} \
                -od {output.dih_hist} \
                -ov {output.dih_v_time} -type dihedral
        """


rule plot_t33b_n1039_chi1_all_heatmaps:
    input:
        rules.make_t33b_n1039_dihs_xvg.output.dih_v_time,
    output:
        "results/{folder}/n1039_dihs/N1039_chi1_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="N1039 chi1 angle",
        vmin=-180,
        vmax=180,
        dist_on=False,
    script:
        "../scripts/plot_single_heatmap.py"


rule plot_t33b_n1039_chi1_all_hist:
    input:
        rules.make_t33b_n1039_dihs_xvg.output.dih_v_time,
    output:
        "results/{folder}/n1039_dihs/N1039_chi1_hist.png",
    params:
        bins = 36,
        xmin = -180,
        xmax = 180,
        xlabel = "chi1 angle (°)",
        ymax = 40000,
    script:
        "../scripts/plot_hist.py"
