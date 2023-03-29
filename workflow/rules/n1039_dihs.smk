# For dihedral angles, e.g. n1039 chi1
# Optional input, since chi1 angles are not relevant to the N1039A mutant
# Data already aggregated, no input functions required


rule make_t33b_chi:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        "runs/{folder}/n1039_dihs.ndx",
    params:
        grp_1="r 1039 & a N CA CB CG",
        grp_2="r 1039 & a CA CB CG OD1",
    shell:
        """
        mkdir -p results/{wildcards.folder}/n1039_dihs/data
        echo -e "{params.grp_1} \n{params.grp_2} \nq" |
        gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_n1039_dihs_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx=rules.make_t33b_chi.output,
    output:
        dih_1_hist = "results/{folder}/n1039_dihs/data/n1039_chi1_hist.xvg",
        dih_1_v_time = "results/{folder}/n1039_dihs/data/n1039_chi1.xvg",
        dih_2_hist = "results/{folder}/n1039_dihs/data/n1039_chi2_hist.xvg",
        dih_2_v_time = "results/{folder}/n1039_dihs/data/n1039_chi2.xvg",
    params:
        ndx_group_1="r_1039_&_N_CA_CB_CG",
        ndx_group_2="r_1039_&_CA_CB_CG_OD1",
    shell:
        """
        echo '{params.ndx_group_1}' |
        gmx angle -f {input.xtc} -n {input.ndx} \
                -od {output.dih_1_hist} \
                -ov {output.dih_1_v_time} -type dihedral
        echo '{params.ndx_group_2}' |
        gmx angle -f {input.xtc} -n {input.ndx} \
                -od {output.dih_2_hist} \
                -ov {output.dih_2_v_time} -type dihedral
        """


rule plot_t33b_n1039_chi_all_heatmaps:
    input:
        "results/{folder}/n1039_dihs/data/n1039_chi{x}.xvg",
    output:
        "results/{folder}/n1039_dihs/N1039_chi{x}_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="N1039 chi{x} angle",
        vmin=-180,
        vmax=180,
        dist_on=False,
    script:
        "../scripts/plot_single_heatmap.py"


rule plot_t33b_n1039_chi_all_hist:
    input:
        "results/{folder}/n1039_dihs/data/n1039_chi{x}.xvg",
    output:
        "results/{folder}/n1039_dihs/N1039_chi{x}_hist.png",
    params:
        bins = 36,
        xmin = -180,
        xmax = 180,
        xlabel = "chi{x} angle (°)",
        ymax = 40000,
    script:
        "../scripts/plot_hist.py"
