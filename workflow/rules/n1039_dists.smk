# For n1039 distances, this is an optional input
# Data already aggregated, no input functions required


rule make_t33b_n1039_ndx:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        n1039_ndx="runs/{folder}/n1039_dist.ndx",
    params:
        grp_1="r 1039 & a ND2",
        grp_2="r 18 & a OH",
        grp_3="r 981 & a OE1",
        grp_4="r 981 & a OE2",
    shell:
        """
        echo -e "{params.grp_1} \n{params.grp_2} \n{params.grp_3} \n{params.grp_4} \nq" |
        gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_n1039_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/n1039_dist.ndx",
    output:
        "results/{folder}/n1039_dists/data/n1039_nd2_k18ac.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1039_&_ND2" plus group "r_18_&_OH"' \
                     -oall {output} -tu ns
        """


rule make_t33b_n1039_e981_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/n1039_dist.ndx",
    output:
        n1039_e981_oe1="results/{folder}/n1039_dists/data/n1039_nd2_e981_oe1.xvg",
        n1039_e981_oe2="results/{folder}/n1039_dists/data/n1039_nd2_e981_oe2.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1039_&_ND2" plus group "r_981_&_OE1"' \
                     -oall {output.n1039_e981_oe1} -tu ns
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1039_&_ND2" plus group "r_981_&_OE2"' \
                     -oall {output.n1039_e981_oe2} -tu ns

        """


rule plot_t33b_n1039_k18ac_heatmap:
    input:
        rules.make_t33b_n1039_xvg.output,
    output:
        "results/{folder}/n1039_dists/N1039_k18ac_dist_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="N1039\u2013K18Ac",
        vmin=2,
        vmax=10,
    script:
        "../scripts/plot_single_heatmap.py"


rule plot_t33b_n1039_k18ac_kde:
    input:
        rules.make_t33b_n1039_xvg.output,
    output:
        "results/{folder}/n1039_dists/N1039_k18ac_dist_kde.png",
    params:
        xmin=2,
        xmax=10,
        xlabel="Distance (Å)",
    script:
        "../scripts/plot_kde.py"


rule plot_t33b_n1039_e981_heatmaps:
    input:
        "results/{folder}/n1039_dists/data/n1039_nd2_e981_{x}.xvg",
    output:
        "results/{folder}/n1039_dists/N1039_E981_{x}_dist_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="N1039\u2013E981 ({x})",
        vmin=2,
        vmax=20,
    script:
        "../scripts/plot_single_heatmap.py"


rule make_t33b_n1039_xvg_clip_10ns:
    input:
        xtc="runs/{folder}/combined_clip_10ns.xtc",
        ndx="runs/{folder}/n1039_dist.ndx",
    output:
        "results/{folder}/n1039_dists/data/clip_n1039_nd2_k18ac.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_1039_&_ND2" plus group "r_18_&_OH"' \
                     -oall {output} -tu ns
        """


rule plot_t33b_n1039_kde_clip_10ns:
    input:
        rules.make_t33b_n1039_xvg_clip_10ns.output,
    output:
        "results/{folder}/n1039_dists/N1039_k18ac_dist_kde_clip_10ns.png",
    params:
        xmin=2,
        xmax=10,
        xlabel="Distance (Å)",
    script:
        "../scripts/plot_kde.py"
