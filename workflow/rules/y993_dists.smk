# Distances to Y993 (mostly for K18Ac)
# Data already aggregated, no input functions required


rule make_t33b_y993_ndx:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        y993_ndx="runs/{folder}/y993_dist.ndx",
    params:
        grp_1="r 993 & a OH",
        grp_2="r 18 & a OH",
    shell:
        """
        mkdir -p results/{wildcards.folder}/y993_dists/data
        echo -e "{params.grp_1} \n{params.grp_2} \nq" |
        gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_y993_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/y993_dist.ndx",
    output:
        "results/{folder}/y993_dists/data/y993_oh_k18ac.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_993_&_OH" plus group "r_18_&_OH"' \
                     -oall {output} -tu ns
        """


rule plot_t33b_y993_all_heatmaps:
    input:
        rules.make_t33b_y993_xvg.output,
    output:
        "results/{folder}/y993_dists/Y993_dist_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="Y993\u2013K18Ac",
        vmin=2,
        vmax=10,
    script:
        "../scripts/plot_single_heatmap.py"


rule plot_t33b_y993_ar_kde:
    input:
        rules.make_t33b_y993_xvg.output,
    output:
        "results/{folder}/y993_dists/Y993_dist_kde.png",
    params:
        xmin=2,
        xmax=10,
        xlabel="Distance (Å)",
    script:
        "../scripts/plot_kde.py"


rule make_t33b_y993_xvg_clip_10ns:
    input:
        xtc="runs/{folder}/combined_clip_10ns.xtc",
        ndx="runs/{folder}/y993_dist.ndx",
    output:
        "results/{folder}/y993_dists/data/clip_y993_oh_k18ac.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'group "r_993_&_OH" plus group "r_18_&_OH"' \
                     -oall {output} -tu ns
        """


rule plot_t33b_y993_kde_clip_10ns:
    input:
        rules.make_t33b_y993_xvg_clip_10ns.output,
    output:
        "results/{folder}/y993_dists/Y993_dist_kde_clip_10ns.png",
    params:
        xmin=2,
        xmax=10,
        xlabel="Distance (Å)",
    script:
        "../scripts/plot_kde.py"
