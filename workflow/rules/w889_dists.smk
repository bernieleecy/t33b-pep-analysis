# For w889 distances, this is an optional input
# Data already aggregated, no input functions required


rule make_t33b_w889_ndx:
    input:
        "runs/{folder}/md_1.tpr",
    output:
        w889_ndx="runs/{folder}/w889_dist.ndx",
    params:
        grp_1="r 889 & a NE1 CD* CE* CG CH2 CZ*",
        grp_2="r 9 & a NZ",
    shell:
        """
        echo -e "{params.grp_1} \n{params.grp_2} \nq" |
        gmx make_ndx -f {input} -o {output}
        """


rule make_t33b_w889_xvg:
    input:
        xtc="runs/{folder}/combined.xtc",
        ndx="runs/{folder}/w889_dist.ndx",
    output:
        "results/{folder}/w889_dists/data/w889_to_k9me3.xvg",
    params:
        prefix="runs/{folder}",
    shell:
        """
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'com of group "r_889_&_NE1_CD*_CE*_CG_CH2_CZ*" plus group "r_9_&_NZ"' \
                     -oall {output} -tu ns
        """


rule plot_t33b_w889_k9me3_heatmap:
    input:
        rules.make_t33b_w889_xvg.output,
    output:
        "results/{folder}/w889_dists/w889_k9me3_dist_heatmap.png",
    params:
        n_runs=N_RUNS,
        title="W889\u2013K9 NZ",
        vmin=3,
        vmax=6,
    script:
        "../scripts/plot_single_heatmap.py"


rule plot_t33b_w889_k9me3_kde:
    input:
        rules.make_t33b_w889_xvg.output,
    output:
        "results/{folder}/w889_dists/w889_k9me3_dist_kde.png",
    params:
        xmin=3,
        xmax=6,
        xlabel="Distance (Å)",
        dist_on=True,
    script:
        "../scripts/plot_kde.py"
