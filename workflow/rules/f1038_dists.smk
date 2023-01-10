# For F1038 distances
# Data already aggregated, no input functions required

rule make_t33b_f1038_ndx:
    input:
        'runs/{folder}/md_1.tpr'
    output:
        'runs/{folder}/f1038_dist.ndx',
    params:
        grp_1 = 'r 1038 & a CD* CE* CG CZ',
        grp_2 = 'r 993 & a CD* CE* CG CZ',
    shell:
        '''
        echo -e "{params.grp_1} \n{params.grp_2} \nq" |
        gmx make_ndx -f {input} -o {output}
        '''


rule make_t33b_f1038_xvgs:
    input:
        xtc = 'runs/{folder}/combined.xtc',
        ndx = 'runs/{folder}/f1038_dist.ndx',
    output:
        'results/{folder}/f1038_dists/data/f1038_y993.xvg',
    params:
        prefix = 'runs/{folder}'
    shell:
        '''
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'com of group "r_1038_&_CD*_CE*_CG_CZ" plus com of group "r_993_&_CD*_CE*_CG_CZ"' \
                     -oall {output} -tu ns
        '''

rule plot_t33b_f1038_y993_heatmap:
    input:
        rules.make_t33b_f1038_xvgs.output
    output:
        'results/{folder}/f1038_dists/F1038_Y993_heatmap.png',
    params:
        vmin = 3.5,
        vmax = 12,
        title = "Y993\u2013F1038",
        n_runs = N_RUNS,
    script:
        '../scripts/plot_single_heatmap.py'

rule plot_t33b_y993_f1038_kde:
    input:
        rules.make_t33b_f1038_xvgs.output
    output:
        'results/{folder}/f1038_dists/F1038_Y993_kde.png',
    params:
        xmin = 3.5,
        xmax = 12,
    script:
        '../scripts/plot_kde.py'


rule make_t33b_f1038_xvgs_clip_10ns:
    input:
        xtc = 'runs/{folder}/combined_clip_10ns.xtc',
        ndx = 'runs/{folder}/f1038_dist.ndx',
    output:
        'results/{folder}/f1038_dists/data/clip_f1038_y993.xvg',
    params:
        prefix = 'runs/{folder}'
    shell:
        '''
        gmx distance -f {input.xtc} -s {params.prefix}/{config[md_tpr]} \
                     -n {input.ndx} \
                     -select 'com of group "r_1038_&_CD*_CE*_CG_CZ" plus com of group "r_993_&_CD*_CE*_CG_CZ"' \
                     -oall {output} -tu ns
        '''

rule plot_t33b_y993_f1038_kde_clip_10ns:
    input:
        f1 = rules.make_t33b_f1038_xvgs_clip_10ns.output
    output:
        'results/{folder}/f1038_dists/F1038_Y993_kde_clip_10ns.png',
    params:
        xmin = 3.5,
        xmax = 12,
    script:
        '../scripts/plot_kde.py'
