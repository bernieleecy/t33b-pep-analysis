# calculate unique hbonds using VMD


rule get_hbonds:
    input:
        xtc="runs/{folder}/combined.xtc",
        gro="runs/{folder}/protein_noh2o.gro",
    output:
        hbonds_time="results/{folder}/vmd_hbonds/data/hbonds_v_time.dat",
        pairs="results/{folder}/vmd_hbonds/data/hbonds_pairs.dat",
        unique="results/{folder}/vmd_hbonds/data/hbonds_unique.dat",
        tcl_pairs="results/{folder}/vmd_hbonds/hbonds_pairs.tcl",
        tcl_unique="results/{folder}/vmd_hbonds/hbonds_unique.tcl",
    params:
        sel_1="resid 883 to 1071",
        sel_2="resid 1 to 20",
        vmd="/Applications/VMD\ 1.9.4a51-arm64-Rev9.app/Contents/MacOS/startup.command",
        dist=3.5,
        angstrom=30,
    shell:
        """
        echo -e 'package require hbonds \
                \nhbonds -sel1 [atomselect top "{params.sel_1}"] -sel2 [atomselect top "{params.sel_2}"] \
                -dist {params.dist} -ang {params.angstrom} \
                -writefile yes -type pair \
                -outfile {output.hbonds_time} -detailout {output.pairs}' > {output.tcl_pairs}
        echo -e 'package require hbonds \
                \nhbonds -sel1 [atomselect top "{params.sel_1}"] -sel2 [atomselect top "{params.sel_2}"] \
                -dist {params.dist} -ang {params.angstrom} \
                -writefile yes -type unique \
                -outfile {output.hbonds_time} -detailout {output.unique}' > {output.tcl_unique}
        {params.vmd} {input.xtc} {input.gro} -dispdev text < {output.tcl_pairs}
        {params.vmd} {input.xtc} {input.gro} -dispdev text < {output.tcl_unique}
        """


rule get_key_bd_pairs:
    """
    Look for hydrogen bonds between the bromodomain and peptide (above 5% occupancy)
    """
    input:
        "results/{folder}/vmd_hbonds/data/hbonds_{type}.dat",
    output:
        "results/{folder}/vmd_hbonds/data/bd_hbonds_above_5_{type}.dat",
    shell:
        """
        cat {input} | tr -d '%' |
        awk ' /[[:alpha:]]{{3}}[1-2][0-9]/ {{ if ($3>5) print }} ' |
        sort -nk3 -r > {output}
        """


rule get_key_phd_pairs:
    """
    Look for hydrogen bonds between the PHD finger and peptide (above 5% occupancy)
    """
    input:
        "results/{folder}/vmd_hbonds/data/hbonds_{type}.dat",
    output:
        "results/{folder}/vmd_hbonds/data/phd_hbonds_above_5_{type}.dat",
    shell:
        """
        cat {input} | tr -d '%' |
        awk ' /[[:alpha:]]{{1}}[[:alnum:]]{{1}}[[:alpha:]]{{1}}[1-9][0-2]?[\-]/ {{ if ($3>5) print }} ' |
        sort -nk3 -r > {output}
        """
