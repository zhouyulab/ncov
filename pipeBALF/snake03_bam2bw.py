import os, sys

configfile: "config.yaml"
workdir: ".."

patients = config["SAMPLES"]["patient"]
ctrls = config["SAMPLES"]["ctrl"]
samples = patients + ctrls

uniq_reads = config["UNIQ_READS"]["patient"]
patient_scales = [1000000 / x for x in uniq_reads]
uniq_reads = config["UNIQ_READS"]["ctrl"]
ctrl_scales = [1000000 / x for x in uniq_reads]

colors = config["COLORS"]

sdirs = ["Forward", "Reverse"]
norms = ["raw", "norm2CPM"][0]

ruleorder: ctrl_norm > patient_norm

rule all:
    input:
        expand("result/track/raw/{patient}.bw", patient=patients),
        expand("result/track/raw/{ctrl}.{sdir}.bw", ctrl=ctrls, sdir=sdirs),
        # expand("result/track/norm2CPM/{patient}.bw", patient=patients),
        # expand("result/track/norm2CPM/{ctrl}.{sdir}.bw", ctrl=ctrls, sdir=sdirs),
        "result/track/raw/trackDb.ra",
        "result/track/norm2CPM/trackDb.ra",

rule patient_bam2bw:
    input:
        bam = "result/mapping/hg38/{patient}/{patient}_rmdup.bam",
        chr_info = config["chromsize"],
    output:
        wig = temp("result/track/raw/{patient}.wig"),
        bw = "result/track/raw/{patient}.bw",
    log: "result/track/raw/{patient}_bam2wig.log",
    params:
        wig_pfx = "result/track/raw/{patient}",
        strand = config["strand_rule"]["patient"],
    shell:"""
        set +u; source ~/miniconda3/bin/activate rseqc; set -u
        echo "input bam file: {input.bam}" > {log}
        bam2wig.py -i {input.bam} -s {input.chr_info} \
            -o {params.wig_pfx} -u >> {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule ctrl_bam2bw:
    input:
        bam = "result/mapping/hg38/{ctrl}/{ctrl}_rmdup.bam",
        chr_info = config["chromsize"],
    output:
        wig_f = temp("result/track/raw/{ctrl}.Forward.wig"),
        wig_r = temp("result/track/raw/{ctrl}.Reverse.wig"),
        bw_f = "result/track/raw/{ctrl}.Forward.bw",
        bw_r = "result/track/raw/{ctrl}.Reverse.bw",
    log: "result/track/raw/{ctrl}_bam2wig.log",
    params:
        wig_pfx = "result/track/raw/{ctrl}",
        strand = config["strand_rule"]["ctrl"],
    shell:"""
        set +u; source ~/miniconda3/bin/activate rseqc; set -u
        echo "input bam file: {input.bam}" > {log}
        bam2wig.py -i {input.bam} -s {input.chr_info} \
            -o {params.wig_pfx} -u --strand='{params.strand}' \
            >> {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule patient_norm:
    input:
        raw = "result/track/raw/{patient}.bw",
    output:
        norm = "result/track/norm2CPM/{patient}.bw",
    params:
        scale = lambda wildcards: patient_scales[patients.index(wildcards.patient)],
    shell: """
        script/bam2bw/norm_bw.py \
            -i {input.raw} -o {output.norm} -s {params.scale}
        """

rule ctrl_norm:
    input:
        raw = "result/track/raw/{ctrl}.{sdir}.bw",
    output:
        norm = "result/track/norm2CPM/{ctrl}.{sdir}.bw",
    params:
        scale = lambda wildcards: ctrl_scales[ctrls.index(wildcards.ctrl)],
    shell: """
        script/bam2bw/norm_bw.py \
            -i {input.raw} -o {output.norm} -s {params.scale}
        """

rule make_trackDB:
    ## LOCAL ##
    input:
        bw = "result/track/{norm}/{track}.bw",
    output:
        txt = temp("result/track/{norm}/{track}.txt"),
    params:
        url = "ftp://202.114.67.213/yangjiayi/balf_cov/rmdup/{norm}/bw/{track}.bw",
        color = lambda wildcards: colors[wildcards.track],
    run:
        with open(output.txt, 'w') as f:
            print("track " + wildcards.track, file=f)
            print("bigDataUrl " + params.url, file=f)
            print("shortLabel " + wildcards.track, file=f)
            print("longLabel " + wildcards.track, file=f)
            print("type bigWig", file=f)
            print("autoScale on", file=f)
            print("alwaysZero on", file=f)
            print("color " + params.color, file=f)
            print("maxHeightPixels 50:50:50", file=f)
            print("", file=f)

rule merge_raw_trackDB:
    ## LOCAL ##
    input:
        txt = expand("result/track/raw/{track}.txt", track=colors.keys()),
    output:
        db = "result/track/raw/trackDb.ra",
    shell: """
        cat {input.txt} > {output.db}
        """

rule merge_norm_trackDB:
    ## LOCAL ##
    input:
        txt = expand("result/track/norm2CPM/{track}.txt", track=colors.keys()),
    output:
        db = "result/track/norm2CPM/trackDb.ra",
    shell: """
        cat {input.txt} > {output.db}
        """