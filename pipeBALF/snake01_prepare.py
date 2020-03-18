import os

configfile: "config.yaml"
workdir: ".."

patients = config["SAMPLES"]["patient"]
ctrls = config["SAMPLES"]["ctrl"]
samples = patients + ctrls
reads = ["R1", "R2"]

rule all:
    input:
        expand("data/{ctrl}_{read}.fq.gz", ctrl=ctrls[1], read=reads),
        expand("result/fastqc/{sample}_{read}/{sample}_{read}_fastqc.zip", 
               sample=samples, read=reads),
        expand("result/cutadapt_fq/{sample}/{sample}_{read}_cutadapt.fq.gz", 
               sample=samples, read=reads),
        expand("result/fastqc/{sample}_{read}/{sample}_{read}_cutadapt_fastqc.zip", 
               sample=samples, read=reads),

rule srr2fq:
    input:
        srr = "data/control/{ctrl}",
    output:
        raw1 = temp("data/control/{ctrl}_split/{ctrl}_1.fastq"),
        raw2 = temp("data/control/{ctrl}_split/{ctrl}_2.fastq"),
        gz1 = "data/{ctrl}_R1.fq.gz",
        gz2 = "data/{ctrl}_R2.fq.gz",
    params:
        outdir = "data/control/{ctrl}_split",
    threads: 7
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        fasterq-dump {input.srr} --split-3 -e {threads} -O {params.outdir}
        pigz -p {threads} -c {output.raw1} > {output.gz1}
        pigz -p {threads} -c {output.raw2} > {output.gz2}
        set +u; conda deactivate; set -u
        """

rule fastqc_raw:
    input:
        fq = "data/{sample}_{read}.fq.gz",
    output:
        qc = "result/fastqc/{sample}_{read}/{sample}_{read}_fastqc.zip",
    params:
        qc_dir = "result/fastqc/{sample}_{read}"
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p {params.qc_dir}
        fastqc -o {params.qc_dir} {input.fq}
        set +u; conda deactivate; set -u
        """

rule cutadapt:
    input:
        fq1 = "data/{sample}_R1.fq.gz",
        fq2 = "data/{sample}_R2.fq.gz",
    output:
        cut1 = "result/cutadapt_fq/{sample}/{sample}_R1_cutadapt.fq.gz",
        cut2 = "result/cutadapt_fq/{sample}/{sample}_R2_cutadapt.fq.gz",
    log: "result/cutadapt_fq/{sample}/{sample}_cutadapt.log",
    params:
        adp = "AGATCGGAAGAG",
        m = 20,
        maxn = 5,
        q = 30,
    threads: 4
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        cutadapt -a {params.adp} -A {params.adp} \
            -o {output.cut1} -p {output.cut2} \
            {input.fq1} {input.fq2} \
            -m {params.m} --max-n {params.maxn} \
            -q {params.q} -j {threads} > {log}
        set +u; conda deactivate; set -u
        """

rule fastqc_aft:
    input:
        fq = "result/cutadapt_fq/{sample}/{sample}_{rep}_cutadapt.fq.gz",
    output:
        qc = "result/fastqc/{sample}_{rep}/{sample}_{rep}_cutadapt_fastqc.zip",
    params:
        qc_dir = "result/fastqc/{sample}_{rep}"
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p {params.qc_dir}
        fastqc -o {params.qc_dir} {input.fq}
        set +u; conda deactivate; set -u
        """
    
