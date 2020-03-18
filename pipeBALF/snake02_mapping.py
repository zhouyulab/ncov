import os

configfile: "config.yaml"
workdir: ".."

patients = config["SAMPLES"]["patient"]
ctrls = config["SAMPLES"]["ctrl"]
samples = patients + ctrls
genos = ["chimeric", "hg38"]

rule all:
    input:
        expand("result/mapping/phix/{patient}/{patient}_Aligned.out.bam", patient=patients),
        expand("result/mapping/hg38/{patient}/{patient}_Aligned.out.bam", patient=patients),
        expand("result/mapping/ncov/{patient}/{patient}_Aligned.out.bam", patient=samples),
        expand("result/mapping/ctrl_hg38/{ctrl}/{ctrl}_Aligned.out.bam", ctrl=ctrls),
        expand("result/mapping/chimeric/{patient}/{patient}_Aligned.out.bam", patient=patients),
        "result/mapping/hg38/link_finish.log",
        expand("result/check_ss/{sample}_strand_specificity.log", sample=samples),
        expand("result/mapping/{geno}/{sample}/{sample}_uniq_sort.bam.bai", geno=genos, sample=samples),
        expand("result/mapping/{geno}/{sample}/{sample}_rmdup.bam", geno=genos, sample=samples),

rule phix_mapping:
    input:
        index_dir = config["star_index"]["phix"],
        fq1 = "result/cutadapt_fq/{patient}/{patient}_R1_cutadapt.fq.gz",
        fq2 = "result/cutadapt_fq/{patient}/{patient}_R2_cutadapt.fq.gz",
    output:
        bam = "result/mapping/phix/{patient}/{patient}_Aligned.out.bam",
        un1 = temp("result/mapping/phix/{patient}/{patient}_Unmapped.out.mate1"),
        un2 = temp("result/mapping/phix/{patient}/{patient}_Unmapped.out.mate2"),
    params:
        bam_pfx = "result/mapping/phix/{patient}/{patient}_",
    threads: 7
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        STAR --runThreadN {threads} \
            --genomeDir {input.index_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.bam_pfx} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx
        set +u; conda deactivate; set -u
        """

rule chimeric_mapping:
    input:
        index_dir = config["star_index"]["chimeric"],
        fq1 = rules.phix_mapping.output.un1,
        fq2 = rules.phix_mapping.output.un2,
    output:
        bam = "result/mapping/chimeric/{patient}/{patient}_Aligned.out.bam",
        un1 = "result/mapping/chimeric/{patient}/{patient}_Unmapped.out.mate1",
        un2 = "result/mapping/chimeric/{patient}/{patient}_Unmapped.out.mate2",
    params:
        bam_pfx = "result/mapping/chimeric/{patient}/{patient}_",
    threads: 7
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        STAR --runThreadN {threads} \
            --genomeDir {input.index_dir} \
            --chimSegmentMin 10 \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.bam_pfx} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx
        set +u; conda deactivate; set -u
        """

rule hg38_mapping:
    input:
        index_dir = config["star_index"]["hg38"],
        gtf = config["gtf"],
        fq1 = rules.phix_mapping.output.un1,
        fq2 = rules.phix_mapping.output.un2,
    output:
        bam = "result/mapping/hg38/{patient}/{patient}_Aligned.out.bam",
        un1 = temp("result/mapping/hg38/{patient}/{patient}_Unmapped.out.mate1"),
        un2 = temp("result/mapping/hg38/{patient}/{patient}_Unmapped.out.mate2"),
        gz1 = "result/mapping/phix/{patient}/{patient}_Unmapped.out.mate1.gz",
        gz2 = "result/mapping/phix/{patient}/{patient}_Unmapped.out.mate2.gz",
    params:
        bam_pfx = "result/mapping/hg38/{patient}/{patient}_",
    threads: 7
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        STAR --runThreadN {threads} \
            --genomeDir {input.index_dir} \
            --sjdbGTFfile {input.gtf} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.bam_pfx} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx
        pigz -p {threads} -c {input.fq1} > {output.gz1}
        pigz -p {threads} -c {input.fq2} > {output.gz2}
        set +u; conda deactivate; set -u
        """

rule whu01_mapping:
    input:
        index_dir = config["star_index"]["ncov"],
        fq1 = rules.hg38_mapping.output.un1,
        fq2 = rules.hg38_mapping.output.un2,
    output:
        bam = "result/mapping/ncov/{patient}/{patient}_Aligned.out.bam",
        un1 = "result/mapping/ncov/{patient}/{patient}_Unmapped.out.mate1",
        un2 = "result/mapping/ncov/{patient}/{patient}_Unmapped.out.mate2",
        gz1 = "result/mapping/hg38/{patient}/{patient}_Unmapped.out.mate1.gz",
        gz2 = "result/mapping/hg38/{patient}/{patient}_Unmapped.out.mate2.gz",
    params:
        bam_pfx = "result/mapping/ncov/{patient}/{patient}_",
    threads: 1
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        STAR --runThreadN {threads} \
            --genomeDir {input.index_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.bam_pfx} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx
        pigz -p {threads} -c {input.fq1} > {output.gz1}
        pigz -p {threads} -c {input.fq2} > {output.gz2}
        set +u; conda deactivate; set -u
        """

rule ctrl_mapping:
    input:
        index_dir = config["star_index"]["hg38"],
        gtf = config["gtf"],
        fq1 = "result/cutadapt_fq/{ctrl}/{ctrl}_R1_cutadapt.fq.gz",
        fq2 = "result/cutadapt_fq/{ctrl}/{ctrl}_R2_cutadapt.fq.gz",
    output:
        bam = "result/mapping/ctrl_hg38/{ctrl}/{ctrl}_Aligned.out.bam",
        un1 = "result/mapping/ctrl_hg38/{ctrl}/{ctrl}_Unmapped.out.mate1",
        un2 = "result/mapping/ctrl_hg38/{ctrl}/{ctrl}_Unmapped.out.mate2",
    params:
        bam_pfx = "result/mapping/ctrl_hg38/{ctrl}/{ctrl}_",
    threads: 7
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        STAR --runThreadN {threads} \
            --genomeDir {input.index_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --sjdbGTFfile {input.gtf} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.bam_pfx} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx
        set +u; conda deactivate; set -u
        """

rule make_link:
    ## LOCAL ##
    input:
    output:
        "result/mapping/hg38/link_finish.log",
    params:
        "result/mapping/hg38",
    shell: """
        cd {params}
        for f in `ls ../ctrl_hg38`; do
            ln -s ../ctrl_hg38/${{f}}
        done
        touch link_finish.log
        """

#result: ss:""
rule check_ss:
    input:
        bed = config["bed"],
        bam = "result/mapping/hg38/{sample}/{sample}_Aligned.out.bam",
    output:
        res = "result/check_ss/{sample}_strand_specificity.log",
    shell:"""
        set +u; source ~/miniconda3/bin/activate rseqc; set -u
        echo "sample name:{wildcards.sample}" > {output.res}
        infer_experiment.py -r {input.bed} -i {input.bam} \
            >> {output.res} 2>&1
        set +u; conda deactivate; set -u
        """
    
rule index_bam:
    input:
        bam = rules.check_ss.input.bam,
    output:
        uniq = "result/mapping/{geno}/{sample}/{sample}_uniq_sort.bam",
        bai = "result/mapping/{geno}/{sample}/{sample}_uniq_sort.bam.bai",
    threads: 4
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        samtools view -@ {threads} -hb -q 255 {input.bam} | \
        samtools sort -@ {threads} -o {output.uniq}
        samtools index -@ {threads} {output.uniq}
        set +u; conda deactivate; set -u
        """

rule picard_rmdup:
    input:
        bam = rules.index_bam.output.uniq,
    output:
        bam = "result/mapping/{geno}/{sample}/{sample}_rmdup.bam",
        bai = "result/mapping/{geno}/{sample}/{sample}_rmdup.bam.bai",
        mark = "result/mapping/{geno}/{sample}/{sample}_rmdup_markdup.txt",
    log: "result/mapping/{geno}/{sample}/{sample}_rmdup.log",
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.mark} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true 2> {log}
        samtools index -@ 1 {output.bam}
        set +u; conda deactivate; set -u
        """

