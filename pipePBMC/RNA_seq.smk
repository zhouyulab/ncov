import os
include: "MetaInfo.smk"

RNA_BASE = os.path.join("analysis", "RNA_seq")

rule RNA_rm_rRNA:
    input:
        fq_R1 = os.path.join("data", "RNA_seq", "{sample}.R1.fq.gz"),
        fq_R2 = os.path.join("data", "RNA_seq", "{sample}.R2.fq.gz"),
    output:
        unmap_R1 = os.path.join(RNA_BASE, "rm_rRNA", "{sample}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_BASE, "rm_rRNA", "{sample}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "rm_rRNA", "{sample}"),
        indx = rRNA_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --outReadsUnmapped Fastx \
    --readFilesCommand gunzip -c
gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
rm -r {params.mapping_dir}
        """

rule RNA_STAR:
    input:
        fq_R1 = rules.RNA_rm_rRNA.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA.output.unmap_R2,
        gtf = REF_GTF,
    output:
        bam = os.path.join(RNA_BASE, "STAR", "{sample}.sort.bam"),
        unmap_R1 = os.path.join(RNA_BASE, "STAR", "{sample}.unmap.R1.fq.gz"),
        unmap_R2 = os.path.join(RNA_BASE, "STAR", "{sample}.unmap.R2.fq.gz"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "STAR", "{sample}"),
        indx = hg38_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --quantMode TranscriptomeSAM \
    --genomeDir {params.indx} \
    --sjdbGTFfile {input.gtf} \
    --sjdbScore 1 \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outReadsUnmapped Fastx \
    --readFilesCommand gunzip -c
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
gzip -c {params.mapping_dir}/Unmapped.out.mate1 > {output.unmap_R1}
gzip -c {params.mapping_dir}/Unmapped.out.mate2 > {output.unmap_R2}
        """

rule RNA_STAR_nCoV:
    input:
        fq_R1 = rules.RNA_STAR.output.unmap_R1,
        fq_R2 = rules.RNA_STAR.output.unmap_R2,
    output:
        bam = os.path.join(RNA_BASE, "STAR_nCoV", "{sample}.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "STAR_nCoV", "{sample}"),
        indx = WHU01_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --readFilesCommand gunzip -c
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
        """

rule RNA_STAR_nCoV_direct:
    input:
        fq_R1 = rules.RNA_rm_rRNA.output.unmap_R1,
        fq_R2 = rules.RNA_rm_rRNA.output.unmap_R2,
    output:
        bam = os.path.join(RNA_BASE, "STAR_nCoV_direct", "{sample}.sort.bam"),
    params:
        mapping_dir = os.path.join(RNA_BASE, "STAR_nCoV_direct", "{sample}"),
        indx = WHU01_INDX,
    threads: 32
    shell:
        """
if [[ -e {params.mapping_dir} ]]; then
    rm -r {params.mapping_dir}
fi
mkdir -p {params.mapping_dir}
STAR --runThreadN {threads} \
    --genomeDir {params.indx} \
    --outFileNamePrefix {params.mapping_dir}/ \
    --readFilesIn {input.fq_R1} {input.fq_R2} \
    --readFilesCommand gunzip -c
samtools sort -@ {threads} -o {output.bam} {params.mapping_dir}/Aligned.out.sam
samtools index {output.bam}
rm {params.mapping_dir}/Aligned.out.sam
        """

rule RNA_rmDup:
    input:
        bam = rules.RNA_STAR.output.bam,
    output:
        bam = os.path.join(RNA_BASE, "rmDup", "{sample}.bam"),
        metrics = os.path.join(RNA_BASE, "rmDup", "{sample}.metrics")
    shell:
        """
picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics}
samtools index {output.bam}
        """

rule RNA_bam2bw:
    input:
        bam = rules.RNA_rmDup.output.bam,
        size = CHROM_SIZE,
    output:
        flag = touch(os.path.join(RNA_BASE, "bw", "{sample}.flag")),
    params:
        job_name = "RNA_bam2bw",
        total_sum = 1000000000,
        work_dir = os.path.join(RNA_BASE, "bw"),
        rule = "1+-,1-+,2++,2--",
    shell:
        """
bam2wig.py -t {params.total_sum} -s {input.size} -i {input.bam} -o {params.work_dir}/{wildcards.sample} -d "{params.rule}"
wigToBigWig -clip {params.work_dir}/{wildcards.sample}.Forward.wig {input.size} {params.work_dir}/{wildcards.sample}.Forward.bw
wigToBigWig -clip {params.work_dir}/{wildcards.sample}.Reverse.wig {input.size} {params.work_dir}/{wildcards.sample}.Reverse.bw
rm {params.work_dir}/{wildcards.sample}.Forward.wig {params.work_dir}/{wildcards.sample}.Reverse.wig
        """

rule RNA_bam2bw_test:
    input:
        expand(rules.RNA_bam2bw.output, sample=ALL_SAMPLEs),

rule RNA_seq:
    input:
        expand(rules.RNA_rmDup.output, sample=ALL_SAMPLEs),
        expand(rules.RNA_STAR_nCoV.output, sample=ALL_SAMPLEs),
        expand(rules.RNA_bam2bw.output, sample=ALL_SAMPLEs),
        expand(rules.RNA_STAR_nCoV_direct.output, sample=ALL_SAMPLEs),


rule RNA_STAR_nCoV_direct_test:
    input:
        expand(rules.RNA_STAR_nCoV_direct.output, sample=ALL_SAMPLEs),
