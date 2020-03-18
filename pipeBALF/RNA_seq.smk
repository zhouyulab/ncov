import os
from collections import defaultdict
include: "MetaInfo.smk"

rule sort_nCoV_bam:
    input:
        bam = "data/RNA_seq_rm_dup/nCoV/{donor}_{rep}_rmdup.bam"
    output:
        bam = os.path.join("analysis", "bam", "nCoV", "{donor}.{rep}.sort.bam"),
    params:
        job_name = "sort_nCoV_bam",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=1",
        log = "log",
        queue = "wo_dell01",
    shell:
        """
samtools sort -@ 24 -o {output.bam} {input.bam}
samtools index {output.bam}
        """

rule sort_Ctrl_bam:
    input:
        bam = "data/RNA_seq_rm_dup/Ctrl/{SRR}_rmdup.bam"
    output:
        bam = os.path.join("analysis", "bam", "Ctrl", "{SRR}.sort.bam"),
    params:
        job_name = "sort_Ctrl_bam",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=1",
        log = "log",
        queue = "wo_dell01",
    shell:
        """
samtools sort -@ 24 -o {output.bam} {input.bam}
samtools index {output.bam}
        """

rule feature_counts_nCoV_worker:
    input:
        bam = rules.sort_nCoV_bam.output.bam,
        gtf = REF_GTF,
    output:
        res = os.path.join("analysis", "feature_count", "nCoV", "{donor}.{rep}.featurecount"),
    params:
        job_name = "feature_counts_nCoV_worker",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=14",
        log = "log",
        queue = "All",
    priority: 100
    shell:
        """
featureCounts -M -T 28 -a {input.gtf} -o {output}.fc.out {input.bam}
cat {output}.fc.out | sed -n '3,$p' | awk '{{FS=OFS="\t"}}{{print $1, $7}}' > {output}
rm {output}.fc.out
        """

rule feature_counts_Ctrl_worker:
    input:
        bam = rules.sort_Ctrl_bam.output.bam,
        gtf = REF_GTF,
    output: 
        res = os.path.join("analysis", "feature_count", "Ctrl", "{SRR}.featurecount"),
    params: 
        job_name = "feature_counts_Ctrl_worker",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=14",
        log = "log",
        queue = "All",
    priority: 100
    shell:
        """
featureCounts -M -T 28 -a {input.gtf} -o {output}.fc.out {input.bam}
cat {output}.fc.out | sed -n '3,$p' | awk '{{FS=OFS="\t"}}{{print $1, $7}}' > {output}
rm {output}.fc.out
        """


rule merge_featurecounts:
    input:
        cov_flag = expand(rules.feature_counts_nCoV_worker.output, donor=CoV_ID, rep=CoV_rep),
        ctrl_flag = expand(rules.feature_counts_Ctrl_worker.output, SRR=Ctrl_ID),
    output:
        merge = os.path.join("analysis", "featureCounts.merge.tsv"),
    params:
        job_name = "merge_featurecounts",
        walltime = "480:00:00",
        nodes = "nodes=1:ppn=1",
        log = "log",
        queue = "all",
        merge_featurecount = "Rscript scripts/merge_featurecount.R",
    run:
        name1_li = list()
        name2_li = list()
        file_li = list()
        for donor in CoV_ID:
            for rep in CoV_rep:
                name1_li.append(donor)
                name2_li.append(rep)
                file_li.append(rules.feature_counts_nCoV_worker.output.res.format(donor=donor, rep=rep))
                
        for SRR in Ctrl_ID:
            name1_li.append("Ctrl")
            name2_li.append(SRR)
            file_li.append(rules.feature_counts_Ctrl_worker.output.res.format(SRR=SRR))
            
        cmd = "{prog} -i {file_li} --name1 {name1_li} --name2 {name2_li} -o {output}".format(prog=params.merge_featurecount, file_li=" ".join(file_li), name1_li=" ".join(name1_li), name2_li=" ".join(name2_li), output=output.merge)
        shell(cmd)
