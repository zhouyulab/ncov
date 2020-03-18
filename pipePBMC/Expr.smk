include: "RNA_seq.smk"

EXPR_BASE = os.path.join("analysis", "expr")

rule RNA_feature_counts:
    input:
        bam = rules.RNA_STAR.output.bam,
        gtf = REF_GTF,
    output:
        res = os.path.join(EXPR_BASE, "feature_count", "{sample}.featurecount"),
    threads: 32
    shell:
        """
featureCounts -M -T {threads} -a {input.gtf} -o {output.res}.fc.out {input.bam}
cat {output.res}.fc.out | sed -n '3,$p' | awk '{{FS=OFS="\t"}}{{print $1, $7}}' > {output.res}
rm {output.res}.fc.out
        """

rule merge_feature_counts:
    input:
        flag = expand(rules.RNA_feature_counts.output, sample=ALL_SAMPLEs),
    output:
        merge = os.path.join(EXPR_BASE, "featureCounts.merge.tsv"),
    params:
        merge_featurecount = "Rscript scripts/merge_featurecount.R",
    run:
        cmd = "{prog} -i {file_li} --name {name_li} -o {output}".format(
    prog=params.merge_featurecount, 
    file_li=" ".join([rules.RNA_feature_counts.output.res.format(sample=sample) for sample in ALL_SAMPLEs]),
    name_li=" ".join(ALL_SAMPLEs),
    output=output.merge
    )
        shell(cmd)


rule expr:
    input:
        rules.merge_feature_counts.output,
