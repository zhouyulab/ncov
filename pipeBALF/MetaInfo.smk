import os
from pybicl.io import iterline
curDir = os.getcwd()

CoV_ID = ["patient1", "patient2"]
CoV_rep = ["rep1", "rep2"]
Ctrl_ID = ["SRR10571724", "SRR10571730", "SRR10571732"]


#################
##             ##
## Genome info ##
##             ##
#################

COMMON_DIR = "/home/wangdehe/common"
GENOME_FA = os.path.join(COMMON_DIR, "UCSC", "genome", "fasta", "hg38.fa")
GENOME_2bit = os.path.join(COMMON_DIR, "UCSC", "genome", "2bit", "hg38.2bit")
CHROM_SIES = os.path.join(COMMON_DIR, "UCSC", "genome", "size", "hg38.sizes")
REF_GFF = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "hg38_gencode.v32.gff3")
REF_GTF = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "hg38_gencode.v32.gtf")
REF_BED = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "hg38_gencode.v32.bed12")
REF_GP = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "hg38_gencode.v32.gp")
REF_GENE_BED = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "hg38_gencode.v32.gene.bed6")
REF_VCF = os.path.join(COMMON_DIR, "UCSC", "ref_gene", "hg38", "common_all.vcf.gz")
STAR_rRNA_INDX = os.path.join(COMMON_DIR, "rRNA", "star_indx", "human")
STAR_INDX = os.path.join(COMMON_DIR, "index", "STAR", "hg38")
BWA_INDX = os.path.join(COMMON_DIR, "index", "bwa", "hg38.fa")
