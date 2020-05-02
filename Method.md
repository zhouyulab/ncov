# Description of DEG analysis

Differentially expressed genes were called using DESeq2 package (v1.26.0). 
Considering that BALFs and PBMCs samples have different sequencing depths, 
we used different empirical parameters for these two data sets. 

For BALF data, we used FC (Fold Change) > 4, padj (adjusted P-value) < 1e-10, and average reads counts in all samples > 10 as the criteria for identifying differentially expressed genes, due to low sequencing depth in this set. 
And for PBMC data with higher sequencing depth, the following thresholds are used: FC = 2, padj = 0.01, and average reads counts in all samples = 100. 

