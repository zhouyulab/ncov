#!/usr/bin/env Rscript
# Install requried R and BioConductor packages

packages <- c(
  "ggplot2", 
  "dplyr",
  "readr", 
  "circlize",
  "readr",
  "gridExtra"
)
cat("Packages to install: ", setdiff(packages, rownames(installed.packages())), "\n")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), dep=TRUE)
}

biopkgs <- c(
  "DESeq2",
  "DelayedArray",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
)

pkg2 <- setdiff(biopkgs, rownames(installed.packages()))
if (length(pkg2) > 0) {
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
}
cat("Packages to install: ", pkg2, "\n")
if (length(pkg2) > 0) {
  BiocManager::install(pkg2)
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("clusterProfiler")

