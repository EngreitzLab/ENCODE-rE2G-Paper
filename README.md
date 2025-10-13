# Analyses for the ENCODE distal regulation companion paper
This repository contains analysis code for analyses presented in the ENCODE distal regulation
companion paper, as well as links to the used benchmarking pipelines.

## Setting up the repository
To obtain a local copy of the repository, including the `crispr_analyses/CRISPR_comparison`
submodule containing the CRISPR benchmarking pipeline version used for analyses in this paper, open
a terminal and run:
```console
git clone --recurse-submodules git@github.com:EngreitzLab/ENCODE-Distal-Regulation-Paper.git
```

## Analysis code
The following subdirectories contain code used to perform analyses shown in the paper:

### crispr_analyses:
Snakemake workflow for CRISPR benchmarking and other related analyses results shown in Main
Figures 1, 2, 4 and 5, Extended Data Figures 1, 2, 3, 7, Supplementary Figures S1.1, S1.2, S7, S8,
S9, S10, as well as Supplementary Tables S2, S9, S10, S14, S15, S16.

### gwas_analyses:
Code used to compute GWAS metrics shown in analyses of the manuscript. GWAS benchmarks were
performed using the GWAS benchmarking pipeline (see below).

### process_dnase_metadata:
Code to process DNase-seq metadata from the ENCODE portal to select input files to generate
ENCODE-rE2G predictions for 1,458 human DNase-seq experiments.

### other analyses:
This subdirectory contains code used in other analyses across the manuscript, including annotating
CRISPR data with chromatin categories, analyzing enhancer synergy and analyses related to 3D
contact, enhancer-gene correlation and promoter classes.


## Benchmarking pipelines
Following pipelines were used to benchmark enhancer-gene predictions against CRISPR, eQTL and GWAS
data:

### CRISPR benchmarks
A copy of the CRISPR benchmarking used to intersect predictions with CRISPR enhancer perturbation
results can be found in the **crispr_analyses** subdirectory or in the original GitHib repository
for general use: https://github.com/EngreitzLab/CRISPR_comparison.

### eQTL benchmarks
All eQTL benchmarking analyses, including figures shown in the paper, were computed using the eQTL
benchmarking pipeline available here: https://github.com/EngreitzLab/eQTLEnrichment.

### GWAS benchmarks
The GWAS benchmarks were performed using the GWAS benchmarking pipeline available here: 
https://github.com/Deylab999MSKCC/e2g-benchmarking. Additional analysis code for GWAS-related
analyses is available in the *gwas_analyses* subdirectory.
