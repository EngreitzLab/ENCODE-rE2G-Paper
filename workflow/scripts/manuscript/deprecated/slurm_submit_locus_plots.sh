#!/usr/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=128G
#SBATCH --output=supplementary_fig_locus_plots.log
#SBATCH --partition=owners

# load R module
ml R/4.0.2

# run R code
Rscript supplementary_fig_locus_plots.R

