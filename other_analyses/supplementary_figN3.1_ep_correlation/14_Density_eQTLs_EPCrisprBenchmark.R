#!/usr/bin/env Rscript

options(width=170)
library(data.table)
library(ks)

band_width <- 10000
header = (0:500)*2000
grid_size = length(header)
max_distance = max(header)

for(group in c("positive","negative","match"))
{
    working_data = fread(file=sprintf("work/GTEx_ge_LCL.eQTLs.%s.data.gz",group), sep='\t', header=T, stringsAsFactors=FALSE)
    data_values  = as.numeric(working_data$distance)
    data_values  = data_values[data_values < max_distance]
    density      = kde(data_values, binned=TRUE, h = band_width, xmin=0, xmax=max_distance, compute.cont=TRUE, bgridsize=grid_size)
    correction   = sum(density$estimate)
    final        = abs(density$estimate/correction)
    together     = data.frame(distance = sprintf("%d",header),density = sprintf("%.7f",final), estimate = sprintf("%.4f",final*length(data_values)))
    fwrite(together, file = sprintf("work/GTEx_ge_LCL.eQTLs.%s.density",group) , sep="\t", quote=F)
    cat("Done",group,"density\n")
}

for(group in c("Pos","Neg","Match"))
{
    working_data = fread(file=sprintf("work/EPCrisprBenchmark_%s_glsCoefficient.data",group), sep='\t', header=T, stringsAsFactors=FALSE)
    data_values  = as.numeric(working_data$distance)
    data_values  = data_values[data_values < max_distance]
    density      = kde(data_values, binned=TRUE, h = band_width, xmin=0, xmax=max_distance, compute.cont=TRUE, bgridsize=grid_size)
    correction   = sum(density$estimate)
    final        = abs(density$estimate/correction)
    together     = data.frame(distance = sprintf("%d",header),density = sprintf("%.7f",final), estimate = sprintf("%.4f",final*length(data_values)))
    fwrite(together, file = sprintf("work/EPCrisprBenchmark_%s_glsCoefficient.density",group) , sep="\t", quote=F)
    cat("Done",group,"density\n")
}