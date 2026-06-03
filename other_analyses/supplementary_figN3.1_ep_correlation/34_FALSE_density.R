#!/usr/bin/env Rscript

options(width=170)
library(data.table)
library(ks)

band_width <- 10000
header = (0:500)*2000
grid_size = length(header)


for(label in c("broad_broad","broad_specific","variable_broad","variable_specific"))
{
    working_data = fread(file=sprintf("work/Classes.%s.data",label), sep='\t', header=T, stringsAsFactors=FALSE)
    cat(sprintf("within %s:\t%d pos %d neg\n",label, nrow(working_data[Regulated == "TRUE"]),nrow(working_data[Regulated == "FALSE"])))
    negative_data = working_data[Regulated == "FALSE"]
    data_values   = as.numeric(negative_data$distance)

    cat("Mean ",label," distance:", mean(data_values),"\n")
    cat("Max  ",label," distance:", max(data_values),"\n")
    cat("Min  ",label," distance:", min(data_values),"\n\n")

    density     = kde(data_values, binned=TRUE, h = band_width, xmin=0, xmax=1000000, compute.cont=TRUE, bgridsize=grid_size)
    correction  = sum(density$estimate)
    final       = abs(density$estimate/correction)
    together    = data.frame(distance = sprintf("%d",header),density = sprintf("%.7f",final), estimate = sprintf("%.4f",final*length(data_values)))

    fwrite(together, file = sprintf("work/Classes.%s.density",label), sep="\t", quote=F)
}
