#!/usr/bin/env Rscript

options(width=170)
options(scipen=999)
library(data.table)
library(scam)

data_matrix = fread(file="work/glsCoefficient.statistics", sep='\t', header=T, stringsAsFactors=FALSE)

distance  = data_matrix$distance
test_data = data.frame(distance = distance)
weights   = c(rep(100,10),rep(50,30),rep(10,60),rep(1,length(distance)-100))

#########################
value = data_matrix$mean
data_table = data.frame(distance = distance, value = value)
cat("Working on mean line smooth\n")
con <- scam(value ~ s( distance, k = 400, bs="mpd"), data = data_table, weights = weights)
predict_data = predict(con,test_data)
smooth_mean = sprintf("%.6f",predict_data)
cat("Done\n")

#########################
value = data_matrix$lower
data_table = data.frame(distance  = distance, value = value)

cat("Working on lower 2.5% line smooth\n")
con <- scam(value ~ s( distance, k = 400, bs="mpd"), data = data_table, weights = weights)
predict_data = predict(con,test_data)
smooth_lower = sprintf("%.6f",predict_data)
cat("Done\n")

#########################
value = data_matrix$higher
data_table = data.frame(distance  = distance, value = value)

cat("Working on higher 97.5% line smooth\n")
con <- scam(value ~ s( distance, k = 400, bs="mpd"), data = data_table, weights = weights)
predict_data = predict(con,test_data)
smooth_higher = sprintf("%.6f",predict_data)
cat("Done\n")

report_data = data.table(distance = sprintf("%.0f",distance), mean = smooth_mean, lower = smooth_lower, higher = smooth_higher)
fwrite(report_data, file = "work/Figure_S3A.smooth", sep="\t", quote=F)

proc.time()
