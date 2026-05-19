#!/usr/bin/env Rscript

library(data.table)
set.seed(127)
options(scipen = 999)

cores = 2
report_file = "work/glsCoefficient.statistics"

work_data = fread(file = "work/glsCoefficient.data.gz", sep='\t', header=F, stringsAsFactors=FALSE, blank.lines.skip=T, tmpdir="./temp", nThread = cores)
colnames(work_data) = c("distance","glsCoefficient")
cat("Load work/glsCoefficient.data.gz\n")

write(paste( c("distance","mean","lower","higher"), collapse="\t"), file = report_file)

for(index in 0:3000000)
{
    start = 1+index*5000
    end   = 1+index*5000+9999
    if(end > nrow(work_data))
	break

    local_data     = work_data[start:end]
    distance       = as.numeric(local_data$distance)
    glsCoefficient = as.numeric(local_data$glsCoefficient)

    mean_glsCoefficient = mean(glsCoefficient)
    mean_distance       = mean(distance)
    analysis = quantile(glsCoefficient,probs = c(0.025, 0.5, 0.975))
    write(paste(c(mean_distance,mean_glsCoefficient,analysis[1],analysis[3]), collapse="\t"), file = report_file, append = T)
    if((start-1)%%1000000 == 0)
	cat(sprintf("Done %8d .. %8d\n",start,end))
}
cat("Done\n")