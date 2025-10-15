#!/bin/bash
#SBATCH --job-name=my_project      # Job name
#SBATCH --output=output_%j.txt     # Output file (%j = job ID)
#SBATCH --error=error_%j.txt       # Error file
#SBATCH --time=2-00:00:00            # Time limit (hh:mm:ss)
#SBATCH --mem=400G                   # Memory per node
#SBATCH --cpus-per-task=2        # Number of CPU cores
#SBATCH --partition=engreitz,akundaje         # Partition name (ask your admin for the right one)


# OVERLAP ENHANCERS WITH ENHANCERS TO OBTAIN NUMBER OF BIOSAMPLES OVERLAPPING ENHANCERS
module load biology
module load bedtools 

bedtools intersect -a /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_ALL_active_enhancers.tsv -b /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_ALL_active_enhancers.tsv \
  -f 0.5 -wa -wb \
| awk 'BEGIN{OFS="\t"} $4!=$8' > overlap_enhancers.bed

# SPLIT FILES INTO CHROMSOMES
awk -F'\t' '{print > $1"_enhancers.tsv"}' /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/overlap_enhancers.bed

# GET NUMBER OF BIOSAMPLES 
awk '{
  key = $1"\t"$2"\t"$3"\t"$4
  biosample = $8
  print key"\t"biosample
}' chr${1}_enhancers.tsv | sort | uniq > chr${1}_temp.txt

awk -F'\t' '{count[$1"\t"$2"\t"$3"\t"$4]++; biosamples[$1"\t"$2"\t"$3"\t"$4][$5]=1}
END {
  for (k in count) {
    n=0
    for (b in biosamples[k]) { n++ }
    print k"\t"n
  }
}' chr${1}_temp.txt > chr${1}_biosample_count_enhancers.txt

# OBTAIN CONCATENATED NUMBER OF ACTIVE GENES/ENHANCERS ACROSS ALL BIOSAMPLES
#python concatenate_active_genes_enhancers.py
#touch output_${1}.txt
#zcat /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_${1}_active_genes_only.tsv.gz | gzip >> /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_ALL_active_genes.tsv.gz 
#zcat /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_${1}_active_enhancers_only.tsv.gz | gzip >> /oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_ALL_active_enhancers.tsv.gz

# THESE SCRIPTS HELP TO GENERATE EXTENDED FIG 6 G&H
#python generate_encodee2g_fig_g_h.py $1 

# THESE SCRIPTS HELP TO GENERATE E2G PROPERTIES ACROSS DIFFERENT FILE FORMATS/PREDICTORS
#generate_encodee2g_e2g_properties.py $1 #generate_graphreg_e2g_properties.py #generate_e2g_properties.py $1 
