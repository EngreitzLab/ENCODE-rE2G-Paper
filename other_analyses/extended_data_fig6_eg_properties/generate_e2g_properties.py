import pandas as pd
import numpy as np
import os, sys

# add biosample
biosample = sys.argv[1]

#biosample="encode_e2g_predictions_K562"
# for single file
#data = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/extended/thresholded_predictions/encode_e2g_predictions_K562_thresholded_predictions.tsv.gz", sep="\t")
data = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ABC/dhs_only/{}/Predictions/abc_predictions_{}_thresholded_predictions.tsv.gz".format(biosample, biosample), sep="\t")
subset = data.loc[data['class']!="promoter"]

# get midpoint 
subset['midpoint'] = 0.5*(subset['start']+subset['end'])
subset['midpoint'] = subset['midpoint'].astype('int')

# log10(Distance to TSS)
subset['distance'] = subset['TargetGeneTSS'] - subset['midpoint']
subset['log10(Distance_to_TSS)'] = np.log10(subset['distance'])
distance_to_tss = np.array(subset['log10(Distance_to_TSS)']+1)
np.savez("ABC/{}_distance_to_tss.npz".format(biosample), distance_to_tss)

# Number of enhancer-gene links
#subset['Num_enhancer_gene_links'] = len(subset)
#
#with open('ENCODE-rE2G_extended/{}_num_enhancer_gene_links.txt'.format(biosample), 'w') as output:
#    output.write(str(subset['Num_enhancer_gene_links'].values[0]))
#
## Number of genes with at least 1 distal enhancer
#data_genes = len(subset.groupby(['TargetGene']))
#with open('ENCODE-rE2G_extended/{}_num_genes_one_distal_enhancer.txt'.format(biosample), 'w') as output:
#    output.write(str(data_genes))
#
## Number of genes per enhancer 
#data_genes_per_enhancer = subset.groupby(['name']).size().values
#np.savez("ENCODE-rE2G_extended/{}_genes_per_enhancer.npz".format(biosample), data_genes_per_enhancer)
#
## Number of enhancers per gene 
#data_enhancer_per_gene = subset.groupby(['TargetGene']).size().values
#np.savez("ENCODE-rE2G_extended/{}_enhancer_per_gene.npz".format(biosample), data_enhancer_per_gene)
#
## Size of enhancer regions 
#subset['size_regions'] = subset['end']-subset['start']
#size_regions = np.array(subset['size_regions'])
#np.savez("ENCODE-rE2G_extended/{}_enhancer_size.npz".format(biosample), size_regions)

