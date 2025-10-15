import pandas as pd
import numpy as np
import os, sys

# add biosample
biosample = sys.argv[1]

# for single file
data = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/EpiMap/processed/thresholded/linking_collated_{}.bed.gz.thresholded.gz".format(biosample), sep="\t", names=["chr","start","end","name","class","TargetGene","TargetGeneEnsemblID","TargetGeneTSS","CellType","Score","unnormalized.Score", "DistanceToTSS"])
tss = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/reference_files/RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS.bed", sep="\t", names=["chr_gene", "start_gene", "end_gene", "name_gene", "score_gene", "strand_gene", "tss_gene"])

merge = data.merge(tss, left_on=['TargetGene'], right_on=['name_gene'])

subset = merge.loc[merge['class']!="promoter"]
#subset = subset.loc[subset['Score']>=0.0045]

# get midpoint 
subset['midpoint'] = 0.5*(subset['start']+subset['end'])
subset['midpoint'] = subset['midpoint'].astype('int')

# log10(Distance to TSS)
subset['distance'] = subset['tss_gene'] - subset['midpoint']
subset['log10(Distance_to_TSS)'] = np.log10(subset['distance'])
distance_to_tss = np.array(subset['log10(Distance_to_TSS)'])
np.savez("EpiMap/{}_distance_to_tss.npz".format(biosample), distance_to_tss)

# Number of enhancer-gene links
subset['Num_enhancer_gene_links'] = len(subset)

with open('EpiMap/{}_num_enhancer_gene_links.txt'.format(biosample), 'w') as output:
    output.write(str(subset['Num_enhancer_gene_links'].values[0]))

# Number of genes with at least 1 distal enhancer
data_genes = len(subset.groupby(['TargetGene']))
with open('EpiMap/{}_num_genes_one_distal_enhancer.txt'.format(biosample), 'w') as output:
    output.write(str(data_genes))

# Number of genes per enhancer 
data_genes_per_enhancer = subset.groupby(['name']).size().values
np.savez("EpiMap/{}_genes_per_enhancer.npz".format(biosample), data_genes_per_enhancer)

# Number of enhancers per gene 
data_enhancer_per_gene = subset.groupby(['TargetGene']).size().values
np.savez("EpiMap/{}_enhancer_per_gene.npz".format(biosample), data_enhancer_per_gene)

# Size of enhancer regions 
subset['size_regions'] = subset['end']-subset['start']
size_regions = np.array(subset['size_regions'])
np.savez("EpiMap/{}_enhancer_size.npz".format(biosample), size_regions)

