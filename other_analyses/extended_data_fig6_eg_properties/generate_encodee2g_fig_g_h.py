import pandas as pd
import numpy as np
import os, sys

# add biosample
biosample = sys.argv[1]

# /oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_8988T_ENCSR000EID_thresholded_predictions.tsv.gz
# for single file
data = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_{}_thresholded_predictions.tsv.gz".format(biosample), sep="\t")

subset = data.loc[data['class']!="promoter"]
subset = subset[['#chr', 'start', 'end']].drop_duplicates()
subset['Biosample'] = biosample

genes = data[['TargetGene']].drop_duplicates()
genes['Biosample'] = biosample
genes.to_csv("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_{}_active_genes_only.tsv.gz".format(biosample), sep="\t", index=False, header=False)

subset.to_csv("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_{}_active_enhancers_only.tsv.gz".format(biosample), sep="\t", index=False, header=False)
