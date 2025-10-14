import pandas as pd
import numpy as np
import os, sys
import pyranges as pr
from pathlib import Path

# 2) Load all biosample enhancer BEDs into one table with a biosample label
#    Assumes files like encode_enhancers_<biosample>.bed
biosample_files = list(Path("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/").glob("encode_e2g_predictions*active_enhancers_only.tsv.gz"))

rows = []
for f in biosample_files:
    biosample = f.stem  # or parse from filename however you like
    df = pd.read_csv(f, sep="\t", header=None, names=["Chromosome","Start","End", "Biosample"])
    rows.append(df)

all_bio = pd.concat(rows, ignore_index=True)
all_bio.to_csv("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_ALL_active_enhancers.tsv.gz", sep="\t", index=False)


# 1) Obtain list of active enhancers for each biosample
# & Obtain list of active genes for each biosample
# READ IN BIOSAMPLE NAME
biosample = sys.argv[1]

# for single file
data = pd.read_csv("/oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_{}_thresholded_predictions.tsv.gz".format(biosample), sep="\t")

subset = data.loc[data['class']!="promoter"]
subset = subset[['#chr', 'start', 'end']].drop_duplicates()
subset['Biosample'] = biosample

genes = data[['TargetGene']].drop_duplicates()
genes['Biosample'] = biosample
genes.to_csv("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_{}_active_genes_only.tsv.gz".format(biosample), sep="\t", index=False, header=False)

subset.to_csv("/oak/stanford/groups/engreitz/Users/kmualim/forRevisions/e2g_properties/IntermediateFiles/encode_e2g_predictions_{}_active_enhancers_only.tsv.gz".format(biosample), sep="\t", index=False, header=False)
