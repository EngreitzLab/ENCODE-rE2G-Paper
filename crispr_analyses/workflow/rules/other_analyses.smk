## Rules to perform some of the other analyses presented in the manuscript


# Main figures -------------------------------------------------------------------------------------    

# get ENCODE-rE2G scores for one gene needed to create locus plots showing scores in many cell types
rule get_e2g_scores_gene:
  input: "resources/encode_re2g_metadata.tsv"
  output: "results/encode_re2g_scores_{gene}_thresholded.rds"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/other_analyses/get_e2g_scores_gene.R"

# create locus plot showing ENCODE-rE2G predictions and CRISPR experimental validated E-G
# interactions for Figure 1b
rule fig1b_prkar2b_locus:
  input:
    k562_e2g_links = config["eg_predictions"]["encode_re2g_thresholded"],
    k562_e2g_ext_links = config["eg_predictions"]["encode_re2g_extended_thresholded"],
    crispr = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz",
    k562_dhs_peaks = "resources/ENCFF185XRG.bed.gz",
    k562_dnase_bigwig = "resources/ENCFF414OGC.bigWig",
    k562_h3k27ac_bigwig = "resources/ENCFF849TDM.bigWig",
    e2g_metadata = "resources/encode_re2g_metadata.tsv",
    e2g_scores_gene = "results/encode_re2g_scores_prkar2b_thresholded.rds"
  output: "results/manuscript/fig1b_prkar2b_locus.html"
  params:
    seed = config["seed"],
    filter_links_locus = True,
    min_enhancers = 1
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"  
  script:
    "../scripts/other_analyses/fig1b_prkar2b_locus.Rmd" 
    
# make locus plot for Figure 3e 
rule fig3e_gwas_locus:
  input:
    e2g_metadata = "resources/encode_re2g_metadata.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/manuscript/fig3e_gwas_locus.html"
  params:
    seed = config["seed"],
    scratch_dir = config["scratch_dir"],
    plot_outside = True
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G" 
  script:
    "../scripts/other_analyses/fig3e_gwas_locus.Rmd"

# create feature locus plot for Figure 4a
rule fig4a_feature_locus_plot:
  input:
    e2g_features = config["eg_predictions"]["encode_re2g"],
    e2g_links = config["eg_predictions"]["encode_re2g_bedpe"],
    feat_config = "resources/e2g_features_config.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/manuscript/plots/fig4a_prkar2b_features.pdf"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/other_analyses/fig4a_feature_locus_plot.R"


# Extended data figures ---------------------------------------------------------------------------- 

# plot correlation between ENCODE-rE2G_Extended features shown in Extended Data Figure 1
rule extended_data_fig1_feature_correlation:
  input:
    features = config["eg_predictions"]["encode_re2g_extended_locov"],
    feature_names = "resources/eg_feature_names.txt"
  output: "results/manuscript/extended_data_fig1_feature_correlation.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/other_analyses/extended_data_fig1_feature_correlation.Rmd"
    
# calculate ENCODE-rE2G properties for Extended Data Figure 6 panels a-c
rule compute_encode_re2g_properties:
  input: config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/encode_re2g_properties.rds"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "48G",
    runtime = "1h"
  script:
    "../scripts/other_analyses/compute_encode_re2g_properties.R"

# plot distributions of ENCODE-rE2G properties for Extended Data Figure 6 panels a-c
rule extended_data_fig6abc_property_plots:
  input: "results/encode_re2g_properties.rds"
  output: "results/manuscript/plots/extended_data_fig6abc_property_plots.pdf"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "48G",
    runtime = "1h"
  script:
    "../scripts/other_analyses/extended_data_fig6abc_property_plots.R"
    

# Supplementary figures ---------------------------------------------------------------------------- 

# make additional GWAS locus plots for Supplementary Figure N5.2
rule supplementary_figN5_additional_gwas_loci:
  input:
    e2g_metadata = "resources/encode_re2g_metadata.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output:
    rs4927708 = "results/manuscript/plots/supplementary_figN5.2_gwas_locus_plots/rs4927708_locus_plot.pdf",
    rs218265  = "results/manuscript/plots/supplementary_figN5.2_gwas_locus_plots/rs218265_locus_plot.pdf",
    rs7599488 = "results/manuscript/plots/supplementary_figN5.2_gwas_locus_plots/rs7599488_locus_plot.pdf"
  params:
    scratch_dir = config["scratch_dir"],
    plot_outside = True
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/other_analyses/supplementary_figN5_additional_gwas_loci.R"    

# create supplementary figure plotting ENCODE-rE2G scores against gene expression metrics
rule supplementary_fig4_encode_re2g_vs_gene_expression:
  input:
    pred = config["eg_predictions"]["encode_re2g_thresholded"],
    expr = "resources/gasperini_cpms.csv.gz"
  output: "results/manuscript/plots/supplementary_fig4_encode_re2g_vs_gene_expression.pdf"
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/other_analyses/supplementary_fig4_encode_re2g_vs_gene_expression.R"

# plot ENCODE-rE2G enhancer activity as function of distance to TSS for Supplementary Figure 5a
rule supplementary_fig5a_activity_vs_distance:
  input:
    pred_files = config["encode_re2g_predictions"]["thresholded"].values(),
    enh_files = config["encode_re2g_enhancer_lists"].values()
  output: "results/manuscript/plots/supplementary_fig5a_activity_vs_distance.pdf"
  threads: 16
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "128G",
    runtime = "6h"
  script:
    "../scripts/other_analyses/supplementary_fig5a_activity_vs_distance.R"
    

# Supplementary tables -----------------------------------------------------------------------------    
    
# create supplementary table with ENCODE-rE2G metadata
rule supplementary_table2_encode_re2g_metadata:
  input:
    e2g_meta = "resources/encode_re2g_metadata.tsv",
    e2g_portal_meta = "resources/encode_re2g_portal_accessions_20250812.tsv",
    extended_assays = "resources/encode_re2g_extended_assays.tsv"
  output:
    dnase_table = "results/manuscript/tables/supplementary_table2_encode_re2g_metadata.csv",
    extended_table = "results/manuscript/tables/supplementary_table2_encode_re2g_extended_metadata.csv",
    combined_table = "results/manuscript/tables/supplementary_table2_encode_re2g_metadata_combined.csv"
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/other_analyses/supplementary_table2_encode_re2g_metadata.R"

# create supplementary table with ENCODE-rE2G summary stats per biosample
rule supplementary_table8_encode_re2g_stats_per_biosample:
  input: config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/supplementary_table8_encode_re2g_stats_per_biosample.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/other_analyses/supplementary_table8_encode_re2g_stats_per_biosample.R"

# create supplementary table with ENCODE-rE2G summary stats per gene
rule supplementary_table9_encode_re2g_stats_per_gene:
  input:
    genes = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    pred = config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/supplementary_table9_encode_re2g_stats_per_gene.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/other_analyses/supplementary_table9_encode_re2g_stats_per_gene.R"
    
# create supplementary table with baseline predictor information
rule supplementary_table16_baseline_predictors:
  input:
    input_files = "resources/baseline_predictors/e2g_baseline_preds_input_files_bam.tsv",
    abc_elements_meta = "resources/metadata_encode_portal_abc_elements.tsv",
    full_pred_meta = "resources/metadata_encode_portal_full_predictions.tsv",
    dhs_synapse_manifest = "resources/baseline_predictors/synapse_manifest_dhs.tsv",
    abc_synapse_manifest = "resources/baseline_predictors/synapse_manifest_abc.tsv"
  output: "results/manuscript/tables/supplementary_table16_baseline_predictors.csv"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
     "../scripts/other_analyses/supplementary_table16_baseline_predictors.R"    
    
