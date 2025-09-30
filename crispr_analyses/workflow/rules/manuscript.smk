## Rules to create manuscript figures and tables

# get ENCODE-rE2G scores for one gene needed to create locus plots showing scores in many cell types
rule get_e2g_scores_gene:
  input: config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv"
  output: "results/manuscript/encode_re2g_scores_{gene}_thresholded.rds"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"
  script:
    "../scripts/manuscript/get_e2g_scores_gene.R"

# create locus plot showing ENCODE-rE2G predictions and CRISPR experimental validated E-G
# interactions for Figure 1b
rule main_fig1b_prkar2b_locus:
  input:
    k562_e2g_links = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_K562_ENCSR000EOT_thresholded_predictions.tsv.gz",
    k562_e2g_ext_links = config["share_dir"] + "/Predictors/ENCODE-rE2G/additional_models/formatted/Extended/K562/encode_e2g_predictions_threshold0.336.tsv.gz",
    crispr = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz",
    k562_dhs_peaks = "resources/ENCFF185XRG.bed.gz",
    k562_dnase_bigwig = "resources/ENCFF414OGC.bigWig",
    k562_h3k27ac_bigwig = "resources/ENCFF849TDM.bigWig",
    e2g_metadata = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv",
    e2g_scores_gene = "results/manuscript/encode_re2g_scores_prkar2b_thresholded.rds"
  output: "results/manuscript/main_fig1b_prkar2b_locus.html"
  params:
    seed = config["seed"],
    filter_links_locus = True,
    min_enhancers = 1
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"  
  script:
    "../scripts/manuscript/main_fig1b_prkar2b_locus.Rmd"    

# perform benchmarking analyses of main predictors against training CRISPR data and make plots for
# Figure 2
rule main_crispr_benchmarks:
  input:
    tss_annot = config["proj_dir"] + "/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["proj_dir"] + "/CRISPR_benchmarks/resources/genomic_features/K562_expressed_genes.tsv",
    perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/performance_summary.txt",
    merged_data = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv"
  output: "results/manuscript/main_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/main_crispr_benchmarks.Rmd"

# additional benchmarks to show performance on different E-G pair classes
rule additional_crispr_benchmarks:
  input:
    merged_data = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv"
  output: "results/manuscript/additional_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G",
    runtime = "3h"
  script:
    "../scripts/manuscript/additional_crispr_benchmarks.Rmd"

# perform indirect effects analyses ## TODO: CLEAN UP 
# rule indirect_crispr_effects:
#   input:
#     cis_rates = config["proj_dir"] + "/CRISPR_indirect_effects/results/encode_analyses/cis_positive_hit_rates.tsv",
#     trans_rates = config["proj_dir"] + "/CRISPR_indirect_effects/results/encode_analyses/trans_positive_hit_rates_imputed.tsv",
#     direct_rates_datasets = config["proj_dir"] + "/CRISPR_indirect_effects/results/encode_analyses/direct_effect_models/direct_rates_per_dataset.tsv",
#     direct_rates_average  = config["proj_dir"] + "/CRISPR_indirect_effects/results/encode_analyses/direct_effect_models/direct_rates_average_across_datasets.tsv",    
#     dataset_ids = config["proj_dir"] + "/CRISPR_indirect_effects/config/crispr_dataset_ids.tsv",
#     direct_effect_models = config["proj_dir"] + "/CRISPR_indirect_effects/results/encode_analyses/direct_effect_models/direct_effects_model.rds",
#     merged_training = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
#     merged_heldout = config["proj_dir"] + "/CRISPR_benchmarks/results/CombinedHeldout/expt_pred_merged_annot.txt.gz",
#   output: "results/manuscript/indirect_effects.html"
#   conda: "../envs/analyses_env.yml"
#   resources:
#     mem = "16G"  
#   script:
#     "../scripts/manuscript/indirect_crispr_effects.Rmd"
    
# perform benchmarking analyses of main predictors against held-out CRISPR data and make plots for
# Figure 2 and supplementary material (need to run 'indirect_crispr_effects' first)
rule heldout_crispr_benchmarks:
  input:
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv",
    merged_training = config["proj_dir"] + "/CRISPR_benchmarks/results/CominedTrainingFiltered/expt_pred_merged_annot.txt.gz",
    merged_heldout = config["proj_dir"] + "/CRISPR_benchmarks/results/CombinedHeldoutFiltered/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/heldout_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G",
    runtime = "3h"
  script:
    "../scripts/manuscript/heldout_crispr_benchmarks.Rmd"

# make locus plot for Figure 4 
rule main_fig4e_gwas_locus:
  input:
    e2g_metadata = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/manuscript/main_fig4e_gwas_locus.html"
  params:
    seed = config["seed"],
    scratch_dir = config["scratch_dir"],
    filter_links_locus = True
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"  
  script:
    "../scripts/manuscript/main_fig4e_gwas_locus.Rmd"

# plot performance of ENCODE-rE2G models with additional assays for main figure 5
rule main_fig5_additional_models:
  input:
    perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/AdditionalModels/performance_summary.txt",
    merged_data = config["proj_dir"] + "/CRISPR_benchmarks/results/AdditionalModels/expt_pred_merged_annot.txt.gz",
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv",
    model_features = "resources/additional_model_features.tsv"
  output: "results/manuscript/main_fig5_additional_models.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/main_fig5_additional_models.Rmd"

# analyses performance differences of additional ENCODE-rE2G models for supplementary materials   
rule extended_data_fig3_additional_models:
  input:
    merged_data = config["proj_dir"] + "/CRISPR_benchmarks/results/AdditionalModelsAnalysis/expt_pred_merged_annot.txt.gz",
    pred_config = config["proj_dir"] + "/CRISPR_benchmarks/config/predictor_config_files/additional_models_pred_config.tsv"
  output: "results/manuscript/extended_data_fig3_additional_models.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/extended_data_fig3_additional_models.Rmd"
    
# make figures for assessing the performance of different chromatin assays in ABC    
rule assaying_enhancer_activity_figures:
  input:
    metadata_bigwig = config["proj_dir"] + "/predictors/enhancer_activity/resources/processed_encode_chromatin_metadata_bigWig.tsv.gz",
    merged_data = config["proj_dir"] + "/CRISPR_benchmarks/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config = config["proj_dir"] + "/predictors/enhancer_activity/results/bigWig/K562/EnhActABC_distal_reg_pred_config_bigWig.tsv",
  output: "results/manuscript/assaying_enhancer_activity_figures.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"  
  script:
    "../scripts/manuscript/assaying_enhancer_activity_figures.Rmd"

## TODO: FIX THIS
# rule crispr_supplementary_note_S1_1:
#   input:
#     gasp_crispr = config["share_dir"] + "/CRISPR_data/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
#     schrai_crispr = config["share_dir"] + "/CRISPR_data/EPCrisprBenchmark_TAPseq_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
#     combined_crispr = config["share_dir"] + "/CRISPR_data/indirect_effects/formatted/EPCrisprBenchmark.combined_training.annotated.tsv.gz",
#     gene_universe = config["proj_dir"] + "/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
#     merged_training = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz"
#   output: "results/manuscript/crispr_supplementary_note_S1.1.html"
#   params:
#     seed = config["seed"]
#   conda: "../envs/analyses_env.yml"
#   resources:
#     mem = "16G"
#   script:
#     "../scripts/manuscript/crispr_supplementary_note_S1_1.Rmd"

# plot correlation between ENCODE-rE2G features
rule extended_data_fig1_feature_correlation:
  input:
    features = config["share_dir"] + "/Predictors/ENCODE-rE2G/additional_models/Extended/K562/rE2G_LOCOV_crispr_predictions.tsv",
    feature_names = "resources/eg_feature_names.txt"
  output: "results/manuscript/extended_data_fig1_feature_correlation.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/extended_data_fig1_feature_correlation.Rmd"
    
rule extended_data_fig2_crispr_benchmarking:
  input:
    tss_annot = config["proj_dir"] + "/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["proj_dir"] + "/CRISPR_benchmarks/resources/genomic_features/K562_expressed_genes.tsv",
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv",
    merged_baseline_dhs = config["proj_dir"] + "/CRISPR_benchmarks/results/BaselinePredsDHS/expt_pred_merged_annot.txt.gz",
    merged_baseline_abc = config["proj_dir"] + "/CRISPR_benchmarks/results/BaselinePredsABC/expt_pred_merged_annot.txt.gz",
    merged_published_preds = config["proj_dir"] + "/CRISPR_benchmarks/results/PublishedPredictors/expt_pred_merged_annot.txt.gz",
    merged_main_preds = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
  output: "results/manuscript/extended_data_fig2_crispr_benchmarking.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/extended_data_fig2_crispr_benchmarking.Rmd"
    
rule figS8_additional_crispr_benchmarks:
  input:
    perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/performance_summary.txt",
    tss_annot = config["proj_dir"] + "/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["proj_dir"] + "/CRISPR_benchmarks/resources/genomic_features/K562_expressed_genes.tsv",
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv",
    merged_nasser = config["proj_dir"] + "/CRISPR_benchmarks/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["proj_dir"] + "/CRISPR_benchmarks/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["proj_dir"] + "/CRISPR_benchmarks/results/Schraivogel2020/expt_pred_merged_annot.txt.gz",
    merged_main_preds = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/figS8_additional_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/figS8_additional_crispr_benchmarks.Rmd"
    
rule figS9_predictors_vs_crispr_effects:
  input:
    pred_config = config["share_dir"] + "/Predictors/benchmarking_pred_config.tsv",
    merged_combined = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    merged_nasser = config["proj_dir"] + "/CRISPR_benchmarks/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["proj_dir"] + "/CRISPR_benchmarks/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["proj_dir"] + "/CRISPR_benchmarks/results/Schraivogel2020/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/figS9_predictors_vs_crispr_effects.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/figS9_predictors_vs_crispr_effects.Rmd"

rule main_fig5a_feature_locus_plot:
  input:
    e2g_features = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/full_predictions/encode_e2g_predictions_K562_ENCSR000EOT_full_predictions.tsv.gz",
    e2g_links = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/thresholded_bedpe/encode_e2g_predictions_K562_ENCSR000EOT_thresholded_predictions.bedpe.gz",
    feat_config = "resources/e2g_features_config.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/manuscript/plots/main_fig5a_prkar2b_features.pdf"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/main_fig5a_feature_locus_plot.R"
    
rule supplementary_fig_locus_plots:
  input:
    genome_annot = "resources/gencode.v29.annotation.gtf.gz",
    crispr_data = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    dhs_peaks = "resources/ENCFF185XRG.bed.gz",
    dnase_bam = config["scratch_dir"] + "/ENCFF860XAE.sorted.bam",
    dnase_bai = config["scratch_dir"] + "/ENCFF860XAE.sorted.bam.bai",
    h3k27ac_bam = config["scratch_dir"] + "/ENCFF790GFL.sorted.bam",
    h3k27ac_bai = config["scratch_dir"] + "/ENCFF790GFL.sorted.bam.bai",
    e2g_links = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_K562_ENCSR000EOT_thresholded_predictions.tsv.gz",
    e2g_ext_links = config["share_dir"] + "/Predictors/ENCODE-rE2G/additional_models/formatted/Extended/K562/encode_e2g_predictions_threshold0.336.tsv.gz"
  output:
    hbe1 = "results/manuscript/plots/supp_crispr_locus_plots/hbe1_locus_plot.pdf",
    myc = "results/manuscript/plots/supp_crispr_locus_plots/myc_locus_plot.pdf",
    myc_zoomed = "results/manuscript/plots/supp_crispr_locus_plots/myc_locus_plot_zoomed.pdf",
    gata1 = "results/manuscript/plots/supp_crispr_locus_plots/gata1_locus_plot.pdf"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "128G"
  script:
    "../scripts/manuscript/supplementary_fig_locus_plots.R"    

rule table_s2_encode_re2g_metadata:
  input:
    e2g_meta = config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/encode_re2g_metadata.tsv",
    e2g_portal_meta = "resources/encode_re2g_portal_accessions_20250812.tsv",
    extended_assays = "resources/encode_re2g_extended_assays.tsv"
  output:
    dnase_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_predictions.csv",
    extended_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_ext_predictions.csv",
    combined_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_predictions_combined.csv"
  conda: "../envs/analyses_env.yml" 
  script:
    "../scripts/manuscript/table_s2_encode_re2g_metadata.R"

rule table_s11_baseline_predictors:
  input:
    input_files = config["proj_dir"] + "/predictors/baseline_predictors/config/e2g_baseline_preds_input_files_bam.tsv",
    abc_elements_meta = "resources/metadata_encode_portal_abc_elements.tsv",
    full_pred_meta = "resources/metadata_encode_portal_full_predictions.tsv",
    dhs_synapse_manifest = config["share_dir"] + "/Predictors/BaselinePredictors/synapse_manifest_dhs.tsv",
    abc_synapse_manifest = config["share_dir"] + "/Predictors/BaselinePredictors/synapse_manifest_abc.tsv"
  output: "results/manuscript/tables/table_S11_baseline_predictors.csv"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "16G"
  script:
     "../scripts/manuscript/table_s11_baseline_predictors.R"
 
rule table_s12_performance_summary:
  input:
    main_perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/performance_summary.txt",
    published_perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/PublishedPredictors/performance_summary.txt",
    baseline_dhs_perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/BaselinePredsDHS/performance_summary.txt",
    baseline_abc_perf_summary = config["proj_dir"] + "/CRISPR_benchmarks/results/BaselinePredsABC/performance_summary.txt",
    merged_enhAct = config["proj_dir"] + "/CRISPR_benchmarks/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config_enhAct = config["proj_dir"] + "/predictors/enhancer_activity/results/bigWig/K562/EnhActABC_distal_reg_pred_config_bigWig.tsv"
  output: "results/manuscript/tables/table_s12_performance_summary.csv"
  params:
    seed = config["seed"],
    compute_enhAct_performance = False
  conda: "../envs/analyses_env.yml"  
  threads: 8
  resources:
    mem = "128G",
    runtime = "6h"
  script:
     "../scripts/manuscript/table_s12_performance_summary.R"
  
rule table_s17_enhancer_activity_files:
  input: config["proj_dir"] + "/predictors/enhancer_activity/resources/processed_encode_chromatin_metadata_bigWig.tsv.gz"
  output: "results/manuscript/tables/table_s17_enhancer_activity_files.csv"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "16G"
  script:
     "../scripts/manuscript/table_s17_enhancer_activity_files.R"

rule table_s18_encode_re2g_stats_per_experiment:
  input: config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/table_s18_encode_re2g_stats_per_experiment.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/table_s18_encode_re2g_stats_per_experiment.R"
  
rule table_s19_encode_re2g_stats_per_gene:
  input: 
    genes = config["proj_dir"] + "/CRISPR_benchmarks/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    pred = config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/table_s19_encode_re2g_stats_per_gene.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/table_s19_encode_re2g_stats_per_gene.R"
