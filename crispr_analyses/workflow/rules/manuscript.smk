## Rules to create manuscript figures and tables

# get ENCODE-rE2G scores for one gene needed to create locus plots showing scores in many cell types
rule get_e2g_scores_gene:
  input: "resources/encode_re2g_metadata.tsv"
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
    k562_e2g_links = config["eg_predictions"]["encode_re2g_thresholded"],
    k562_e2g_ext_links = config["eg_predictions"]["encode_re2g_extended_thresholded"],
    crispr = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz",
    k562_dhs_peaks = "resources/ENCFF185XRG.bed.gz",
    k562_dnase_bigwig = "resources/ENCFF414OGC.bigWig",
    k562_h3k27ac_bigwig = "resources/ENCFF849TDM.bigWig",
    e2g_metadata = "resources/encode_re2g_metadata.tsv",
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
    tss_annot = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["crispr_benchmark_dir"] + "/resources/genomic_features/K562_expressed_genes.tsv",
    perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    merged_data = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv"
  output: "results/manuscript/main_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/main_crispr_benchmarks.Rmd"
    
# perform benchmarking analyses of main predictors against held-out CRISPR data and make plots for
# Figure 2 and Extended Data Figure 3
rule heldout_crispr_benchmarks:
  input:
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_training = config["crispr_benchmark_dir"] + "/results/CominedTrainingFiltered/expt_pred_merged_annot.txt.gz",
    merged_heldout = config["crispr_benchmark_dir"] + "/results/CombinedHeldoutFiltered/expt_pred_merged_annot.txt.gz"
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
    
# make locus plot for Figure 4e 
rule main_fig4e_gwas_locus:
  input:
    e2g_metadata = "resources/encode_re2g_metadata.tsv",
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
    
# create feature locus plot for Figure 5a
rule main_fig5a_feature_locus_plot:
  input:
    e2g_features = config["eg_predictions"]["encode_re2g"],
    e2g_links = config["eg_predictions"]["encode_re2g_bedpe"],
    feat_config = "resources/e2g_features_config.tsv",
    genome_annot = "resources/gencode.v29.annotation.gtf.gz"
  output: "results/manuscript/plots/main_fig5a_prkar2b_features.pdf"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/main_fig5a_feature_locus_plot.R"
    
# make figures for assessing the performance of different chromatin assays in ABC shown in Figure 5d
# and Figure S10
rule assaying_enhancer_activity_figures:
  input:
    metadata_bigwig = "resources/assaying_enhancer_activity/processed_encode_chromatin_metadata_bigWig.tsv.gz",
    merged_data = config["crispr_benchmark_dir"] + "/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/EnhActABC_distal_reg_pred_config_bigWig.tsv",
  output: "results/manuscript/assaying_enhancer_activity_figures.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"  
  script:
    "../scripts/manuscript/assaying_enhancer_activity_figures.Rmd"    
    
# plot performance of ENCODE-rE2G models with additional assays for Figure 5
rule main_fig5_additional_models:
  input:
    perf_summary = config["crispr_benchmark_dir"] + "/results/AdditionalModels/performance_summary.txt",
    merged_data = config["crispr_benchmark_dir"] + "/results/AdditionalModels/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
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
    "../scripts/manuscript/extended_data_fig1_feature_correlation.Rmd"

# CRISPR benchmarking analyses to create Extended Data Figure 2    
rule extended_data_fig2_crispr_benchmarking:
  input:
    tss_annot = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["crispr_benchmark_dir"] + "/resources/genomic_features/K562_expressed_genes.tsv",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_baseline_dhs = config["crispr_benchmark_dir"] + "/results/BaselinePredsDHS/expt_pred_merged_annot.txt.gz",
    merged_baseline_abc = config["crispr_benchmark_dir"] + "/results/BaselinePredsABC/expt_pred_merged_annot.txt.gz",
    merged_published_preds = config["crispr_benchmark_dir"] + "/results/PublishedPredictors/expt_pred_merged_annot.txt.gz",
    merged_main_preds = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
  output: "results/manuscript/extended_data_fig2_crispr_benchmarking.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/extended_data_fig2_crispr_benchmarking.Rmd"

# additional benchmarks to show performance on different CRISPR E-G pair classes shown in Extended
# Data Figure 3
rule extended_data_fig3a_perf_eg_classes:
  input:
    merged_data = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv"
  output: "results/manuscript/extended_data_fig3a_perf_eg_classes.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G",
    runtime = "3h"
  script:
    "../scripts/manuscript/extended_data_fig3a_perf_eg_classes.Rmd"

# analyses performance differences of additional ENCODE-rE2G models shown in Extended Data Figure 3   
rule extended_data_fig3_additional_models:
  input:
    merged_data = config["crispr_benchmark_dir"] + "/results/AdditionalModelsAnalysis/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/additional_models_analysis_pred_config.tsv"
  output: "results/manuscript/extended_data_fig3_additional_models.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/extended_data_fig3_additional_models.Rmd"

# plot ENCODE-rE2G enhancer activity as function of distance to TSS for extended figure 6    
rule extended_data_fig6i_activity_vs_distance:
  input: 
    pred_files = config["encode_re2g_predictions"]["thresholded"].values(),
    enh_files = config["encode_re2g_enhancer_lists"].values()
  output: "results/manuscript/plots/extended_data_fig6i_activity_vs_distance.pdf"
  threads: 16
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "128G",
    runtime = "6h"
  script:
    "../scripts/manuscript/extended_data_fig6i_activity_vs_distance.R"

# create extended data figure 7 showing comparison to other enhancer definitions
rule extended_data_fig7_enhancer_definitions:
  input:
    e2g_ccres = "results/encode_re2g/encode_re2g_k562_enhancers.ENCFF286VQG.bed.gz",
    e2g_chromhmm = "results/encode_re2g/encode_re2g_k562_enhancers.ENCFF963KIA.bed.gz",
    crispr_ccres = "results/crispr/combined_crispr_k562_enhancers.ENCFF286VQG.bed.gz",
    crispr_chromhmm = "results/crispr/combined_crispr_k562_enhancers.ENCFF963KIA.bed.gz"
  output: "results/manuscript/plots/extended_data_fig7_enhancer_definitions.pdf"
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/extended_data_fig7_enhancer_definitions.R"
    
# create CRISPR training data properties figure for Note S1
rule crispr_supplementary_note_S1_1:
  input:
    gasp_crispr = "CRISPR_comparison_resources/crispr_data/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
    schrai_crispr = "CRISPR_comparison_resources/crispr_data/EPCrisprBenchmark_Schraivogel2020_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
    combined_crispr = "CRISPR_comparison/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    gene_universe = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    merged_training = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/crispr_supplementary_note_S1.1.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/crispr_supplementary_note_S1.1.Rmd"

# perform indirect effects analyses including probability of direct effects as function of distance
# to TSS shown in Figure S1.2b
rule indirect_crispr_effects_analyses:
  input:
    cis_rates = "resources/indirect_effects/cis_positive_hit_rates.tsv",
    trans_rates = "resources/indirect_effects/trans_positive_hit_rates_imputed.tsv",
    direct_rates_datasets = "resources/indirect_effects/direct_effect_models/direct_rates_per_dataset.tsv",
    direct_rates_average  = "resources/indirect_effects/direct_effect_models/direct_rates_average_across_datasets.tsv",
    dataset_ids = "resources/indirect_effects/crispr_dataset_ids.tsv",
    direct_effect_models = "resources/indirect_effects/direct_effect_models/direct_effects_model.rds",
    merged_training = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    merged_heldout = config["crispr_benchmark_dir"] + "/results/CombinedHeldout/expt_pred_merged_annot.txt.gz",
  output: "results/manuscript/indirect_effects_analyses.html"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/indirect_crispr_effects_analyses.Rmd"
    
# perform additional CRISPR benchmarks shown in Fig S8    
rule figS7_additional_crispr_benchmarks:
  input:
    perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    tss_annot = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["crispr_benchmark_dir"] + "/resources/genomic_features/K562_expressed_genes.tsv",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_nasser = config["crispr_benchmark_dir"] + "/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["crispr_benchmark_dir"] + "/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["crispr_benchmark_dir"] + "/results/Schraivogel2020/expt_pred_merged_annot.txt.gz",
    merged_main_preds = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/figS7_additional_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/figS7_additional_crispr_benchmarks.Rmd"

# plot correlation betwee predictor scores and CRISPR effect sizes     
rule figS8_predictors_vs_crispr_effects:
  input:
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_combined = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    merged_nasser = config["crispr_benchmark_dir"] + "/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["crispr_benchmark_dir"] + "/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["crispr_benchmark_dir"] + "/results/Schraivogel2020/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/figS8_predictors_vs_crispr_effects.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/manuscript/figS8_predictors_vs_crispr_effects.Rmd"

# create locus plots for CRISPR example loci shown in supplementary figures     
rule figS9_crispr_locus_plots:
  input:
    genome_annot = "resources/gencode.v29.annotation.gtf.gz",
    crispr_data = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    dhs_peaks = "resources/ENCFF185XRG.bed.gz",
    dnase_bam = config["scratch_dir"] + "/ENCFF860XAE.sorted.bam",
    dnase_bai = config["scratch_dir"] + "/ENCFF860XAE.sorted.bam.bai",
    h3k27ac_bam = config["scratch_dir"] + "/ENCFF790GFL.sorted.bam",
    h3k27ac_bai = config["scratch_dir"] + "/ENCFF790GFL.sorted.bam.bai",
    e2g_links = config["eg_predictions"]["encode_re2g_thresholded"],
    e2g_ext_links = config["eg_predictions"]["encode_re2g_extended_thresholded"]
  output:
    hbe1 = "results/manuscript/plots/figS9_crispr_locus_plots/hbe1_locus_plot.pdf",
    myc = "results/manuscript/plots/figS9_crispr_locus_plots/myc_locus_plot.pdf",
    myc_zoomed = "results/manuscript/plots/figS9_crispr_locus_plots/myc_locus_plot_zoomed.pdf",
    gata1 = "results/manuscript/plots/figS9_crispr_locus_plots/gata1_locus_plot.pdf"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "128G"
  script:
    "../scripts/manuscript/figS9_crispr_locus_plots.R"

# create supplementary figure plotting ENCODE-rE2G scores against gene expression metrics
rule figS11_encode_re2g_vs_gene_expression:
  input:
    pred = config["eg_predictions"]["encode_re2g_thresholded"],
    expr = "resources/gasperini_cpms.csv.gz"
  output: "results/manuscript/plots/figS11_encode_re2g_vs_gene_expression.pdf"
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/figS11_encode_re2g_vs_gene_expression.R"

# create supplementary table with ENCODE-rE2G metadata
rule table_s2_encode_re2g_metadata:
  input:
    e2g_meta = "resources/encode_re2g_metadata.tsv",
    e2g_portal_meta = "resources/encode_re2g_portal_accessions_20250812.tsv",
    extended_assays = "resources/encode_re2g_extended_assays.tsv"
  output:
    dnase_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_predictions.csv",
    extended_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_ext_predictions.csv",
    combined_table = "results/manuscript/tables/table_s2_metadata_encode_re2g_predictions_combined.csv"
  conda: "../envs/analyses_env.yml" 
  script:
    "../scripts/manuscript/table_s2_encode_re2g_metadata.R"

# create supplementary table with baseline predictor information
rule table_s9_baseline_predictors:
  input:
    input_files = "resources/baseline_predictors/e2g_baseline_preds_input_files_bam.tsv",
    abc_elements_meta = "resources/metadata_encode_portal_abc_elements.tsv",
    full_pred_meta = "resources/metadata_encode_portal_full_predictions.tsv",
    dhs_synapse_manifest = "resources/baseline_predictors/synapse_manifest_dhs.tsv",
    abc_synapse_manifest = "resources/baseline_predictors/synapse_manifest_abc.tsv"
  output: "results/manuscript/tables/table_s9_baseline_predictors.csv"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "16G"
  script:
     "../scripts/manuscript/table_s9_baseline_predictors.R"

# create supplementary table with CRISPR benchmarking performance
rule table_s10_performance_summary:
  input:
    main_perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    published_perf_summary = config["crispr_benchmark_dir"] + "/results/PublishedPredictors/performance_summary.txt",
    baseline_dhs_perf_summary = config["crispr_benchmark_dir"] + "/results/BaselinePredsDHS/performance_summary.txt",
    baseline_abc_perf_summary = config["crispr_benchmark_dir"] + "/results/BaselinePredsABC/performance_summary.txt",
    merged_enhAct = config["crispr_benchmark_dir"] + "/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config_enhAct = "CRISPR_comparison_resources/pred_config_files/EnhActABC_distal_reg_pred_config_bigWig.tsv"
  output: "results/manuscript/tables/table_s10_performance_summary.csv"
  params:
    seed = config["seed"],
    compute_enhAct_performance = True
  conda: "../envs/analyses_env.yml"  
  threads: 8
  resources:
    mem = "128G",
    runtime = "6h"
  script:
     "../scripts/manuscript/table_s10_performance_summary.R"

# create supplementary table with information on enhancer activity ABC models input assays  
rule table_s15_enhancer_activity_files:
  input: "resources/assaying_enhancer_activity/processed_encode_chromatin_metadata_bigWig.tsv.gz"
  output: "results/manuscript/tables/table_s15_enhancer_activity_files.csv"
  conda: "../envs/analyses_env.yml"  
  resources:
    mem = "16G"
  script:
     "../scripts/manuscript/table_s15_enhancer_activity_files.R"

# create supplementary table with ENCODE-rE2G summary stats per biosample
rule table_s16_encode_re2g_stats_per_experiment:
  input: config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/table_s16_encode_re2g_stats_per_experiment.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/table_s16_encode_re2g_stats_per_experiment.R"

# create supplementary table with ENCODE-rE2G summary stats per gene
rule table_s17_encode_re2g_stats_per_gene:
  input: 
    genes = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    pred = config["encode_re2g_predictions"]["thresholded"].values()
  output: "results/manuscript/tables/table_s17_encode_re2g_stats_per_gene.csv"
  threads: 10
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/manuscript/table_s17_encode_re2g_stats_per_gene.R"
