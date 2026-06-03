## Rules to perform CRISPR benchmarking analyses presented in the paper


# Main figures -------------------------------------------------------------------------------------    

# perform benchmarking analyses of main predictors against combined K562 (trainint) CRISPR data and
# make plots for Figure 2
rule fig2bc_main_crispr_benchmarks:
  input:
    tss_annot = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["crispr_benchmark_dir"] + "/resources/genomic_features/K562_expressed_genes.tsv",
    perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    merged_data = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv"
  output: "results/manuscript/fig2bc_main_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G",
    runtime = "3h"
  script:
    "../scripts/crispr_benchmarking/fig2bc_main_crispr_benchmarks.Rmd"
    
# perform benchmarking analyses of main predictors against held-out CRISPR data and make plots for
# Figure 2d and Extended Data Figure 3e-g
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
    "../scripts/crispr_benchmarking/heldout_crispr_benchmarks.Rmd"
    
# make figures for assessing the performance of different chromatin assays in ABC shown in Figure 4d
# and Supplementary Figure 6
rule assaying_enhancer_activity:
  input:
    metadata_bigwig = "resources/assaying_enhancer_activity/processed_encode_chromatin_metadata_bigWig.tsv.gz",
    merged_data = config["crispr_benchmark_dir"] + "/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/EnhActABC_distal_reg_pred_config_bigWig.tsv",
  output: "results/manuscript/assaying_enhancer_activity.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"  
  script:
    "../scripts/crispr_benchmarking/assaying_enhancer_activity.Rmd"
    
# plot performance of ENCODE-rE2G models with additional assays for Figure 4f,g
rule main_fig4fg_additional_models:
  input:
    perf_summary = config["crispr_benchmark_dir"] + "/results/AdditionalModels/performance_summary.txt",
    merged_data = config["crispr_benchmark_dir"] + "/results/AdditionalModels/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    model_features = "resources/additional_model_features.tsv"
  output: "results/manuscript/fig4fg_additional_models.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/crispr_benchmarking/fig4fg_additional_models.Rmd"
    

# Extended data figures ----------------------------------------------------------------------------
    
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
    mem = "16G",
    runtime = "3h"
  script:
    "../scripts/crispr_benchmarking/extended_data_fig2_crispr_benchmarking.Rmd"

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
    "../scripts/crispr_benchmarking/extended_data_fig3a_perf_eg_classes.Rmd"

# analyze performance differences of additional ENCODE-rE2G models shown in Extended Data Figure 3
# panels b-d
rule extended_data_fig3bcd_additional_models:
  input:
    merged_data = config["crispr_benchmark_dir"] + "/results/AdditionalModelsAnalysis/expt_pred_merged_annot.txt.gz",
    pred_config = "CRISPR_comparison_resources/pred_config_files/additional_models_analysis_pred_config.tsv"
  output: "results/manuscript/extended_data_fig3bcd_additional_models.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/crispr_benchmarking/extended_data_fig3bcd_additional_models.Rmd"
    
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
    "../scripts/crispr_benchmarking/extended_data_fig7_enhancer_definitions.R"    


# Supplementary figures ---------------------------------------------------------------------------- 

# create CRISPR training data properties for Supplementary Figure N1.1 Note (Supplementary Note 1)
rule supplementary_note_N1_crispr_datasets:
  input:
    gasp_crispr = "CRISPR_comparison_resources/crispr_data/EPCrisprBenchmark_Gasperini2019_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
    schrai_crispr = "CRISPR_comparison_resources/crispr_data/EPCrisprBenchmark_Schraivogel2020_0.13gStd_0.8pwrAt15effect_GRCh38.tsv.gz",
    combined_crispr = "CRISPR_comparison/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    gene_universe = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    merged_training = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/supplementary_note_N1_crispr_datasets.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/crispr_benchmarking/supplementary_note_N1_crispr_datasets.Rmd"

# perform indirect effects analyses including probability of direct effects as function of distance
# to TSS shown in Supplementary Figure N1.2b (Supplementary Note 1)
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
    "../scripts/crispr_benchmarking/indirect_crispr_effects_analyses.Rmd"

# create locus plots for CRISPR example loci shown in Supplementary Figure 1
rule supplementary_fig1_crispr_locus_plots:
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
    hbe1 = "results/manuscript/plots/supplementary_fig1_crispr_locus_plots/hbe1_locus_plot.pdf",
    myc = "results/manuscript/plots/supplementary_fig1_crispr_locus_plots/myc_locus_plot.pdf",
    myc_zoomed = "results/manuscript/plots/supplementary_fig1_crispr_locus_plots/myc_locus_plot_zoomed.pdf",
    gata1 = "results/manuscript/plots/supplementary_fig1_crispr_locus_plots/gata1_locus_plot.pdf"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "128G"
  script:
    "../scripts/crispr_benchmarking/supplementary_fig1_crispr_locus_plots.R"

# perform additional CRISPR benchmarks shown in Supplementary Figure 2
rule supplementary_fig2_crispr_benchmarks:
  input:
    perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    tss_annot = config["crispr_benchmark_dir"] + "/resources/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed",
    expr_genes = config["crispr_benchmark_dir"] + "/resources/genomic_features/K562_expressed_genes.tsv",
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_nasser = config["crispr_benchmark_dir"] + "/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["crispr_benchmark_dir"] + "/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["crispr_benchmark_dir"] + "/results/Schraivogel2020/expt_pred_merged_annot.txt.gz",
    merged_main_preds = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/supplementary_fig2_crispr_benchmarks.html"
  params:
    seed = config["seed"],
    bootstrap_iterations = 10000
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G",
    runtime = "6h"
  script:
    "../scripts/crispr_benchmarking/supplementary_fig2_crispr_benchmarks.Rmd"

# plot correlation betwee predictor scores and CRISPR effect sizes shown in Supplementary Figure 3
rule supplementary_fig3_predictors_vs_crispr_effects:
  input:
    pred_config = "CRISPR_comparison_resources/pred_config_files/main_benchmarking_pred_config.tsv",
    merged_combined = config["crispr_benchmark_dir"] + "/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    merged_nasser = config["crispr_benchmark_dir"] + "/results/Nasser2021/expt_pred_merged_annot.txt.gz",
    merged_gasperini = config["crispr_benchmark_dir"] + "/results/Gasperini2019/expt_pred_merged_annot.txt.gz",
    merged_schraivogel = config["crispr_benchmark_dir"] + "/results/Schraivogel2020/expt_pred_merged_annot.txt.gz"
  output: "results/manuscript/supplementary_fig3_predictors_vs_crispr_effects.html"
  params:
    seed = config["seed"]
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
    "../scripts/crispr_benchmarking/supplementary_fig3_predictors_vs_crispr_effects.Rmd"


# Supplementary tables -----------------------------------------------------------------------------    

# create supplementary table with CRISPR benchmarking performance
rule supplementary_table3_crispr_performance_summary:
  input:
    main_perf_summary = config["crispr_benchmark_dir"] + "/results/MainPredictors/performance_summary.txt",
    published_perf_summary = config["crispr_benchmark_dir"] + "/results/PublishedPredictors/performance_summary.txt",
    baseline_dhs_perf_summary = config["crispr_benchmark_dir"] + "/results/BaselinePredsDHS/performance_summary.txt",
    baseline_abc_perf_summary = config["crispr_benchmark_dir"] + "/results/BaselinePredsABC/performance_summary.txt",
    merged_enhAct = config["crispr_benchmark_dir"] + "/results/EnhancerActivityBigWig/expt_pred_merged_annot.txt.gz",
    pred_config_enhAct = "CRISPR_comparison_resources/pred_config_files/EnhActABC_distal_reg_pred_config_bigWig.tsv"
  output: "results/manuscript/tables/supplementary_table3_crispr_performance_summary.csv"
  params:
    seed = config["seed"],
    compute_enhAct_performance = True
  conda: "../envs/analyses_env.yml"
  threads: 8
  resources:
    mem = "128G",
    runtime = "6h"
  script:
     "../scripts/crispr_benchmarking/supplementary_table3_crispr_performance_summary.R"
     
# create supplementary table with information on enhancer activity ABC models input assays
rule supplementary_table17_enhancer_activity_files:
  input: "resources/assaying_enhancer_activity/processed_encode_chromatin_metadata_bigWig.tsv.gz"
  output: "results/manuscript/tables/supplementary_table17_enhancer_activity_files.csv"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "16G"
  script:
     "../scripts/crispr_benchmarking/supplementary_table17_enhancer_activity_files.R"
