
rule indirect_effects:
  input:
    perf_training = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/performance_summary.txt",
    cis_ratio = "heldout_crispr_indirect_effects/results/cis_positive_hit_ratio.tsv",
    trans_ratio = "heldout_crispr_indirect_effects/results/trans_positive_hit_ratio.tsv",
    crispr_training = config["share_dir"] + "/CRISPR_data/element_classes/EPCrisprBenchmark_combined_training_element_classes.tsv.gz",
    crispr_heldout = config["share_dir"] + "/CRISPR_data/element_classes/EPCrisprBenchmark_combined_heldout_element_classes.tsv.gz",
    merged_training = config["proj_dir"] + "/CRISPR_benchmarks/results/MainPredictors/expt_pred_merged_annot.txt.gz",
    merged_heldout = config["proj_dir"] + "/CRISPR_benchmarks/results/CombinedHeldout/expt_pred_merged_annot.txt.gz",
    pred_config = config["proj_dir"] + "/CRISPR_benchmarks/config/predictor_config_files/benchmarking_pred_config_with_alpha.tsv"
  output: "results/manuscript/indirect_effects.html"
  conda: "../envs/analyses_env.yml"
  resources:
    mem = "32G"  
  script:
    "../scripts/WIP/indirect_effects.Rmd"
