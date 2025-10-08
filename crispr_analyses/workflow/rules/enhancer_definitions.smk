## Rules to compute inputs for Extended Data Figure 7

# download file containing ENCODE cCREs
rule download_ccres:
  output: "resources/ENCFF286VQG.bed.gz"
  params:
     url = config["download_urls"]["cCREs"]
  conda: "../envs/analyses_env.yml"
  shell:
    "wget -O {output} {params.url}"
  
# download file containing ENCODE chromHMM chromatin states
rule download_chromhmm:
  output: "resources/ENCFF963KIA.bed.gz"
  params:
     url = config["download_urls"]["chromHMM"]
  conda: "../envs/analyses_env.yml"
  shell:
    "wget -O {output} {params.url}"
    
# extract enhancers from combined CRISPR data
rule extract_crispr_enhancers:
  input: config["share_dir"] + "/CRISPR_data/indirect_effects/formatted/EPCrisprBenchmark.combined_training.annotated.tsv.gz"
  output: temp("results/enhancer_definitions/crispr/combined_crispr_k562_enhancers.bed.gz")
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/enhancer_definitions/extract_crispr_enhancers.R"

# extract distal enhancers from E-G predictions
rule extract_enhancers_predictions:
  input: config["share_dir"] + "/Predictors/ENCODE-rE2G/dhs_only/thresholded_predictions/encode_e2g_predictions_K562_ENCSR000EOT_thresholded_predictions.tsv.gz"
  output: temp("results/enhancer_definitions/encode_re2g/encode_re2g_k562_enhancers.bed.gz")
  conda: "../envs/analyses_env.yml"
  shell:
    "zcat {input} | sed 1d | "
    """awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, $5}}' | """
    "sort | uniq | tr -d '\"' | gzip > {output}"

# overlap enhancers with element universe file
rule overlap_enhancers:
  input:
    enh = "results/enhancer_definitions/{predictor}/{sample}_enhancers.bed.gz",
    univ = "resources/{univ}.bed.gz"
  output: "results/{predictor}/{sample}_enhancers.{univ}.bed.gz"
  conda: "../envs/analyses_env.yml"
  shell:
    "bedtools intersect -a {input.enh} -b {input.univ} -loj | gzip > {output}"
