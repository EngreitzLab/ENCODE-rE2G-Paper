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
  input: "CRISPR_comparison/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
  output: temp("results/enhancer_definitions/crispr/combined_crispr_k562_enhancers.bed.gz")
  conda: "../envs/analyses_env.yml"
  script:
    "../scripts/enhancer_definitions/extract_crispr_enhancers.R"

# extract distal enhancers from E-G predictions
rule extract_enhancers_predictions:
  input: config["eg_predictions"]["encode_re2g_thresholded"]
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
  resources:
    mem = "8G"
  shell:
    "bedtools intersect -a {input.enh} -b {input.univ} -loj | gzip > {output}"
