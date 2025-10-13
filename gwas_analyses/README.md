# GWAS analysis

This directory contains the code that was used to perform GWAS related presented in the "Linking
noncoding variants to target genes and cell types" section.

## Analysis scripts
Following scripts were used to intersect GWAS variants with enhancer-gene predicitons and to 
compute GWAS metrics used in analyses.
- **create_ALL_EG_predictions.R**: create enhancer-gene predictions pooled across all biosamples
(focusing on thresholded links)
- **create_BLD_EG_predictions.R**: create enhancer-gene predictions pooled across blood-related
biosamples (focusing on thresholded links)
- **top10_biosamples_ENCODE_E2G_enrichment.R**: For each GWAS trait, compute enriched representation
of enhancers linked to genes for each biosample, which can be used to select the top 10 biosamples
for each trait
- **PoPS_enrichment_top10_biosamples.R**: PoPS score enrichment for genes linked to finemapped
variants for that trait by ENCODE-rE2G in the top 10 globally enriched biosamples from the previous
code, compared to all biosamples. 
- **finemap_variants_unique_topbiosample.R**: For each trait, the number of fine-mapped variants
(PIP >  0.1) that are mapped to a gene by ENCODE-rEG in (a) any of the 1458 biosamples, (b) top 1
globally enriched biosample, and (c) top 10 globally enriched biosamples
- **ENCODE_E2G_PoPS_match_finemapvariant.R**: For each trait, the number of fine-mapped variants
(PIP	>  0.1) that are mapped to a gene by ENCODE-rE2G, PoPS and ENCODE-rE2G + PoPS
- **ENCODE_E2G_PoPS_match_variant.R**: Cataloguing the entire set of common+low-freq variants that
are mapped by ENCODE-rE2G + PoPS
- **ENCODE_E2G_PoPS_match_CS.R**: For each trait,  the number of GWAS fine-mapping credible sets
that are mapped to genes by ENCODE-rE2G + PoPS
- **Variant_PoPS_1KG.R**: mapping all variants to genes with highest PoPS score among proximal
genes. 
- **GTEx_ENCODE_E2G_preds.R**: mapping GTEx fine-mapped eQTLs to genes using ENCODE-rE2G and then
measuring enrichment 

## Variant enrichment heatmap
Following scripts (with example output) are used to compute the variant enrichment shown as a
heatmap in Figure 4d:

**Step 1. ENCODE_E2G_1KG.R**: Create a variant-level linked gene map with scores for ENCODE-E2G for
each biosample (cell type) - 1000G variants:
```
seqnames	start	end	width	strand	snp_name	seqnames	start	end	width	strand	gene	score Biosample BP_hg38
chr10		134388	134390	3	*	rs182548202	chr10		134080	135557	1478	*	ZMYND11	0.999986273564113 encode_e2g_predictions_8988T_ENCSR000EID_thresholded_predictions.tsv.gz 134389
chr10		134946	134948	3	*	rs146152255	chr10		134080	135557	1478	*	ZMYND11	0.999986273564113 encode_e2g_predictions_8988T_ENCSR000EID_thresholded_predictions.tsv.gz 134947
```

**Step 2. Variant_PoPS_1KG.R**: Create a variant-level linked gene map for each of the UKBB traits
(PoPS approach): includes PoPS score and rank of genes in 1Mb around variant (we do this for all
1000G variants and onlyt use gene-level PoPS scores of linked genes):
```
seqnames	     start    end    width	strand	snp_name	seqnames	start	end	width	strand	pops_score pops_rank gene ensembl BP_hg19
chr1		     11007    11009  3		*	rs575272151	chr1		-430909	570008	1000918	*	-0.384677940707002   15875	  OR4F5 ENSG00000186092 11008
chr1		     11007    11009  3		*	rs575272151	chr1		-132360	868634	1000995	*	-0.232227604026806   12884	  OR4F29		ENSG00000235249 11008
```

**Step 3. ENCODE_E2G_top_two_genes.R**: Obtain top two ENCODE-E2G genes for each variant in each
ENCODE bio-sample:
```
SNP	Gene1	Gene2
rs1001179	CAT	NA
rs1002226	KCNJ11	NCR3LG1
rs1003483	IGF2	NA
```
NA corresponds to the fact that no other gene crossed the threshold set for the ENCODE-E2G
predictions.

**Step 4. PoPS_top_two_genes.R**:	Obtain top two PoPS scored genes for each variant within 1Mb
region around the variant:
```
SNP    Gene1  Gene2
rs575272151   OR4F29	OR4F5
rs544419019   OR4F29	OR4F5
rs540538026   OR4F29	OR4F5
rs62635286    OR4F29	OR4F5
```

**Step 5. Enrichment_calc_GWAS.R & Recall_calc_GWAS.R**: Run enrichment and recall analysis
for ENCODE-E2G+PoPS (using top 2 of each method and seeing if there is overlap between two).

**Step 6. process_heatmap_matrix.R**: Processing the matrix to plot for heatmap based on the
enrichment and recall scores computed in the previous step.
