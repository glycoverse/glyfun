# Glycan-Centric GO Gene Set Enrichment Analysis

Performs glycan-centric Gene Ontology (GO) Gene Set Enrichment Analysis
(GSEA). Instead of traditional protein-centric enrichment, this function
links specific glycan traits (e.g., core-fucosylation, sialylation) to
functional annotations. It ranks proteins within each glycan trait
separately, then compares enriched biological functions across traits by
running GSEA for each trait.

## Usage

``` r
enrich_gc_gsea_go(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE
)
```

## Arguments

- dea_res:

  Differential analysis result. Can be one of:

  - Result from
    [`glystats::gly_limma()`](https://glycoverse.github.io/glystats/reference/gly_limma.html)
    (two groups),
    [`glystats::gly_ttest()`](https://glycoverse.github.io/glystats/reference/gly_ttest.html),
    or
    [`glystats::gly_wilcox()`](https://glycoverse.github.io/glystats/reference/gly_wilcox.html),
    called on an
    [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
    of "traitproteomics" type.

  - A tibble with the following columns:

    - `protein`: Uniprot ID of proteins

    - `trait`: A glycosylation trait (e.g. "TFc" for proportion of
      core-fucosylated glycans)

    - `site`: The glycosylation site.

    - `p_val`: p-values, preferably adjusted p-values

    - `log2fc`: log2 of fold change

- rank_by:

  Criteria for ranking proteins. One of the following:

  - "log2fc": log2 fold change with signs

  - "abs_log2fc": absolute log2 fold change

  - "log10p": negative log10 p-value

  - "signed_log10p" (default): log10 p-value with signs of log2 fold
    change

  - "log2fc_log10p": log2 fold change multiplied by negative log10
    p-value

- aggr:

  Aggregation method for combining multiple scores across different
  traits and sites for the same protein. One of "median", "mean", or
  "max". Defaults to "median".

- orgdb:

  An OrgDb object. Passed to `OrgDb` of downstream enrichment function.

- ont:

  Ontology type. Passed to `ont` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
  "BP", "MF", "CC", or "ALL". Defaults to "MF".

- p_adj_method:

  P-value adjustment method. One of "holm", "hochberg", "hommel",
  "bonferroni", "BH", "BY", "fdr", "none". Passed to `pAdjustMethod` of
  downstream enrichment function. Defaults to "BH".

- p_cutoff:

  P-value cutoff to filter significant terms. Passed to `pvalueCutoff`
  of downstream enrichment function. Defaults to 0.05.

- min_gs_size:

  Minimal size of each gene set for analyzing. Gene sets with fewer
  genes than this threshold will be excluded. Passed to `minGSSize` of
  downstream enrichment function. Defaults to 10.

- max_gs_size:

  Maximum size of each gene set for analyzing. Gene sets with more genes
  than this threshold will be excluded. Passed to `maxGSSize` of
  downstream enrichment function. Defaults to 500.

- exponent:

  Weight of each step. Passed to `exponent` of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html).
  Defaults to 1.

- eps:

  Epsilon for calculating p-values. Passed to `eps` of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html).
  Defaults to 1e-10.

- seed:

  Logical indicating whether to set a random seed for reproducibility.
  Passed to `seed` of
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html).
  Defaults to `FALSE`.

## Value

A clusterProfiler `compareClusterResult` object. It can be readily
converted to a tibble with
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
or visualized with `clusterProfiler` functions like
[`clusterProfiler::dotplot()`](https://rdrr.io/pkg/clusterProfiler/man/reexports.html).

## What is glycan-centric enrichment?

In traditional glycoproteomics data analysis, we usually perform
differential expression analysis (DEA) on glycoforms, extract proteins
that have dysregulated glycosylation, then perform functional enrichment
(e.g. GO) on these proteins. This is what `enrich_xxx()` functions do
(e.g.
[`enrich_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gsea_go.md)).

`enrich_gc_xxx()` functions differ in that they link specific glycan
traits with functional annotations. Instead of answering the question
"Which functions are enriched in dysregulated glycoproteins?",
`enrich_gc_xxx()` answers questions like "Which functions are enriched
in proteins ranked highly for core-fucosylation changes?" Higher
specificity, deeper insights. By focusing on distinct glycan motifs, it
helps you pinpoint the functional relevance of specific glycosylation
changes.

## How it ranks proteins

GSEA requires a ranked list of proteins as input. This function first
splits the DEA result by glycan trait, then ranks proteins within each
trait separately. For each trait, it applies the same ranking and
aggregation logic as
[`enrich_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gsea_go.md),
producing one ranked protein list per trait. Those trait-specific ranked
lists are then compared with
[`clusterProfiler::compareCluster()`](https://rdrr.io/pkg/clusterProfiler/man/compareCluster.html)
using its formula interface.

## Common usage pattern

A common pattern of using this function is:

    # 1. Use `glydet` to calculate derived traits or motif quantification.
    trait_exp <- derive_traits(exp)  # or `quantify_motifs()`

    # 2. Perform differential analysis with `glystats`.
    dea_res <- gly_ttest(trait_exp)

    # 3. Use this function.
    go_res <- enrich_gc_gsea_go(dea_res)  # or other `enrich_gc_gsea_xxx()` functions

## See also

[`clusterProfiler::compareCluster()`](https://rdrr.io/pkg/clusterProfiler/man/compareCluster.html),
[`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
