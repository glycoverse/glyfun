# Glycan-Centric Disease Ontology (DO) Over Representation Analysis

Performs glycan-centric Disease Ontology (DO) Over-Representation
Analysis (ORA). Instead of traditional protein-centric enrichment, this
function links specific glycan traits to disease associations. It helps
answer questions like "Which diseases are enriched in proteins with a
specific dysregulated glycan motif?", by grouping differential analysis
results by glycan traits and computing disease enrichment for each
trait.

## Usage

``` r
enrich_gc_ora_do(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ont = "HDO",
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
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

- dea_p_cutoff:

  P-value cutoff for statistical significance. Defaults to 0.05. For
  `glystats` result input, adjusted p-values are used.

- dea_log2fc_cutoff:

  Log2 fold change cutoff statistical significance. A length-2 numeric
  vector, being negative and positive boundaries, respectively. For
  example, `c(-1, 1)` means "log2fc \< -1 or log2fc \> 1", and
  `c(-Inf, 1)` means "log2fc \> 1". Defaults to `c(-1, 1)`.

- ont:

  One of "HDO" (Human Disease Ontology), "MPO" (Mammalian Phenotype
  Ontology), or "VDO" (Vector Disease Ontology). Passed to `ont` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).
  Defaults to "HDO".

- organism:

  "hsa" (Homo sapiens) or "mmu" (Mus musculus). Passed to `organism` of
  [`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html).
  Defaults to "hsa".

- universe:

  Background genes Uniprot IDs, directly passed to `universe` of
  downstream enrichment function. If `NULL` (default), all genes in the
  data will be used. Another common pattern is to use all detected
  proteins as backgroud genes. You can use
  [`detected_universe()`](https://glycoverse.github.io/glyfun/dev/reference/detected_universe.md)
  to help you.

- p_adj_method:

  P-value adjustment method. One of "holm", "hochberg", "hommel",
  "bonferroni", "BH", "BY", "fdr", "none". Passed to `pAdjustMethod` of
  downstream enrichment function. Defaults to "BH".

- p_cutoff:

  P-value cutoff to filter significant terms. Passed to `pvalueCutoff`
  of downstream enrichment function. Defaults to 0.05.

- q_cutoff:

  Q-value (FDR) cutoff to filter significant terms. Passed to
  `qvalueCutoff` of downstream enrichment function. Defaults to 0.2.

- min_gs_size:

  Minimal size of each gene set for analyzing. Gene sets with fewer
  genes than this threshold will be excluded. Passed to `minGSSize` of
  downstream enrichment function. Defaults to 10.

- max_gs_size:

  Maximum size of each gene set for analyzing. Gene sets with more genes
  than this threshold will be excluded. Passed to `maxGSSize` of
  downstream enrichment function. Defaults to 500.

## Value

A clusterProfiler `compareClusterResult` object with additional `glyfun`
classes. It can be readily converted to a tibble with
[`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html),
or visualized with `clusterProfiler` functions like
[`clusterProfiler::dotplot()`](https://rdrr.io/pkg/clusterProfiler/man/reexports.html).

## What is glycan-centric enrichment?

In traditional glycoproteomics data analysis, we usually perform
differential expression analysis (DEA) on glycoforms, extract proteins
that have dysregulated glycosylation, then perform functional enrichment
(e.g. GO) on these proteins. This is what `enrich_xxx()` functions do
(e.g.
[`enrich_ora_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_ora_go.md)).

`enrich_gc_xxx()` functions differ in that they link specific glycan
traits with functional annotations. Instead of answering the question
"Which functions are enriched in dysregulated glycoproteins?",
`enrich_gc_xxx()` answers questions like “Which functions are enriched
in proteins with dysregulated core-fucosylation?” Higher specificity,
deeper insights. By focusing on distinct glycan motifs, it helps you
pinpoint the functional relevance of specific glycosylation changes.

## Common usage pattern

A common pattern of using this function is:

    # 1. Use `glydet` to calculate derived traits or motif quantification.
    trait_exp <- derive_traits(exp)  # or `quantify_motifs()`

    # 2. Perform differential analysis with `glystats`.
    dea_res <- gly_ttest(trait_exp)

    # 3. Use this function.
    go_res <- enrich_gc_ora_go(dea_res)  # or other `enrich_gc_xxx()` functions

## See also

[`clusterProfiler::compareCluster()`](https://rdrr.io/pkg/clusterProfiler/man/compareCluster.html),
[`DOSE::enrichDO()`](https://rdrr.io/pkg/DOSE/man/enrichDO.html)
