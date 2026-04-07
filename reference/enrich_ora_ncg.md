# Network of Cancer Genes (NCG) Over Representation Analysis

Performs Network of Cancer Genes (NCG) Over-Representation Analysis
(ORA) on glycoproteins with dysregulated glycosylation.

## Usage

``` r
enrich_ora_ncg(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
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

    - `log2FC`: log2 of fold change

- dea_p_cutoff:

  P-value cutoff for statistical significance. Defaults to 0.05. For
  `glystats` result input, adjusted p-values are used.

- dea_log2fc_cutoff:

  Log2 fold change cutoff statistical significance. A length-2 numeric
  vector, being negative and positive boundaries, respectively. For
  example, `c(-1, 1)` means "log2FC \< -1 or log2FC \> 1", and
  `c(-Inf, 1)` means "log2FC \> 1". Defaults to `c(-1, 1)`.

- universe:

  Background genes Uniprot IDs, directly passed to `universe` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
  If `NULL` (default), all genes in the data will be used. Another
  common pattern is to use all detected proteins as backgroud genes. You
  can use
  [`detected_universe()`](https://glycoverse.github.io/glyfun/reference/detected_universe.md)
  to help you.

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

## Value

A list with two elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: NCG cancer gene set ID

  - `description`: Cancer type or gene set description

  - `gene_ratio`: Ratio of genes in the set to total genes in the input

  - `bg_ratio`: Ratio of genes in the set to total genes in the
    background

  - `rich_factor`: Proportion of the set's total background genes found
    in the input

  - `fold_enrichment`: Ratio of `gene_ratio` to `bg_ratio` (magnitude of
    enrichment)

  - `z_score`: Directional trend of regulation (positive for up,
    negative for down)

  - `p_val`: Raw p-value from hypergeometric test

  - `p_adj`: Adjusted p-value

  - `q_val`: Q-value (FDR)

  - `gene_id`: Gene IDs in the set (separated by "/")

  - `count`: Number of genes in the set

- `raw_result`: The raw clusterProfiler `enrichResult` object The list
  has classes `glyfun_ora_ncg_res`, `glyfun_ora_res`, and `glyfun_res`.

## Common usage pattern

A common pattern of using this function is:

    # 1. Perform differential analysis with `glystats`.
    dea_res <- gly_ttest(exp)

    # 2. Use this function.
    go_res <- enrich_gc_ora_go(dea_res)  # or other glyfun functions

## See also

[`DOSE::enrichNCG()`](https://rdrr.io/pkg/DOSE/man/enrichNCG.html)
