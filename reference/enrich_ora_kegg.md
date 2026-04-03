# KEGG Over Representation Analysis

Performs KEGG pathway Over-Representation Analysis (ORA) on
glycoproteins with dysregulated glycosylation.

## Usage

``` r
enrich_ora_kegg(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  organism = "hsa",
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

- organism:

  KEGG organism code. Passed to `organism` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).
  Defaults to "hsa" (Homo sapiens). Common codes: "hsa" (human), "mmu"
  (mouse), "rno" (rat).

- universe:

  Background genes. If a character vector, directly passed to `universe`
  of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).
  You can also provide a
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  object with "glycoproteomics" type. In this case all detected proteins
  in this experiment will be extracted and passed to
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- p_adj_method:

  Passed to `pAdjustMethod` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- p_cutoff:

  Passed to `pvalueCutoff` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

- q_cutoff:

  Passed to `qvalueCutoff` of
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html).

## Value

A list with two elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `id`: KEGG pathway ID

  - `description`: Pathway description

  - `gene_ratio`: Ratio of genes in the pathway to total genes in the
    input

  - `bg_ratio`: Ratio of genes in the pathway to total genes in the
    background

  - `rich_factor`: Proportion of the pathway's total background genes
    found in the input

  - `fold_enrichment`: Ratio of `gene_ratio` to `bg_ratio` (magnitude of
    enrichment)

  - `z_score`: Directional trend of regulation (positive for up,
    negative for down)

  - `p_val`: Raw p-value from hypergeometric test

  - `p_adj`: Adjusted p-value

  - `q_val`: Q-value (FDR)

  - `gene_id`: Gene IDs in the pathway (separated by "/")

  - `count`: Number of genes in the pathway

- `raw_result`: The raw clusterProfiler `enrichResult` object The list
  has classes `glyfun_ora_kegg_res`, `glyfun_ora_res`, and `glyfun_res`.

## Common usage pattern

A common pattern of using this function is:

    # 1. Perform differential analysis with `glystats`.
    dea_res <- gly_ttest(exp)

    # 2. Use this function.
    go_res <- enrich_gc_ora_go(dea_res)  # or other glyfun functions

## See also

[`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html)
