# Glycan-Centric Reactome Pathway Over Representation Analysis

Performs glycan-centric Reactome pathway Over-Representation Analysis
(ORA). Instead of traditional protein-centric enrichment, this function
links specific glycan traits to biological pathways. It helps answer
questions like "Which Reactome pathways are enriched in proteins with a
specific dysregulated glycan motif?", by grouping differential analysis
results by glycan traits and computing pathway enrichment for each
trait.

## Usage

``` r
enrich_gc_ora_reactome(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  orgdb = "org.Hs.eg.db",
  organism = "human",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
)
```

## Value

A list with two elements:

- `tidy_result`: A tibble with enrichment results containing the
  following columns:

  - `trait`: Glycan trait

  - `id`: Reactome pathway ID

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

- `raw_result`: The raw clusterProfiler clusterProfResult object The
  list has classes `glyfun_gc_ora_reactome_res`, `glyfun_gc_ora_res`,
  and `glyfun_res`.

## What is glycan-centric enrichment?

In traditional glycoproteomics data analysis, we usually perform
differential expression analysis (DEA) on glycoforms, extract proteins
that have dysregulated glycosylation, then perform functional enrichment
(e.g. GO) on these proteins. This is what `enrich_xxx()` functions do
(e.g.
[`enrich_ora_go()`](https://glycoverse.github.io/glyfun/reference/enrich_ora_go.md)).

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
[`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html)
