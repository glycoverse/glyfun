# Unified functional enrichment for glycoproteomics

Runs protein-centric or glycan-trait-centric over-representation
analysis (ORA) or gene set enrichment analysis (GSEA) through one entry
point.

Protein-centric enrichment analyzes proteins with dysregulated
glycosylation. Glycan-trait-centric enrichment analyzes protein sets or
protein rankings stratified by a glycan trait; it does not directly
enrich glycan structures.

## Usage

``` r
enrich_glyco(
  dea_res,
  method = c("ora", "gsea"),
  database = c("GO", "KEGG", "Reactome", "WikiPathways", "DO", "NCG"),
  scope = c("protein", "glycan_trait"),
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  rank_by = "signed_log10p",
  aggr = "median",
  orgdb = NULL,
  ont = NULL,
  organism = NULL,
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE,
  ...,
  by = NULL
)
```

## Arguments

- dea_res:

  Differential analysis result. It can be a supported `glystats` result
  or a data frame accepted by the selected existing `enrich_*()` family.

- method:

  Enrichment method: `"ora"` or `"gsea"`.

- database:

  Annotation database: `"GO"`, `"KEGG"`, `"Reactome"`, `"WikiPathways"`,
  `"DO"`, or `"NCG"`.

- scope:

  Analysis scope: `"protein"` for protein-centric enrichment or
  `"glycan_trait"` for glycan-trait-centric enrichment.

- dea_p_cutoff:

  P-value cutoff used to select proteins for ORA.

- dea_log2fc_cutoff:

  Log2 fold-change cutoffs used to select proteins for ORA. The value
  must contain a negative and a positive boundary.

- rank_by:

  Ranking criterion used for GSEA. See
  [`enrich_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gsea_go.md).

- aggr:

  Aggregation method used to combine scores for the same protein in
  GSEA. See
  [`enrich_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gsea_go.md).

- orgdb:

  An `OrgDb` object or its package name. It defaults to `"org.Hs.eg.db"`
  for GO and NCG.

- ont:

  Ontology used by GO or DO. It defaults to `"MF"` for GO and `"HDO"`
  for DO.

- organism:

  Organism passed to the selected database. Defaults are `"hsa"` for
  KEGG and DO, `"human"` for Reactome, and `"Homo sapiens"` for
  WikiPathways.

- universe:

  Background proteins for ORA.

- p_adj_method:

  P-value adjustment method passed to the enrichment backend.

- p_cutoff:

  P-value cutoff passed to the enrichment backend.

- q_cutoff:

  Q-value cutoff passed to the enrichment backend for ORA.

- min_gs_size:

  Minimum gene-set size.

- max_gs_size:

  Maximum gene-set size.

- exponent:

  GSEA weighting exponent.

- eps:

  Epsilon used for GSEA p-value calculation.

- seed:

  Whether to set a random seed in GSEA.

- ...:

  Additional arguments passed to the selected enrichment backend.

- by:

  Glycan-trait grouping column when `scope = "glycan_trait"`. See
  [`enrich_gc_ora_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gc_ora_go.md)
  for accepted values and automatic inference. Must be `NULL` for
  `scope = "protein"`.

## Value

A clusterProfiler enrichment result. ORA returns an `enrichResult` or
`compareClusterResult`, and GSEA returns a `gseaResult` or
`compareClusterResult`.

## See also

[`enrich_ora_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_ora_go.md),
[`enrich_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gsea_go.md),
[`enrich_gc_ora_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gc_ora_go.md),
[`enrich_gc_gsea_go()`](https://glycoverse.github.io/glyfun/dev/reference/enrich_gc_gsea_go.md)

## Examples

``` r
if (FALSE) { # \dontrun{
enrich_glyco(dea_res, method = "ora", database = "GO")
enrich_glyco(
  dea_res,
  method = "gsea",
  database = "Reactome",
  scope = "glycan_trait"
)
} # }
```
