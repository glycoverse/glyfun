# glyfun 0.2.0

## New features

* `enrich_glyco()` provides a unified entry point for protein-centric and glycan-trait-centric ORA and GSEA across GO, KEGG, Reactome, WikiPathways, DO, and NCG. (#5)
* Glycan-centric ORA and GSEA functions gain `by` to select the glycan-trait grouping column, including through `enrich_glyco(scope = "glycan_trait")`; automatic selection preserves existing defaults and reports its choice when multiple candidates are present. (#3)

## Minor improvements and bug fixes

* ORA and effect-dependent GSEA reject unavailable log2 fold changes instead of silently producing empty input; significance-only GSEA accepts results without log2 fold changes and diagnoses invalid or excluded scores. (#8)
* Protein-centric GSEA functions return `NULL` with a diagnostic message instead of failing when an enrichment backend returns no result.

# glyfun 0.1.3

* Documentation and vignettes now use `GlycoproteomicSE` inputs and `SummarizedExperiment` accessors for Stage II of glycoverse/glyexp#15. (#4)

# glyfun 0.1.2

* The glyfun vignette now runs with both legacy `glyexp::experiment()` and current `SummarizedExperiment` containers.

# glyfun 0.1.1

* `detected_universe()` now accepts `glyexp::GlycoproteomicSE` objects while retaining support for `glyexp::experiment()` objects. (#1)

# glyfun 0.1.0

* First release.
