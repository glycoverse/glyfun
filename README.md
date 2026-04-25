
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyfun <a href="https://glycoverse.github.io/glyfun/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe
version](https://glycoverse.r-universe.dev/glyfun/badges/version)](https://glycoverse.r-universe.dev/glyfun)
[![R-CMD-check](https://github.com/glycoverse/glyfun/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glyfun/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glyfun/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glyfun)
<!-- badges: end -->

`glyfun` provides functional enrichment analysis for glycoproteomics
data. It provides two sets of functions, answering different questions:

-   `enrich_xxx()`: “Which functions are enriched for proteins with
    dysregulated glycosylation?”
-   `enrich_gc_xxx()`: “Which functions are enriched for dyregulated
    glycan traits?”

`gc` is for “glycan-centric”, which means connecting glycan traits
(e.g. core-fucosylation) with functional annotations.

Both Over Representation Analysis (ORA) and Gene Set Enrichment Analysis
(GSEA) are supported, with common databases including GO, KEGG,
Reactome, WikiPathways, DO (Disease Ontology), and NCG (Network of
Cancer Genes).

## Installation

The glyfun package is not a core glycoverse package. You need to install
it individually even if you have installed the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse).

You can install the latest release of glyfun from
[r-universe](https://glycoverse.r-universe.dev/glyfun)
(**recommended**):

``` r
# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::pkg_install("glyfun")
```

Or from [GitHub](https://github.com/glycoverse/glyfun):

``` r
pak::pkg_install("glycoverse/glyfun@*release")
```

Or install the development version (NOT recommended):

``` r
pak::pkg_install("glycoverse/glyfun")
```

**Note:** Tips and troubleshooting for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) are also
applicable here: [Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

## Documentation

-   🚀 Get started:
    [Here](https://glycoverse.github.io/glyfun/articles/glyfun.html)
-   📚 Reference:
    [Here](https://glycoverse.github.io/glyfun/reference/index.html)

## Role in `glycoverse`

`glyfun` integrates seamlessly with `glydet` (for calculating
site-specific derived traits and motif quantification) and `glystats`
(for differential expression analysis). As the final piece of the
glycoproteomics analytical puzzle, `glyfun` bridges the gap between raw
glycan quantification and biological insight, enabling researchers to
move beyond data processing to functional glycomining.

## Example

``` r
library(glycoverse)
library(glyfun)

exp <- auto_clean(real_experiment)

# GO ORA for proteins with dysregulated glycosylation
dea_res <- gly_limma(exp)
enrich_ora_go(dea_res)

# Glycan-centric GO ORA for derived traits
trait_exp <- derive_traits(exp)
trait_dea_res <- gly_limma(trait_exp)
enrich_gc_ora_go(trait_dea_res)
```
