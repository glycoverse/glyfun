# Perform ORA

Perform ORA

## Usage

``` r
.ora(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
)
```

## Arguments

- dea_res:

  DEA result from glystats.

- enrich_fun:

  An enrichment function.

- result_class:

  A string of the concrete result class.

- dea_p_cutoff:

  P-value cutoff to define statistical significance.

- dea_log2fc_cutoff:

  Log2FC cutoffs to define statistical significance.

- ...:

  Parameters passed to downstream enrichment package.
