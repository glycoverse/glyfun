# Helper function to prepare the `universe` parameter

This function extracts all detected proteins in a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
or a `glystats` result. It can be readily passed to the `universe`
parameter of all `glyfun` functions.

## Usage

``` r
detected_universe(x)
```

## Arguments

- x:

  A
  [`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
  or a `glystats` result.

## Value

A character vector of protein UniProt IDs.

## Examples

``` r
library(glyexp)
universe <- detected_universe(real_experiment)
length(universe)
#> [1] 4262
universe[1:5]
#> [1] "P08185" "P04196" "P04196" "P04196" "P10909"
```
