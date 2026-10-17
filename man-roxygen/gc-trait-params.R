#' @param dea_res A supported `glystats` differential analysis result, or a data
#'   frame with `protein` (UniProt IDs), `site`, `p_val`, `log2fc`, and a
#'   glycan-trait grouping column selected by `by`.
#' @param by A single column name in the tidy DEA result to use for glycan-trait
#'   grouping. Any one-dimensional atomic vector (including character or factor
#'   labels) is accepted. With `NULL`, `glystats` results use the first
#'   available column in this order: `glycan_structure`, `glycan_composition`,
#'   `trait`, `motif`. Data frames retain their legacy preference for `trait`,
#'   then use the same order if `trait` is absent. When multiple candidates are
#'   present, a message identifies the selected column and explains how to
#'   override it with `by`.
