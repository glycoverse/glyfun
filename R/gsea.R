#' GO Gene Set Enrichment Analysis
#'
#' @description
#' Performs Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA)
#' on glycoproteins with dysregulated glycosylation.
#'
#' @details
#' # How it ranks proteins
#'
#' GSEA requires a ranked list of proteins as input.
#' This function ranks proteins based on the median absolute log2 fold change across all traits and sites.
#' This reflects the overall glycosylation dysregulation degree of each glycoprotein.
#' You can use `rank_by` to specify other ranking criteria, such as p-values or signed log2 fold changes.
#' You can also use `aggr` to specify how to aggregate multiple scores for the same protein across different traits and sites.
#'
#' # Common usage pattern
#'
#' A common pattern of using this function is:
#'
#' ```r
#' # 1. Perform differential analysis with `glystats`.
#' dea_res <- gly_ttest(exp)
#'
#' # 2. Use this function.
#' go_res <- enrich_gsea_go(dea_res)  # or `enrich_gsea_xxx()` functions
#' ```
#'
#' @inheritParams enrich_ora_go
#' @param rank_by Criteria for ranking proteins.
#'   One of the following:
#'   - "log2fc": log2 fold change with signs
#'   - "abs_log2fc": absolute log2 fold change
#'   - "log10p": negative log10 p-value
#'   - "signed_log10p" (default): log10 p-value with signs of log2 fold change
#'   - "log2fc_log10p": log2 fold change multiplied by negative log10 p-value
#' @param aggr Aggregation method for combining multiple scores across different traits and sites for the same protein.
#'   One of "median", "mean", or "max". Defaults to "median".
#' @param exponent Weight of each step. Passed to `exponent` of [clusterProfiler::gseGO()].
#'   Defaults to 1.
#' @param eps Epsilon for calculating p-values. Passed to `eps` of [clusterProfiler::gseGO()].
#'   Defaults to 1e-10.
#' @param seed Logical indicating whether to set a random seed for reproducibility.
#'   Passed to `seed` of [clusterProfiler::gseGO()]. Defaults to `TRUE`.
#'
#' @return A clusterProfiler `gseaResult` object.
#'   It can be readily converted to a tibble with [tibble::as_tibble()],
#'   or visualized with `clusterProfiler` functions like [clusterProfiler::ridgeplot()].
#'
#' @seealso [clusterProfiler::gseGO()]
#' @export
enrich_gsea_go <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = TRUE
) {
  orgdb <- .prepare_orgdb(orgdb)
  .gsea(
    dea_res,
    enrich_fun = clusterProfiler::gseGO,
    result_class = "glyfun_gsea_go",
    rank_by = rank_by,
    aggr = aggr,
    OrgDb = orgdb,
    keyType = "UNIPROT",
    ont = ont,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed
  )
}

#' KEGG Gene Set Enrichment Analysis
#'
#' @description
#' Performs KEGG pathway Gene Set Enrichment Analysis (GSEA)
#' on glycoproteins with dysregulated glycosylation.
#'
#' @inheritSection enrich_gsea_go How it ranks proteins
#' @inheritSection enrich_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_go
#' @param organism KEGG organism code. Defaults to "hsa" (Homo sapiens).
#'   See [clusterProfiler::gseKEGG()] for details.
#'
#' @return A clusterProfiler `gseaResult` object.
#'   It can be readily converted to a tibble with [tibble::as_tibble()],
#'   or visualized with `clusterProfiler` functions like [clusterProfiler::ridgeplot()].
#'
#' @seealso [clusterProfiler::gseKEGG()]
#' @export
enrich_gsea_kegg <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  organism = "hsa",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = TRUE
) {
  .gsea(
    dea_res,
    enrich_fun = clusterProfiler::gseKEGG,
    result_class = "glyfun_gsea_kegg",
    rank_by = rank_by,
    aggr = aggr,
    keyType = "uniprot",
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed
  )
}

.gsea <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  UseMethod(".gsea")
}

.gsea.data.frame <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  pro_fun <- function(dea_res) {
    .prepare_pro_list(dea_res, rank_by, aggr)
  }
  .gsea_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_fun = pro_fun
  )
}

.gsea.glystats_res <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  pro_fun <- function(dea_res) {
    dea_res |>
      glystats::get_tidy_result() |>
      .prepare_pro_list(rank_by, aggr)
  }
  .gsea_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_fun = pro_fun
  )
}

#' Prepare a named vector of protein scores for GSEA.
#'
#' This function requires `df` to have a `protein` column and a `score` column.
#' It calculates the median score for each protein, sorts them in descending order,
#' and returns a named vector with protein names as names and median scores as values.
#' @noRd
.prepare_pro_list <- function(df, rank_by, aggr) {
  p_col <- if ("p_adj" %in% colnames(df)) "p_adj" else "p_val"
  scores <- .calcu_rank_scores(df, rank_by)
  aggr_fun <- .gsea_aggr_fun(aggr)
  df |>
    dplyr::mutate(score = scores) |>
    dplyr::summarise(
      score = aggr_fun(.data$score),
      .by = tidyselect::all_of("protein")
    ) |>
    dplyr::arrange(dplyr::desc(.data$score)) |>
    dplyr::select(tidyselect::all_of(c("protein", "score"))) |>
    tibble::deframe()
}

.calcu_rank_scores <- function(df, rank_by) {
  p_col <- if ("p_adj" %in% colnames(df)) "p_adj" else "p_val"
  switch(
    rank_by,
    log2fc = df$log2fc,
    abs_log2fc = abs(df$log2fc),
    log10p = -log10(df[[p_col]]),
    signed_log10p = sign(df$log2fc) * (-log10(df[[p_col]])),
    log2fc_log10p = df$log2fc * (-log10(df[[p_col]])),
    cli::cli_abort("Invalid rank_by method: {.val {rank_by}}")
  )
}

.gsea_aggr_fun <- function(aggr) {
  switch(
    aggr,
    median = stats::median,
    mean = mean,
    max = max,
    cli::cli_abort("Invalid aggr method: {.val {aggr}}")
  )
}

.gsea_impl <- function(
  dea_res,
  enrich_fun,
  result_class,
  bitr_orgdb = NULL,
  ...,
  pro_fun = NULL,
  uniprot_to_entrez = FALSE
) {
  # Argument validation
  if (uniprot_to_entrez && is.null(bitr_orgdb)) {
    cli::cli_abort(
      "{.arg bitr_orgdb} must be provided when {.arg uniprot_to_entrez} is TRUE."
    )
  }
  .check_dea_res(dea_res)

  # Performing enrichment
  proteins <- pro_fun(dea_res)
  if (uniprot_to_entrez) {
    names(proteins) <- .uniprot_to_entrez(
      names(proteins),
      bitr_orgdb,
      drop_na = FALSE
    )
    proteins <- proteins[!is.na(names(proteins))]
  }
  suppressWarnings(
    res <- suppressPackageStartupMessages(
      rlang::exec(enrich_fun, proteins, ...)
    )
  )
  if (nrow(res) == 0) {
    cli::cli_alert_warning("No terms were enriched. `NULL` will be returned.")
    return(NULL)
  }

  res
}
