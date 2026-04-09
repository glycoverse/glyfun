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
#' go_res <- enrich_gsea_go(dea_res)
#' ```
#'
#' @inheritParams enrich_ora_go
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

.gsea <- function(
  dea_res,
  enrich_fun,
  result_class,
  bitr_orgdb = NULL,
  ...
) {
  UseMethod(".gsea")
}

.gsea.data.frame <- function(
  dea_res,
  enrich_fun,
  result_class,
  bitr_orgdb = NULL,
  ...
) {
  .gsea_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_fun = .prepare_pro_list
  )
}

.gsea.glystats_res <- function(
  dea_res,
  enrich_fun,
  result_class,
  bitr_orgdb = NULL,
  ...
) {
  pro_fun <- function(dea_res) {
    dea_res |>
      glystats::get_tidy_result() |>
      .prepare_pro_list()
  }
  .gsea_impl(
    dea_res,
    enrich_fun = clusterProfiler::gseGO,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_fun = pro_fun
  )
}

.prepare_pro_list <- function(df) {
  df |>
    dplyr::summarise(score = median(abs(.data$log2fc)), .by = .data$protein) |>
    dplyr::arrange(dplyr::desc(.data$score)) |>
    dplyr::select(.data$protein, .data$score) |>
    tibble::deframe()
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
    names(proteins) <- .uniprot_to_entrez(names(proteins), bitr_orgdb, drop_na = FALSE)
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