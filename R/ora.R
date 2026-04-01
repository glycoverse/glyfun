#' Glycan-Centric GO Over Representation Analysis
#'
#' @description
#' This function first groups the proteins in `dea_res` according to `by`,
#' then performs GO ORA analysis with `clusterProfiler::compareCluster()` and
#' `clusterProfiler::enrichGO()`
#'
#' @details
#' # What is glycan-centric enrichment?
#'
#' In traditional glycoproteomics data analysis,
#' we usually perform differential expression analysis (DEA) on glycoforms,
#' extract proteins that have dysregulated glycosylation,
#' then perform functional enrichment (e.g. GO) on these proteins.
#' This is what all the enrichment functions in glystats do (e.g. `glystats::gly_enrich_go()`).
#'
#' `glyfun` functions differ in that they link specific glycan traits with functional annotations.
#' Instead of answering the question
#' "Which functions are enriched in dysregulated glycoproteins?",
#' `glyfun` answers questions like
#' “Which functions are enriched in proteins with dysregulated core-fucosylation?”
#' Higher specificity, deeper insights. By focusing on distinct glycan motifs,
#' glyfun helps you pinpoint the functional relevance of specific glycosylation changes.
#'
#' # Common usage pattern
#'
#' A common pattern of using this function is:
#'
#' ```r
#' # 1. Use `glydet` to calculate derived traits or motif quantification.
#' trait_exp <- derive_traits(exp)  # or `quantify_motifs()`
#'
#' # 2. Perform differential analysis with `glystats`.
#' dea_res <- gly_ttest(trait_exp)
#'
#' # 3. Use this function.
#' go_res <- gc_ora_go(dea_res)
#' ```
#'
#' @param dea_res Differential analysis result. Can be one of:
#'   - Result from [glystats::gly_limma()] (two groups), [glystats::gly_ttest()], or [glystats::gly_wilcox()],
#'     called on an [glyexp::experiment()] of "traitproteomics" type.
#'   - A tibble with the following columns:
#'     - `protein`: Uniprot ID of proteins
#'     - `trait`: A glycosylation trait (e.g. "TFc" for proportion of core-fucosylated glycans)
#'     - `p_val`: p-values, preferably adjusted p-values
#'     - `log2FC`: log2 of fold change
#' @param by A column to group the proteins by.
#'   - If `dea_res` is a [glyexp::experiment()]: the column name in `var_info` of the experiment.
#'   - If `dea_res` is a tibble: the column name in the tibble (defaults to "trait").
#' @param dea_p_cutoff P-value cutoff for statistical significance. Defaults to 0.05.
#'   For `glystats` result input, adjusted p-values are used.
#' @param dea_log2fc_cutoff Log2 fold change cutoff statistical significance.
#'   A length-2 numeric vector, being negative and positive boundaries, respectively.
#'   For example, `c(-1, 1)` means "log2FC < -1 or log2FC > 1", and `c(-Inf, 1)` means "log2FC > 1".
#'   Defaults to `c(-1, 1)`.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::enrichGO()].
#' @param keytype Passed to `keyType` of [clusterProfiler::enrichGO()].
#' @param ont Passed to `ont` of [clusterProfiler::enrichGO()]. "BP", "MF", "CC", or "ALL". Defaults to "MF".
#' @param universe Background genes. If a character vector, directly passed to `universe` of [clusterProfiler::enrichGO()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   `universe` of [clusterProfiler::enrichGO()].
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichGO()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichGO()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichGO()].
#'
#' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `trait`: Glycan trait
#'    - `id`: Term ID
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `rich_factor`: Proportion of the term's total background genes found in the input
#'    - `fold_enrichment`: Ratio of `gene_ratio` to `bg_ratio` (magnitude of enrichment)
#'    - `z_score`: Directional trend of regulation (positive for up, negative for down)
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_val`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler clusterProfResult object
#' The list has classes `glystats_go_ora_res` and `glystats_res`.
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::enrichGO()]
#' @export
gc_ora_go <- function(
  dea_res,
  by = NULL,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  basic_class <- class(dea_res)[[1]]
  supported_classes <- c("glystats_limma_res", "glystats_ttest_res", "glystats_wilcox_res")
  if (!basic_class %in% supported_classes) {
    cli::cli_abort(c(
      "Unsupported input class for {.arg dea_res}.",
      "i" = "Expected: {.cls {supported_classes}}",
      "x" = "Got: {.cls {basic_class}}"
    ))
  }
  .gc_ora(
    dea_res,
    enrich_fun = "enrichGO",
    result_class = "glyfun_ora_go_res",
    by = by,
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    orgdb = orgdb,
    ont = ont,
    universe = universe,
    p_adj_method = p_adj_method,
    p_cutoff = p_cutoff,
    q_cutoff = q_cutoff
  )
}

.gc_ora <- function(
  dea_res,
  enrich_fun,
  result_class,
  by = NULL,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
) {
  by <- .process_by_arg_glystats(dea_res, by)
  .check_p_cutoff_arg(dea_p_cutoff)
  .check_log2fc_cutoff_arg(dea_log2fc_cutoff)

  protein_list <- dea_res |>
    glystats::get_tidy_result() |>
    dplyr::filter(
      .data$p_adj < dea_p_cutoff,
      .data$log2fc < dea_log2fc_cutoff[[1]] | .data$log2fc > dea_log2fc_cutoff[[2]]
    ) |>
    dplyr::summarise(proteins = list(unique(.data$protein)), .by = tidyselect::all_of(by)) |>
    tibble::deframe()

  n_traits <- length(names(protein_list))
  cli::cli_alert_info("Enriching for {.val {n_traits}} glycan traits... (This can take long)")

  suppressWarnings(
    ck <- .call_compare_cluster(
      protein_list,
      fun = enrich_fun,
      keytype = "UNIPROT",
      ...
    )
  )

  if (is.null(ck)) {
    tidy_res <- tibble::tibble(
      trait = character(),
      id = character(),
      description = character(),
      gene_ratio = character(),
      bg_ratio = character(),
      rich_factor = numeric(),
      fold_enrichment = numeric(),
      z_score = numeric(),
      p_val = numeric(),
      p_adj = numeric(),
      q_val = numeric(),
      gene_id = character(),
      count = integer()
    )
  } else {
    tidy_res <- tibble::as_tibble(ck) |>
      janitor::clean_names() |>
      dplyr::rename(tidyselect::all_of(c(
        "trait" = "cluster",
        "p_val" = "pvalue",
        "p_adj" = "p_adjust",
        "q_val" = "qvalue"
      )))
  }

  res <- list(tidy_result = tidy_res, raw_result = ck)
  structure(res, class = c(result_class, "glyfun_ora_res", "glyfun_res"))
}