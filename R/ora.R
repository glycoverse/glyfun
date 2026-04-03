#' GO Over Representation Analysis
#'
#' @description
#' Performs Gene Ontology (GO) Over-Representation Analysis (ORA)
#' on glycoproteins with dysregulated glycosylation.
#'
#' @details
#' # Common usage pattern
#'
#' A common pattern of using this function is:
#'
#' ```r
#' # 1. Perform differential analysis with `glystats`.
#' dea_res <- gly_ttest(exp)
#'
#' # 2. Use this function.
#' go_res <- enrich_gc_ora_go(dea_res)  # or other glyfun functions
#' ```
#'
#' @param dea_res Differential analysis result. Can be one of:
#'   - Result from [glystats::gly_limma()] (two groups), [glystats::gly_ttest()], or [glystats::gly_wilcox()],
#'     called on an [glyexp::experiment()] of "traitproteomics" type.
#'   - A tibble with the following columns:
#'     - `protein`: Uniprot ID of proteins
#'     - `trait`: A glycosylation trait (e.g. "TFc" for proportion of core-fucosylated glycans)
#'     - `site`: The glycosylation site.
#'     - `p_val`: p-values, preferably adjusted p-values
#'     - `log2FC`: log2 of fold change
#' @param dea_p_cutoff P-value cutoff for statistical significance. Defaults to 0.05.
#'   For `glystats` result input, adjusted p-values are used.
#' @param dea_log2fc_cutoff Log2 fold change cutoff statistical significance.
#'   A length-2 numeric vector, being negative and positive boundaries, respectively.
#'   For example, `c(-1, 1)` means "log2FC < -1 or log2FC > 1", and `c(-Inf, 1)` means "log2FC > 1".
#'   Defaults to `c(-1, 1)`.
#' @param orgdb Passed to `OrgDb` of [clusterProfiler::enrichGO()].
#' @param ont Passed to `ont` of [clusterProfiler::enrichGO()]. "BP", "MF", "CC", or "ALL". Defaults to "MF".
#' @param universe Background genes. If a character vector, directly passed to `universe` of [clusterProfiler::enrichGO()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   [clusterProfiler::enrichGO()].
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichGO()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichGO()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichGO()].
#'
#' #' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
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
#'  - `raw_result`: The raw clusterProfiler `enrichResult`` object
#' The list has classes `glyfun_gc_ora_go_res`, `glyfun_gc_ora_res`, and `glyfun_res`.
#'
#' @seealso [clusterProfiler::enrichGO()]
#' @export
enrich_ora_go <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  orgdb = "org.Hs.eg.db",
  ont = "MF",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  .check_dea_res_type(dea_res)
  .ora(
    dea_res,
    enrich_fun = clusterProfiler::enrichGO,
    result_class = "glyfun_ora_go_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    keyType = "UNIPROT",
    OrgDb = orgdb,
    ont = ont,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff
  )
}

#' KEGG Over Representation Analysis
#'
#' @description
#' Performs KEGG pathway Over-Representation Analysis (ORA)
#' on glycoproteins with dysregulated glycosylation.
#'
#' @inheritSection enrich_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_go
#' @param organism KEGG organism code. Passed to `organism` of [clusterProfiler::enrichKEGG()].
#'   Defaults to "hsa" (Homo sapiens). Common codes: "hsa" (human), "mmu" (mouse), "rno" (rat).
#' @param universe Background genes. If a character vector, directly passed to `universe` of [clusterProfiler::enrichKEGG()].
#'   You can also provide a [glyexp::experiment()] object with "glycoproteomics" type.
#'   In this case all detected proteins in this experiment will be extracted and passed to
#'   [clusterProfiler::enrichKEGG()].
#' @param p_adj_method Passed to `pAdjustMethod` of [clusterProfiler::enrichKEGG()].
#' @param p_cutoff Passed to `pvalueCutoff` of [clusterProfiler::enrichKEGG()].
#' @param q_cutoff Passed to `qvalueCutoff` of [clusterProfiler::enrichKEGG()].
#'
#' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `id`: KEGG pathway ID
#'    - `description`: Pathway description
#'    - `gene_ratio`: Ratio of genes in the pathway to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the pathway to total genes in the background
#'    - `rich_factor`: Proportion of the pathway's total background genes found in the input
#'    - `fold_enrichment`: Ratio of `gene_ratio` to `bg_ratio` (magnitude of enrichment)
#'    - `z_score`: Directional trend of regulation (positive for up, negative for down)
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_val`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the pathway (separated by "/")
#'    - `count`: Number of genes in the pathway
#'  - `raw_result`: The raw clusterProfiler `enrichResult` object
#' The list has classes `glyfun_ora_kegg_res`, `glyfun_ora_res`, and `glyfun_res`.
#'
#' @seealso [clusterProfiler::enrichKEGG()]
#' @export
enrich_ora_kegg <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  .check_dea_res_type(dea_res)
  .ora(
    dea_res,
    enrich_fun = clusterProfiler::enrichKEGG,
    result_class = "glyfun_ora_kegg_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    keyType = "uniprot",
    organism = organism,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff
  )
}

#' Perform ORA
#' @param dea_res DEA result from glystats.
#' @param enrich_fun An enrichment function.
#' @param result_class A string of the concrete result class.
#' @param dea_p_cutoff P-value cutoff to define statistical significance.
#' @param dea_log2fc_cutoff Log2FC cutoffs to define statistical significance.
#' @param ... Parameters passed to downstream enrichment package.
.ora <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
) {
  .check_p_cutoff_arg(dea_p_cutoff)
  .check_log2fc_cutoff_arg(dea_log2fc_cutoff)

  proteins <- dea_res |>
    glystats::get_tidy_result() |>
    dplyr::filter(
      .data$p_adj < dea_p_cutoff,
      .data$log2fc < dea_log2fc_cutoff[[1]] |
        .data$log2fc > dea_log2fc_cutoff[[2]]
    ) |>
    dplyr::pull(.data$protein)

  suppressWarnings(res <- rlang::exec(enrich_fun, proteins, ...))

  if (is.null(res)) {
    cli::cli_alert_warning("No terms were enriched. `NULL` will be returned.")
    return(NULL)
  }
  tidy_res <- tibble::as_tibble(res) |>
    janitor::clean_names() |>
    dplyr::rename(tidyselect::all_of(c(
      "p_val" = "pvalue",
      "p_adj" = "p_adjust",
      "q_val" = "qvalue"
    )))
  res <- list(tidy_result = tidy_res, raw_result = res)
  structure(res, class = c(result_class, "glyfun_ora_res", "glyfun_res"))
}
