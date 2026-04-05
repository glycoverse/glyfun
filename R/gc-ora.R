#' Glycan-Centric GO Over Representation Analysis
#'
#' @description
#' Performs glycan-centric Gene Ontology (GO) Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits (e.g., core-fucosylation, sialylation) to functional annotations.
#' It identifies which biological functions are significantly enriched in proteins
#' exhibiting specific glycosylation changes, grouping the differential analysis
#' results by trait before performing ORA.
#'
#' @details
#' # What is glycan-centric enrichment?
#'
#' In traditional glycoproteomics data analysis,
#' we usually perform differential expression analysis (DEA) on glycoforms,
#' extract proteins that have dysregulated glycosylation,
#' then perform functional enrichment (e.g. GO) on these proteins.
#' This is what `enrich_xxx()` functions do (e.g. [enrich_ora_go()]).
#'
#' `enrich_gc_xxx()` functions differ in that they link specific glycan traits with functional annotations.
#' Instead of answering the question
#' "Which functions are enriched in dysregulated glycoproteins?",
#' `enrich_gc_xxx()` answers questions like
#' “Which functions are enriched in proteins with dysregulated core-fucosylation?”
#' Higher specificity, deeper insights. By focusing on distinct glycan motifs,
#' it helps you pinpoint the functional relevance of specific glycosylation changes.
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
#' go_res <- enrich_gc_ora_go(dea_res)  # or other `enrich_gc_xxx()` functions
#' ```
#'
#' @inheritParams enrich_ora_go
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
#' The list has classes `glyfun_gc_ora_go_res`, `glyfun_gc_ora_res`, and `glyfun_res`.
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::enrichGO()]
#' @export
enrich_gc_ora_go <- function(
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
  .gc_ora(
    dea_res,
    enrich_fun = "enrichGO",
    result_class = "glyfun_gc_ora_go_res",
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

#' Glycan-Centric KEGG Over Representation Analysis
#'
#' @description
#' Performs glycan-centric KEGG pathway Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits to biological pathways. It helps answer questions like "Which
#' pathways are enriched in proteins with a specific dysregulated glycan motif?",
#' by grouping differential analysis results by glycan traits and computing
#' pathway enrichment for each trait.
#'
#' @inheritSection enrich_gc_ora_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_kegg
#' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `trait`: Glycan trait
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
#'  - `raw_result`: The raw clusterProfiler clusterProfResult object
#' The list has classes `glyfun_gc_ora_kegg_res`, `glyfun_gc_ora_res`, and `glyfun_res`.
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::enrichKEGG()]
#' @export
enrich_gc_ora_kegg <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2
) {
  .gc_ora(
    dea_res,
    enrich_fun = "enrichKEGG",
    result_class = "glyfun_gc_ora_kegg_res",
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

#' Perform Glycan-Centric ORA
#' @param dea_res DEA result from glystats or a tibble.
#' @param enrich_fun An enrichment function name (string).
#' @param result_class A string of the concrete result class.
#' @param dea_p_cutoff P-value cutoff to define statistical significance.
#' @param dea_log2fc_cutoff Log2FC cutoffs to define statistical significance.
#' @param ... Parameters passed to downstream enrichment package.
#' @noRd
.gc_ora <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
) {
  UseMethod(".gc_ora")
}

.gc_ora.data.frame <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
) {
  pro_list_fun <- function(dea_res) {
    dea_res |>
      dplyr::filter(
        .data$p_val < dea_p_cutoff,
        .data$log2FC < dea_log2fc_cutoff[[1]] |
          .data$log2FC > dea_log2fc_cutoff[[2]]
      ) |>
      dplyr::summarise(
        proteins = list(unique(.data$protein)),
        .by = "trait"
      ) |>
      tibble::deframe()
  }
  .gc_ora_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    ...,
    pro_list_fun = pro_list_fun
  )
}

.gc_ora.glystats_res <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...
) {
  pro_list_fun <- function(dea_res) {
    tidy_dea_res <- glystats::get_tidy_result(dea_res)
    by <- dplyr::case_when(
      "glycan_structure" %in% colnames(tidy_dea_res) ~ "glycan_structure",
      "glycan_composition" %in% colnames(tidy_dea_res) ~ "glycan_composition",
      "trait" %in% colnames(tidy_dea_res) ~ "trait",
      "motif" %in% colnames(tidy_dea_res) ~ "motif"
    )
    if (is.na(by)) {
      required_cols <- c(
        "glycan_structure",
        "glycan_composition",
        "trait",
        "motif"
      )
      cli::cli_abort(c(
        "Cannot determine glycan traits.",
        "i" = "At least one of these columns is needed: {.field {required_cols}}",
        "i" = "Did you accidentally set {.arg add_info} to `FALSE` when performing DEA with {.pkg glystats}?"
      ))
    }

    tidy_dea_res |>
      dplyr::filter(
        .data$p_adj < dea_p_cutoff,
        .data$log2fc < dea_log2fc_cutoff[[1]] |
          .data$log2fc > dea_log2fc_cutoff[[2]]
      ) |>
      dplyr::summarise(
        proteins = list(unique(.data$protein)),
        .by = tidyselect::all_of(by)
      ) |>
      tibble::deframe()
  }
  .gc_ora_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    ...,
    pro_list_fun = pro_list_fun
  )
}

#' The implementation template of gc_ora functions
#'
#' The only difference between different `.gc_ora` methods is how to extract protein lists.
#' Other operations, like argument validation, performing enrichment,
#' and packaging the result list, are exactly the same.
#' This function uses a `pro_list_fun` parameter to enable caller functions
#' provide their custom protein list extraction logic.
#'
#' @inheritParams .gc_ora
#' @param pro_list_fun A function with signature `function(dea_res)` that returns a
#'   named list where names are trait names and values are character vectors of protein Uniprot IDs.
#' @noRd
.gc_ora_impl <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ...,
  pro_list_fun = NULL
) {
  # Argument validation
  .check_dea_res(dea_res)
  .check_p_cutoff_arg(dea_p_cutoff)
  .check_log2fc_cutoff_arg(dea_log2fc_cutoff)

  # Performing enrichment
  protein_list <- pro_list_fun(dea_res)

  n_traits <- length(names(protein_list))
  cli::cli_alert_info(
    "Enriching for {.val {n_traits}} glycan traits... (This can take long)"
  )

  suppressWarnings(
    ck <- rlang::exec(
      clusterProfiler::compareCluster,
      protein_list,
      fun = enrich_fun,
      ...
    )
  )

  if (is.null(ck)) {
    cli::cli_alert_warning("No terms were enriched. `NULL` will be returned.")
    return(NULL)
  }
  tidy_res <- tibble::as_tibble(ck) |>
    janitor::clean_names() |>
    dplyr::rename(tidyselect::all_of(c(
      "trait" = "cluster",
      "p_val" = "pvalue",
      "p_adj" = "p_adjust",
      "q_val" = "qvalue"
    )))
  res <- list(tidy_result = tidy_res, raw_result = ck)
  structure(res, class = c(result_class, "glyfun_gc_ora_res", "glyfun_res"))
}
