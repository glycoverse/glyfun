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
#' @return A clusterProfiler `compareClusterResult` object with additional `glyfun` classes.
#'   It can be readily converted to a tibble with [tibble::as_tibble()],
#'   or visualized with `clusterProfiler` functions like [clusterProfiler::dotplot()].
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
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  orgdb <- .prepare_orgdb(orgdb)
  .gc_ora(
    dea_res,
    enrich_fun = clusterProfiler::enrichGO,
    result_class = "glyfun_gc_ora_go_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    keyType = "UNIPROT",
    OrgDb = orgdb,
    ont = ont,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
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
#' @inherit enrich_gc_ora_go return
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
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  .gc_ora(
    dea_res,
    enrich_fun = clusterProfiler::enrichKEGG,
    result_class = "glyfun_gc_ora_kegg_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    keyType = "uniprot",
    organism = organism,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )
}

#' Glycan-Centric Reactome Pathway Over Representation Analysis
#'
#' @description
#' Performs glycan-centric Reactome pathway Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits to biological pathways. It helps answer questions like "Which
#' Reactome pathways are enriched in proteins with a specific dysregulated glycan motif?",
#' by grouping differential analysis results by glycan traits and computing
#' pathway enrichment for each trait.
#'
#' @inheritSection enrich_gc_ora_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_reactome
#' @inherit enrich_gc_ora_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [ReactomePA::enrichPathway()]
#' @export
enrich_gc_ora_reactome <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  organism = "human",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  rlang::check_installed("ReactomePA")
  orgdb <- .reactome_orgdb(organism)
  .gc_ora(
    dea_res,
    enrich_fun = ReactomePA::enrichPathway,
    result_class = "glyfun_gc_ora_reactome_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    bitr_orgdb = orgdb,
    organism = organism,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric WikiPathways Over Representation Analysis
#'
#' @description
#' Performs glycan-centric WikiPathways Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits to biological pathways. It helps answer questions like "Which
#' WikiPathways are enriched in proteins with a specific dysregulated glycan motif?",
#' by grouping differential analysis results by glycan traits and computing
#' pathway enrichment for each trait.
#'
#' @inheritSection enrich_gc_ora_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_wp
#' @inherit enrich_gc_ora_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::enrichWP()]
#' @export
enrich_gc_ora_wp <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  organism = "Homo sapiens",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  orgdb <- .wp_orgdb(organism)
  .gc_ora(
    dea_res,
    enrich_fun = clusterProfiler::enrichWP,
    result_class = "glyfun_gc_ora_wp_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    bitr_orgdb = orgdb,
    organism = organism,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric Disease Ontology (DO) Over Representation Analysis
#'
#' @description
#' Performs glycan-centric Disease Ontology (DO) Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits to disease associations. It helps answer questions like "Which
#' diseases are enriched in proteins with a specific dysregulated glycan motif?",
#' by grouping differential analysis results by glycan traits and computing
#' disease enrichment for each trait.
#'
#' @inheritSection enrich_gc_ora_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_do
#' @inherit enrich_gc_ora_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [DOSE::enrichDO()]
#' @export
enrich_gc_ora_do <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  ont = "HDO",
  organism = "hsa",
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  rlang::check_installed("DOSE")
  orgdb <- .do_orgdb(organism)
  .gc_ora(
    dea_res,
    enrich_fun = DOSE::enrichDO,
    result_class = "glyfun_gc_ora_do_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    bitr_orgdb = orgdb,
    ont = ont,
    organism = organism,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric Network of Cancer Genes (NCG) Over Representation Analysis
#'
#' @description
#' Performs glycan-centric Network of Cancer Genes (NCG) Over-Representation Analysis (ORA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits to cancer gene associations. It helps answer questions like "Which
#' cancer gene sets are enriched in proteins with a specific dysregulated glycan motif?",
#' by grouping differential analysis results by glycan traits and computing
#' cancer gene enrichment for each trait.
#'
#' @inheritSection enrich_gc_ora_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_ora_go Common usage pattern
#'
#' @inheritParams enrich_ora_ncg
#' @inherit enrich_gc_ora_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [DOSE::enrichNCG()]
#' @export
enrich_gc_ora_ncg <- function(
  dea_res,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500
) {
  rlang::check_installed("DOSE")
  orgdb <- .prepare_orgdb("org.Hs.eg.db")
  .gc_ora(
    dea_res,
    enrich_fun = DOSE::enrichNCG,
    result_class = "glyfun_gc_ora_ncg_res",
    dea_p_cutoff = dea_p_cutoff,
    dea_log2fc_cutoff = dea_log2fc_cutoff,
    bitr_orgdb = orgdb,
    universe = universe,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    uniprot_to_entrez = TRUE
  )
}

#' Perform Glycan-Centric ORA
#' @param dea_res DEA result from glystats or a tibble.
#' @param enrich_fun An enrichment function.
#' @param result_class A string of the concrete result class.
#' @param dea_p_cutoff P-value cutoff to define statistical significance.
#' @param dea_log2fc_cutoff Log2FC cutoffs to define statistical significance.
#' @param bitr_orgdb OrgDb object for bitr conversion.
#' @param ... Parameters passed to downstream enrichment package.
#' @noRd
.gc_ora <- function(
  dea_res,
  enrich_fun,
  result_class,
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  bitr_orgdb = NULL,
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
  bitr_orgdb = NULL,
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
    bitr_orgdb = bitr_orgdb,
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
  bitr_orgdb = NULL,
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
    bitr_orgdb = bitr_orgdb,
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
  universe = NULL,
  bitr_orgdb = NULL,
  ...,
  pro_list_fun = NULL,
  uniprot_to_entrez = FALSE
) {
  # Argument validation
  if (uniprot_to_entrez && is.null(bitr_orgdb)) {
    cli::cli_abort(
      "{.arg bitr_orgdb} must be provided when {.arg uniprot_to_entrez} is TRUE."
    )
  }
  .check_dea_res(dea_res)
  .check_p_cutoff_arg(dea_p_cutoff)
  .check_log2fc_cutoff_arg(dea_log2fc_cutoff)

  # Performing enrichment
  protein_list <- pro_list_fun(dea_res)

  # Convert Uniprot IDs to Entrez IDs if needed (for ReactomePA)
  if (uniprot_to_entrez) {
    protein_list <- .uniprot_to_entrez_prolist(protein_list, bitr_orgdb)
    if (!is.null(universe)) {
      universe <- .uniprot_to_entrez(universe, bitr_orgdb)
    }
  }

  n_traits <- length(names(protein_list))
  cli::cli_alert_info(
    "Enriching for {.val {n_traits}} glycan traits... (This can take long)"
  )

  suppressWarnings(
    res <- suppressPackageStartupMessages(
      rlang::exec(
        clusterProfiler::compareCluster,
        protein_list,
        fun = enrich_fun,
        universe = universe,
        ...
      )
    )
  )

  if (is.null(res)) {
    cli::cli_alert_warning("No terms were enriched. `NULL` will be returned.")
  }
  res
}
