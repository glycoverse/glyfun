#' Glycan-Centric GO Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA).
#' Instead of traditional protein-centric enrichment, this function links specific
#' glycan traits (e.g., core-fucosylation, sialylation) to functional annotations.
#' It ranks proteins within each glycan trait separately, then compares enriched
#' biological functions across traits by running GSEA for each trait.
#'
#' @details
#' # What is glycan-centric enrichment?
#'
#' In traditional glycoproteomics data analysis,
#' we usually perform differential expression analysis (DEA) on glycoforms,
#' extract proteins that have dysregulated glycosylation,
#' then perform functional enrichment (e.g. GO) on these proteins.
#' This is what `enrich_xxx()` functions do (e.g. [enrich_gsea_go()]).
#'
#' `enrich_gc_xxx()` functions differ in that they link specific glycan traits with functional annotations.
#' Instead of answering the question
#' "Which functions are enriched in dysregulated glycoproteins?",
#' `enrich_gc_xxx()` answers questions like
#' "Which functions are enriched in proteins ranked highly for core-fucosylation changes?"
#' Higher specificity, deeper insights. By focusing on distinct glycan motifs,
#' it helps you pinpoint the functional relevance of specific glycosylation changes.
#'
#' # How it ranks proteins
#'
#' GSEA requires a ranked list of proteins as input.
#' This function first splits the DEA result by glycan trait, then ranks proteins
#' within each trait separately. For each trait, it applies the same ranking and
#' aggregation logic as [enrich_gsea_go()], producing one ranked protein list per trait.
#' Those trait-specific ranked lists are then compared with
#' [clusterProfiler::compareCluster()] using its formula interface.
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
#' go_res <- enrich_gc_gsea_go(dea_res)  # or other `enrich_gc_gsea_xxx()` functions
#' ```
#'
#' @inheritParams enrich_gsea_go
#' @return A clusterProfiler `compareClusterResult` object.
#'   It can be readily converted to a tibble with [tibble::as_tibble()],
#'   or visualized with `clusterProfiler` functions like [clusterProfiler::dotplot()].
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::gseGO()]
#' @export
enrich_gc_gsea_go <- function(
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
  seed = FALSE
) {
  orgdb <- .prepare_orgdb(orgdb)
  .gc_gsea(
    dea_res,
    enrich_fun = clusterProfiler::gseGO,
    result_class = "glyfun_gc_gsea_go_res",
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

#' Glycan-Centric KEGG Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric KEGG pathway Gene Set Enrichment Analysis (GSEA).
#' It ranks proteins within each glycan trait separately, then compares enriched
#' pathways across traits by running GSEA for each trait.
#'
#' @inheritSection enrich_gc_gsea_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_gsea_go How it ranks proteins
#' @inheritSection enrich_gc_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_kegg
#' @inherit enrich_gc_gsea_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::gseKEGG()]
#' @export
enrich_gc_gsea_kegg <- function(
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
  seed = FALSE
) {
  .gc_gsea(
    dea_res,
    enrich_fun = clusterProfiler::gseKEGG,
    result_class = "glyfun_gc_gsea_kegg_res",
    rank_by = rank_by,
    aggr = aggr,
    organism = organism,
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

#' Glycan-Centric Reactome Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric Reactome pathway Gene Set Enrichment Analysis (GSEA).
#' It ranks proteins within each glycan trait separately, then compares enriched
#' pathways across traits by running GSEA for each trait.
#'
#' @inheritSection enrich_gc_gsea_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_gsea_go How it ranks proteins
#' @inheritSection enrich_gc_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_reactome
#' @inherit enrich_gc_gsea_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [ReactomePA::gsePathway()]
#' @export
enrich_gc_gsea_reactome <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  organism = "human",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE
) {
  rlang::check_installed("ReactomePA")
  orgdb <- .reactome_orgdb(organism)
  .gc_gsea(
    dea_res,
    enrich_fun = ReactomePA::gsePathway,
    result_class = "glyfun_gc_gsea_reactome_res",
    rank_by = rank_by,
    aggr = aggr,
    bitr_orgdb = orgdb,
    organism = organism,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric WikiPathways Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric WikiPathways Gene Set Enrichment Analysis (GSEA).
#' It ranks proteins within each glycan trait separately, then compares enriched
#' pathways across traits by running GSEA for each trait.
#'
#' @inheritSection enrich_gc_gsea_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_gsea_go How it ranks proteins
#' @inheritSection enrich_gc_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_wp
#' @inherit enrich_gc_gsea_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::gseWP()]
#' @export
enrich_gc_gsea_wp <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  organism = "Homo sapiens",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE
) {
  orgdb <- .wp_orgdb(organism)
  .gc_gsea(
    dea_res,
    enrich_fun = clusterProfiler::gseWP,
    result_class = "glyfun_gc_gsea_wp_res",
    rank_by = rank_by,
    aggr = aggr,
    bitr_orgdb = orgdb,
    organism = organism,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric Disease Ontology (DO) Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric Disease Ontology (DO) Gene Set Enrichment Analysis (GSEA).
#' It ranks proteins within each glycan trait separately, then compares enriched
#' disease terms across traits by running GSEA for each trait.
#'
#' @inheritSection enrich_gc_gsea_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_gsea_go How it ranks proteins
#' @inheritSection enrich_gc_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_do
#' @inherit enrich_gc_gsea_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [DOSE::gseDO()]
#' @export
enrich_gc_gsea_do <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  ont = "HDO",
  organism = "hsa",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE
) {
  rlang::check_installed("DOSE")
  orgdb <- .do_orgdb(organism)
  .gc_gsea(
    dea_res,
    enrich_fun = DOSE::gseDO,
    result_class = "glyfun_gc_gsea_do_res",
    rank_by = rank_by,
    aggr = aggr,
    bitr_orgdb = orgdb,
    ont = ont,
    organism = organism,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed,
    uniprot_to_entrez = TRUE
  )
}

#' Glycan-Centric Network of Cancer Genes (NCG) Gene Set Enrichment Analysis
#'
#' @description
#' Performs glycan-centric Network of Cancer Genes (NCG) Gene Set Enrichment Analysis (GSEA).
#' It ranks proteins within each glycan trait separately, then compares enriched
#' cancer gene sets across traits by running GSEA for each trait.
#'
#' @inheritSection enrich_gc_gsea_go What is glycan-centric enrichment?
#' @inheritSection enrich_gc_gsea_go How it ranks proteins
#' @inheritSection enrich_gc_gsea_go Common usage pattern
#'
#' @inheritParams enrich_gsea_ncg
#' @inherit enrich_gc_gsea_go return
#'
#' @seealso [clusterProfiler::compareCluster()], [DOSE::gseNCG()]
#' @export
enrich_gc_gsea_ncg <- function(
  dea_res,
  rank_by = "signed_log10p",
  aggr = "median",
  p_adj_method = "BH",
  p_cutoff = 0.05,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE
) {
  rlang::check_installed("DOSE")
  orgdb <- .prepare_orgdb("org.Hs.eg.db")
  .gc_gsea(
    dea_res,
    enrich_fun = DOSE::gseNCG,
    result_class = "glyfun_gc_gsea_ncg_res",
    rank_by = rank_by,
    aggr = aggr,
    bitr_orgdb = orgdb,
    pAdjustMethod = p_adj_method,
    pvalueCutoff = p_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    exponent = exponent,
    eps = eps,
    seed = seed,
    uniprot_to_entrez = TRUE
  )
}

#' Perform Glycan-Centric GSEA
#' @param dea_res DEA result from glystats or a tibble.
#' @param enrich_fun An enrichment function.
#' @param result_class A string of the concrete result class.
#' @param rank_by Criteria for ranking proteins.
#' @param aggr Aggregation method for combining multiple scores for the same protein.
#' @param bitr_orgdb OrgDb object for bitr conversion.
#' @param ... Parameters passed to downstream enrichment package.
#' @noRd
.gc_gsea <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  UseMethod(".gc_gsea")
}

.gc_gsea.data.frame <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  pro_df_fun <- function(dea_res) {
    .prepare_gc_gsea_df(
      dea_res,
      by = "trait",
      rank_by = rank_by,
      aggr = aggr
    )
  }
  .gc_gsea_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_df_fun = pro_df_fun
  )
}

.gc_gsea.glystats_res <- function(
  dea_res,
  enrich_fun,
  result_class,
  rank_by,
  aggr,
  bitr_orgdb = NULL,
  ...
) {
  pro_df_fun <- function(dea_res) {
    tidy_dea_res <- glystats::get_tidy_result(dea_res)
    by <- .infer_gc_trait_col(tidy_dea_res)
    .prepare_gc_gsea_df(
      tidy_dea_res,
      by = by,
      rank_by = rank_by,
      aggr = aggr
    )
  }
  .gc_gsea_impl(
    dea_res,
    enrich_fun = enrich_fun,
    result_class = result_class,
    bitr_orgdb = bitr_orgdb,
    ...,
    pro_df_fun = pro_df_fun
  )
}

#' Prepare trait-specific ranking data for glycan-centric GSEA
#' @param df A tidy DEA result.
#' @param by The column name representing glycan traits.
#' @param rank_by Criteria for ranking proteins.
#' @param aggr Aggregation method for combining multiple scores for the same protein.
#' @returns A tibble with `gene`, `score`, and `trait` columns.
#' @noRd
.prepare_gc_gsea_df <- function(df, by, rank_by, aggr) {
  trait_dfs <- df |>
    dplyr::summarise(
      data = list(dplyr::pick(dplyr::everything())),
      .by = tidyselect::all_of(by)
    ) |>
    dplyr::rename(trait = tidyselect::all_of(by))

  purrr::map2_dfr(trait_dfs$trait, trait_dfs$data, function(trait, trait_df) {
    pro_list <- .prepare_pro_list(trait_df, rank_by, aggr)
    tibble::tibble(
      gene = names(pro_list),
      score = unname(pro_list),
      trait = trait
    )
  })
}

#' The implementation template of gc_gsea functions
#'
#' The only difference between different `.gc_gsea` methods is how to extract
#' trait-specific ranked protein lists. Other operations, like argument validation,
#' ID conversion, performing enrichment, and packaging the result list, are exactly
#' the same. This function uses a `pro_df_fun` parameter to enable caller functions
#' provide their custom ranking data extraction logic.
#'
#' @inheritParams .gc_gsea
#' @param pro_df_fun A function with signature `function(dea_res)` that returns a
#'   tibble with `gene`, `score`, and `trait` columns.
#' @noRd
.gc_gsea_impl <- function(
  dea_res,
  enrich_fun,
  result_class,
  bitr_orgdb = NULL,
  ...,
  pro_df_fun = NULL,
  uniprot_to_entrez = FALSE
) {
  if (uniprot_to_entrez && is.null(bitr_orgdb)) {
    cli::cli_abort(
      "{.arg bitr_orgdb} must be provided when {.arg uniprot_to_entrez} is TRUE."
    )
  }
  .check_dea_res(dea_res)

  pro_df <- pro_df_fun(dea_res)
  if (uniprot_to_entrez) {
    pro_df <- pro_df |>
      dplyr::mutate(
        gene = .uniprot_to_entrez(.data$gene, bitr_orgdb, drop_na = FALSE)
      ) |>
      dplyr::filter(!is.na(.data$gene))
  }

  n_traits <- dplyr::n_distinct(pro_df$trait)
  cli::cli_alert_info(
    "Enriching for {.val {n_traits}} glycan traits... (This can take long)"
  )

  suppressWarnings(
    res <- suppressPackageStartupMessages(
      rlang::exec(
        clusterProfiler::compareCluster,
        stats::as.formula("gene | score ~ trait"),
        data = pro_df,
        fun = enrich_fun,
        ...
      )
    )
  )

  if (is.null(res)) {
    cli::cli_alert_warning("No terms were enriched. `NULL` will be returned.")
  }
  res
}
