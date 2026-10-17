#' Unified functional enrichment for glycoproteomics
#'
#' @description
#' Runs protein-centric or glycan-trait-centric over-representation analysis
#' (ORA) or gene set enrichment analysis (GSEA) through one entry point.
#'
#' Protein-centric enrichment analyzes proteins with dysregulated glycosylation.
#' Glycan-trait-centric enrichment analyzes protein sets or protein rankings
#' stratified by a glycan trait; it does not directly enrich glycan structures.
#'
#' @param dea_res Differential analysis result. It can be a supported
#'   `glystats` result or a data frame accepted by the selected existing
#'   `enrich_*()` family.
#' @param method Enrichment method: `"ora"` or `"gsea"`.
#' @param database Annotation database: `"GO"`, `"KEGG"`, `"Reactome"`,
#'   `"WikiPathways"`, `"DO"`, or `"NCG"`.
#' @param scope Analysis scope: `"protein"` for protein-centric enrichment or
#'   `"glycan_trait"` for glycan-trait-centric enrichment.
#' @param dea_p_cutoff P-value cutoff used to select proteins for ORA.
#' @param dea_log2fc_cutoff Log2 fold-change cutoffs used to select proteins
#'   for ORA. The value must contain a negative and a positive boundary.
#' @param rank_by Ranking criterion used for GSEA. See [enrich_gsea_go()].
#' @param aggr Aggregation method used to combine scores for the same protein
#'   in GSEA. See [enrich_gsea_go()].
#' @param orgdb An `OrgDb` object or its package name. It defaults to
#'   `"org.Hs.eg.db"` for GO and NCG.
#' @param ont Ontology used by GO or DO. It defaults to `"MF"` for GO and
#'   `"HDO"` for DO.
#' @param organism Organism passed to the selected database. Defaults are
#'   `"hsa"` for KEGG and DO, `"human"` for Reactome, and `"Homo sapiens"`
#'   for WikiPathways.
#' @param universe Background proteins for ORA.
#' @param p_adj_method P-value adjustment method passed to the enrichment
#'   backend.
#' @param p_cutoff P-value cutoff passed to the enrichment backend.
#' @param q_cutoff Q-value cutoff passed to the enrichment backend for ORA.
#' @param min_gs_size Minimum gene-set size.
#' @param max_gs_size Maximum gene-set size.
#' @param exponent GSEA weighting exponent.
#' @param eps Epsilon used for GSEA p-value calculation.
#' @param seed Whether to set a random seed in GSEA.
#' @param by Glycan-trait grouping column when `scope = "glycan_trait"`.
#'   See [enrich_gc_ora_go()] for accepted values and automatic inference.
#'   Must be `NULL` for `scope = "protein"`.
#' @param ... Additional arguments passed to the selected enrichment backend.
#'
#' @return A clusterProfiler enrichment result. ORA returns an `enrichResult`
#'   or `compareClusterResult`, and GSEA returns a `gseaResult` or
#'   `compareClusterResult`.
#'
#' @examples
#' \dontrun{
#' enrich_glyco(dea_res, method = "ora", database = "GO")
#' enrich_glyco(
#'   dea_res,
#'   method = "gsea",
#'   database = "Reactome",
#'   scope = "glycan_trait"
#' )
#' }
#'
#' @seealso [enrich_ora_go()], [enrich_gsea_go()], [enrich_gc_ora_go()],
#'   [enrich_gc_gsea_go()]
#' @export
enrich_glyco <- function(
  dea_res,
  method = c("ora", "gsea"),
  database = c("GO", "KEGG", "Reactome", "WikiPathways", "DO", "NCG"),
  scope = c("protein", "glycan_trait"),
  dea_p_cutoff = 0.05,
  dea_log2fc_cutoff = c(-1, 1),
  rank_by = "signed_log10p",
  aggr = "median",
  orgdb = NULL,
  ont = NULL,
  organism = NULL,
  universe = NULL,
  p_adj_method = "BH",
  p_cutoff = 0.05,
  q_cutoff = 0.2,
  min_gs_size = 10,
  max_gs_size = 500,
  exponent = 1,
  eps = 1e-10,
  seed = FALSE,
  ...,
  by = NULL
) {
  method <- rlang::arg_match(method)
  database <- rlang::arg_match(database)
  scope <- rlang::arg_match(scope)
  if (scope == "protein" && !is.null(by)) {
    cli::cli_abort(
      "{.arg by} is only supported when {.code scope = \"glycan_trait\"}."
    )
  }

  backend <- .enrich_glyco_backend(
    method = method,
    database = database,
    orgdb = orgdb,
    ont = ont,
    organism = organism
  )

  helper <- switch(
    paste(method, scope, sep = "."),
    ora.protein = .ora,
    gsea.protein = .gsea,
    ora.glycan_trait = .gc_ora,
    gsea.glycan_trait = .gc_gsea
  )

  helper_args <- list(
    dea_res = dea_res,
    enrich_fun = backend$enrich_fun,
    bitr_orgdb = backend$bitr_orgdb,
    uniprot_to_entrez = backend$uniprot_to_entrez
  )
  if (scope == "glycan_trait") {
    helper_args$by <- by
  }

  if (method == "ora") {
    helper_args <- c(
      helper_args,
      list(
        dea_p_cutoff = dea_p_cutoff,
        dea_log2fc_cutoff = dea_log2fc_cutoff,
        universe = universe,
        pAdjustMethod = p_adj_method,
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )
    )
  } else {
    helper_args <- c(
      helper_args,
      list(
        rank_by = rank_by,
        aggr = aggr,
        pAdjustMethod = p_adj_method,
        pvalueCutoff = p_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        exponent = exponent,
        eps = eps,
        seed = seed
      )
    )
  }

  do.call(
    helper,
    c(helper_args, backend$backend_args, list(...))
  )
}

.enrich_glyco_backend <- function(method, database, orgdb, ont, organism) {
  if (is.null(orgdb)) {
    orgdb <- "org.Hs.eg.db"
  }

  if (is.null(ont)) {
    ont <- switch(database, GO = "MF", DO = "HDO", NULL)
  }

  if (is.null(organism)) {
    organism <- switch(
      database,
      KEGG = "hsa",
      Reactome = "human",
      WikiPathways = "Homo sapiens",
      DO = "hsa",
      NULL
    )
  }

  switch(
    database,
    GO = {
      orgdb <- .prepare_orgdb(orgdb)
      list(
        enrich_fun = switch(
          method,
          ora = clusterProfiler::enrichGO,
          gsea = clusterProfiler::gseGO
        ),
        bitr_orgdb = NULL,
        uniprot_to_entrez = FALSE,
        backend_args = list(
          keyType = "UNIPROT",
          OrgDb = orgdb,
          ont = ont
        )
      )
    },
    KEGG = list(
      enrich_fun = switch(
        method,
        ora = clusterProfiler::enrichKEGG,
        gsea = clusterProfiler::gseKEGG
      ),
      bitr_orgdb = NULL,
      uniprot_to_entrez = FALSE,
      backend_args = list(keyType = "uniprot", organism = organism)
    ),
    Reactome = {
      rlang::check_installed("ReactomePA")
      list(
        enrich_fun = switch(
          method,
          ora = ReactomePA::enrichPathway,
          gsea = ReactomePA::gsePathway
        ),
        bitr_orgdb = .reactome_orgdb(organism),
        uniprot_to_entrez = TRUE,
        backend_args = list(organism = organism)
      )
    },
    WikiPathways = {
      list(
        enrich_fun = switch(
          method,
          ora = clusterProfiler::enrichWP,
          gsea = clusterProfiler::gseWP
        ),
        bitr_orgdb = .wp_orgdb(organism),
        uniprot_to_entrez = TRUE,
        backend_args = list(organism = organism)
      )
    },
    DO = {
      rlang::check_installed("DOSE")
      list(
        enrich_fun = switch(
          method,
          ora = DOSE::enrichDO,
          gsea = DOSE::gseDO
        ),
        bitr_orgdb = .do_orgdb(organism),
        uniprot_to_entrez = TRUE,
        backend_args = list(ont = ont, organism = organism)
      )
    },
    NCG = {
      rlang::check_installed("DOSE")
      orgdb <- .prepare_orgdb(orgdb)
      list(
        enrich_fun = switch(
          method,
          ora = DOSE::enrichNCG,
          gsea = DOSE::gseNCG
        ),
        bitr_orgdb = orgdb,
        uniprot_to_entrez = TRUE,
        backend_args = list()
      )
    }
  )
}
