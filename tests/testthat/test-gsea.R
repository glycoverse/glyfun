skip_if_not_installed("clusterProfiler")

.mock_gsea_dea_df <- function() {
  tibble::tibble(
    protein = c("P01308", "P01308", "P04637", "P42345"),
    site = c("s1", "s2", "s1", "s1"),
    trait = c("t1", "t1", "t1", "t1"),
    p_val = c(0.01, 0.04, 0.001, 0.2),
    log2fc = c(2, -1, -3, 0.5)
  )
}

test_that("enrich_gsea_go forwards args to gseGO through all layers", {
  captured <- NULL
  mock_gse_go <- function(
    geneList,
    OrgDb = NULL,
    keyType = NULL,
    ont = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    exponent = NULL,
    eps = NULL,
    seed = NULL,
    ...
  ) {
    captured <<- list(
      geneList = geneList,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    tibble::tibble(ID = "GO:0008150")
  }

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(gseGO = mock_gse_go, .package = "clusterProfiler")

  result <- suppressMessages(
    suppressWarnings(
      enrich_gsea_go(
        .mock_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        orgdb = "org.Hs.eg.db",
        ont = "BP",
        p_adj_method = "BY",
        p_cutoff = 0.01,
        min_gs_size = 3,
        max_gs_size = 99,
        exponent = 2,
        eps = 1e-8,
        seed = FALSE
      )
    )
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(captured$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$keyType, "UNIPROT")
  expect_identical(captured$ont, "BP")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(unname(captured$geneList), c(2, 0.5, -3))
  expect_equal(names(captured$geneList), c("P01308", "P42345", "P04637"))
})

test_that("enrich_gsea_reactome forwards args to gsePathway through all layers", {
  skip_if_not_installed("ReactomePA")

  captured <- NULL
  mock_gse_pathway <- function(
    geneList,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    exponent = NULL,
    eps = NULL,
    seed = NULL,
    ...
  ) {
    captured <<- list(
      geneList = geneList,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    tibble::tibble(ID = "R-HSA-123456")
  }

  local_mocked_bindings(
    gsePathway = mock_gse_pathway,
    .package = "ReactomePA"
  )
  local_mocked_bindings(
    .reactome_orgdb = function(organism) paste0("MOCK_REACTOME_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = TRUE) {
      expect_identical(drop_na, FALSE)
      expect_identical(orgdb, "MOCK_REACTOME_human")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gsea_reactome(
        .mock_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        organism = "human",
        p_adj_method = "BY",
        p_cutoff = 0.01,
        min_gs_size = 3,
        max_gs_size = 99,
        exponent = 2,
        eps = 1e-8,
        seed = FALSE
      )
    )
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(captured$organism, "human")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(unname(captured$geneList), c(2, 0.5, -3))
  expect_equal(
    names(captured$geneList),
    c("ENTREZ_P01308", "ENTREZ_P42345", "ENTREZ_P04637")
  )
})

test_that("enrich_gsea_wp forwards args to gseWP through all layers", {
  captured <- NULL
  mock_gse_wp <- function(
    geneList,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    exponent = NULL,
    eps = NULL,
    seed = NULL,
    ...
  ) {
    captured <<- list(
      geneList = geneList,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    tibble::tibble(ID = "WP123")
  }

  local_mocked_bindings(gseWP = mock_gse_wp, .package = "clusterProfiler")
  local_mocked_bindings(
    .wp_orgdb = function(organism) paste0("MOCK_WP_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = TRUE) {
      expect_identical(drop_na, FALSE)
      expect_identical(orgdb, "MOCK_WP_Homo sapiens")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gsea_wp(
        .mock_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        organism = "Homo sapiens",
        p_adj_method = "BY",
        p_cutoff = 0.01,
        min_gs_size = 3,
        max_gs_size = 99,
        exponent = 2,
        eps = 1e-8,
        seed = FALSE
      )
    )
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(captured$organism, "Homo sapiens")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(unname(captured$geneList), c(2, 0.5, -3))
  expect_equal(
    names(captured$geneList),
    c("ENTREZ_P01308", "ENTREZ_P42345", "ENTREZ_P04637")
  )
})

test_that("enrich_gsea_do forwards args to gseDO through all layers", {
  skip_if_not_installed("DOSE")

  captured <- NULL
  mock_gse_do <- function(
    geneList,
    ont = NULL,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    exponent = NULL,
    eps = NULL,
    seed = NULL,
    ...
  ) {
    captured <<- list(
      geneList = geneList,
      ont = ont,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    tibble::tibble(ID = "DOID:1234")
  }

  local_mocked_bindings(gseDO = mock_gse_do, .package = "DOSE")
  local_mocked_bindings(
    .do_orgdb = function(organism) paste0("MOCK_DO_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = TRUE) {
      expect_identical(drop_na, FALSE)
      expect_identical(orgdb, "MOCK_DO_hsa")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gsea_do(
        .mock_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        ont = "HDO",
        organism = "hsa",
        p_adj_method = "BY",
        p_cutoff = 0.01,
        min_gs_size = 3,
        max_gs_size = 99,
        exponent = 2,
        eps = 1e-8,
        seed = FALSE
      )
    )
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(captured$ont, "HDO")
  expect_identical(captured$organism, "hsa")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(unname(captured$geneList), c(2, 0.5, -3))
  expect_equal(
    names(captured$geneList),
    c("ENTREZ_P01308", "ENTREZ_P42345", "ENTREZ_P04637")
  )
})

test_that("enrich_gsea_ncg forwards args to gseNCG through all layers", {
  skip_if_not_installed("DOSE")

  captured <- NULL
  mock_gse_ncg <- function(
    geneList,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    exponent = NULL,
    eps = NULL,
    seed = NULL,
    ...
  ) {
    captured <<- list(
      geneList = geneList,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    tibble::tibble(ID = "NCG1")
  }

  local_mocked_bindings(gseNCG = mock_gse_ncg, .package = "DOSE")
  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = TRUE) {
      expect_identical(drop_na, FALSE)
      expect_identical(orgdb, "MOCK_org.Hs.eg.db")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gsea_ncg(
        .mock_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        p_adj_method = "BY",
        p_cutoff = 0.01,
        min_gs_size = 3,
        max_gs_size = 99,
        exponent = 2,
        eps = 1e-8,
        seed = FALSE
      )
    )
  )

  expect_s3_class(result, "tbl_df")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(unname(captured$geneList), c(2, 0.5, -3))
  expect_equal(
    names(captured$geneList),
    c("ENTREZ_P01308", "ENTREZ_P42345", "ENTREZ_P04637")
  )
})

test_that(".prepare_pro_list uses p_adj for p-based rank metrics when available", {
  df <- tibble::tibble(
    protein = c("P1", "P1", "P2"),
    site = c("s1", "s2", "s1"),
    trait = c("t1", "t1", "t1"),
    p_val = c(0.2, 1e-6, 0.001),
    p_adj = c(0.01, 0.5, 0.001),
    log2fc = c(2, -1, -3)
  )

  res <- glyfun:::.prepare_pro_list(
    df,
    rank_by = "signed_log10p",
    aggr = "median"
  )

  expected <- c(P1 = median(c(2, -0.30103)), P2 = -3)
  expect_equal(res[["P1"]], expected[["P1"]], tolerance = 1e-6)
  expect_equal(res[["P2"]], expected[["P2"]], tolerance = 1e-6)
  expect_equal(names(res), c("P1", "P2"))
})

test_that(".calcu_rank_scores supports all rank_by options", {
  df <- tibble::tibble(
    protein = c("P1", "P2"),
    site = c("s1", "s2"),
    trait = c("t1", "t1"),
    p_val = c(0.01, 0.1),
    log2fc = c(2, -1)
  )

  expect_equal(glyfun:::.calcu_rank_scores(df, "log2fc"), c(2, -1))
  expect_equal(glyfun:::.calcu_rank_scores(df, "abs_log2fc"), c(2, 1))
  expect_equal(
    glyfun:::.calcu_rank_scores(df, "log10p"),
    c(-log10(0.01), -log10(0.1))
  )
  expect_equal(
    glyfun:::.calcu_rank_scores(df, "signed_log10p"),
    c(2, -1)
  )
  expect_equal(
    glyfun:::.calcu_rank_scores(df, "log2fc_log10p"),
    c(4, -1)
  )
  expect_error(
    glyfun:::.calcu_rank_scores(df, "not_a_rank_metric"),
    "Invalid rank_by method"
  )
})

test_that(".gsea_aggr_fun returns aggregation function and errors on invalid input", {
  x <- c(1, 2, 10)

  expect_identical(glyfun:::.gsea_aggr_fun("median")(x), stats::median(x))
  expect_identical(glyfun:::.gsea_aggr_fun("mean")(x), mean(x))
  expect_identical(glyfun:::.gsea_aggr_fun("max")(x), max(x))
  expect_error(glyfun:::.gsea_aggr_fun("invalid"), "Invalid aggr method")
})

test_that(".gsea_impl converts UniProt names when uniprot_to_entrez is TRUE", {
  captured <- NULL
  mock_enrich <- function(geneList, ...) {
    captured <<- list(geneList = geneList, dots = list(...))
    tibble::tibble(ID = "GO:0003674")
  }
  local_mocked_bindings(
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = FALSE) {
      expect_identical(orgdb, "MOCK_ORGDB")
      expect_identical(drop_na, FALSE)
      c("ENTREZ_1", NA_character_)[seq_along(uniprot)]
    },
    .package = "glyfun"
  )

  dea_res <- .mock_gsea_dea_df()
  result <- suppressMessages(
    suppressWarnings(
      glyfun:::.gsea_impl(
        dea_res = dea_res,
        enrich_fun = mock_enrich,
        result_class = "unused",
        bitr_orgdb = "MOCK_ORGDB",
        uniprot_to_entrez = TRUE,
        pro_fun = function(dea_res) {
          expect_equal(nrow(dea_res), 4)
          stats::setNames(c(2, -1), c("P01308", "P04637"))
        }
      )
    )
  )

  expect_equal(unname(captured$geneList), 2)
  expect_equal(names(captured$geneList), "ENTREZ_1")
  expect_s3_class(result, "tbl_df")
})

test_that(".gsea_impl collapses duplicate Entrez IDs after conversion", {
  captured <- NULL
  mock_enrich <- function(geneList, ...) {
    captured <<- list(geneList = geneList, dots = list(...))
    tibble::tibble(ID = "R-HSA-123456")
  }
  local_mocked_bindings(
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = FALSE) {
      expect_identical(orgdb, "MOCK_ORGDB")
      expect_identical(drop_na, FALSE)
      dplyr::recode(
        uniprot,
        P01308 = "ENTREZ_DUP",
        P42345 = "ENTREZ_UNIQUE",
        P04637 = "ENTREZ_DUP"
      )
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      glyfun:::.gsea_impl(
        dea_res = .mock_gsea_dea_df(),
        enrich_fun = mock_enrich,
        result_class = "unused",
        bitr_orgdb = "MOCK_ORGDB",
        aggr = "max",
        uniprot_to_entrez = TRUE,
        pro_fun = function(dea_res) {
          expect_equal(nrow(dea_res), 4)
          stats::setNames(c(2, 0.5, -3), c("P01308", "P42345", "P04637"))
        }
      )
    )
  )

  expect_equal(unname(captured$geneList), c(2, 0.5))
  expect_equal(names(captured$geneList), c("ENTREZ_DUP", "ENTREZ_UNIQUE"))
  expect_s3_class(result, "tbl_df")
})

test_that(".gsea_impl errors when uniprot_to_entrez is TRUE but bitr_orgdb is missing", {
  expect_error(
    glyfun:::.gsea_impl(
      dea_res = .mock_gsea_dea_df(),
      enrich_fun = function(...) tibble::tibble(ID = "GO:1"),
      result_class = "unused",
      uniprot_to_entrez = TRUE,
      pro_fun = function(dea_res) {
        stats::setNames(dea_res$log2fc, dea_res$protein)
      }
    ),
    "bitr_orgdb"
  )
})

test_that(".gsea_impl returns NULL when no terms are enriched", {
  result <- suppressMessages(
    suppressWarnings(
      glyfun:::.gsea_impl(
        dea_res = .mock_gsea_dea_df(),
        enrich_fun = function(...) tibble::tibble(),
        result_class = "unused",
        pro_fun = function(dea_res) {
          stats::setNames(dea_res$log2fc, dea_res$protein)
        }
      )
    )
  )
  expect_null(result)
})

test_that(".uniprot_to_entrez preserves input length and order for duplicated UniProt IDs", {
  local_mocked_bindings(
    bitr = function(geneID, fromType, toType, OrgDb) {
      expect_equal(geneID, c("P1", "P2"))
      expect_identical(fromType, "UNIPROT")
      expect_identical(toType, "ENTREZID")
      expect_identical(OrgDb, "MOCK_ORGDB")
      tibble::tibble(
        UNIPROT = c("P1", "P2"),
        ENTREZID = c("E1", "E2")
      )
    },
    .package = "clusterProfiler"
  )

  res <- suppressMessages(
    suppressWarnings(
      glyfun:::.uniprot_to_entrez(
        c("P1", "P1", "P2"),
        orgdb = "MOCK_ORGDB",
        drop_na = FALSE
      )
    )
  )

  expect_equal(res, c("E1", "E1", "E2"))
})
