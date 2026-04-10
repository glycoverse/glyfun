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

test_that(".gsea.data.frame prepares a ranked protein list and calls .gsea_impl", {
  captured <- NULL
  sentinel <- tibble::tibble(ID = "GO:0003674")
  mock_gsea_impl <- function(
    dea_res,
    enrich_fun,
    result_class,
    bitr_orgdb = NULL,
    ...,
    pro_fun = NULL,
    uniprot_to_entrez = FALSE
  ) {
    captured <<- list(
      proteins = pro_fun(dea_res),
      enrich_fun = enrich_fun,
      result_class = result_class,
      bitr_orgdb = bitr_orgdb,
      dots = list(...),
      uniprot_to_entrez = uniprot_to_entrez
    )
    sentinel
  }
  local_mocked_bindings(.gsea_impl = mock_gsea_impl, .package = "glyfun")

  result <- glyfun:::.gsea(
    dea_res = .mock_gsea_dea_df(),
    enrich_fun = function(...) NULL,
    result_class = "mock_result_class",
    rank_by = "log2fc",
    aggr = "max",
    test_arg = "forwarded"
  )

  expect_identical(result, sentinel)
  expect_equal(unname(captured$proteins), c(2, 0.5, -3))
  expect_equal(names(captured$proteins), c("P01308", "P42345", "P04637"))
  expect_identical(captured$result_class, "mock_result_class")
  expect_identical(captured$dots$test_arg, "forwarded")
  expect_false(captured$uniprot_to_entrez)
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

  res <- glyfun:::.prepare_pro_list(df, rank_by = "signed_log10p", aggr = "median")

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
