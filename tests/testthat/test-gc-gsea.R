skip_if_not_installed("clusterProfiler")

.mock_gc_gsea_dea_df <- function() {
  tibble::tibble(
    protein = c("P01308", "P01308", "P04637", "P42345", "P00533", "P00533"),
    site = c("s1", "s2", "s1", "s1", "s1", "s2"),
    trait = c("trait_A", "trait_A", "trait_A", "trait_B", "trait_B", "trait_B"),
    p_val = c(0.01, 0.04, 0.001, 0.2, 0.03, 0.001),
    log2fc = c(2, -1, -3, 0.5, 1.5, 2.5)
  )
}

for (database in c("go", "kegg", "reactome", "wp", "do", "ncg")) {
  test_that(paste("explicit by controls glycan-centric GSEA for", database), {
    if (database == "reactome") {
      skip_if_not_installed("ReactomePA")
    }
    if (database %in% c("do", "ncg")) {
      skip_if_not_installed("DOSE")
    }

    captured <- NULL
    local_mocked_bindings(
      compareCluster = function(formula, data, ...) {
        captured <<- list(...)
        data
      },
      .package = "clusterProfiler"
    )
    local_mocked_bindings(
      .prepare_orgdb = function(...) "MOCK_ORGDB",
      .reactome_orgdb = function(...) "MOCK_ORGDB",
      .wp_orgdb = function(...) "MOCK_ORGDB",
      .do_orgdb = function(...) "MOCK_ORGDB",
      .uniprot_to_entrez = function(uniprot, ...) uniprot,
      .package = "glyfun"
    )
    df <- tibble::tibble(
      protein = c("P1", "P2", "P3", "P4"),
      site = "s1",
      trait = "T",
      glycan_structure = "S",
      glycan_composition = c("C1", "C2", "C1", "C2"),
      p_val = c(0.001, 0.001, 0.2, 0.001),
      p_adj = c(0.001, 0.2, 0.2, 0.001),
      log2fc = c(2, -2, 3, 0.5)
    )
    glystats_res <- structure(
      list(tidy_result = df, raw_result = list()),
      class = c("glystats_ttest_res", "glystats_res")
    )
    expected <- tibble::tibble(
      gene = c("P3", "P1", "P4", "P2"),
      score = c(3, 2, 0.5, -2),
      trait = c("C1", "C1", "C2", "C2")
    )
    enrich <- get(paste0("enrich_gc_gsea_", database))
    for (dea_res in list(df, glystats_res)) {
      expect_equal(
        suppressMessages(enrich(
          dea_res,
          by = "glycan_composition",
          rank_by = "log2fc"
        )),
        expected
      )
      expect_null(captured$by)
    }
  })
}

test_that("glycan-centric GSEA preserves automatic defaults for both input types", {
  local_mocked_bindings(
    compareCluster = function(formula, data, ...) data,
    .package = "clusterProfiler"
  )
  df <- .mock_gc_gsea_dea_df()
  df$glycan_structure <- "S"
  glystats_res <- structure(
    list(tidy_result = dplyr::mutate(df, p_adj = p_val), raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )
  auto_df <- suppressMessages(enrich_gc_gsea_kegg(df, rank_by = "log2fc"))
  auto_glystats <- suppressMessages(enrich_gc_gsea_kegg(
    glystats_res,
    rank_by = "log2fc"
  ))
  expect_equal(auto_df$trait, c("trait_A", "trait_A", "trait_B", "trait_B"))
  expect_equal(auto_glystats$trait, rep("S", 4))
  expect_equal(
    auto_df,
    suppressMessages(enrich_gc_gsea_kegg(df, by = "trait", rank_by = "log2fc"))
  )
  expect_equal(
    auto_glystats,
    suppressMessages(enrich_gc_gsea_kegg(
      glystats_res,
      by = "glycan_structure",
      rank_by = "log2fc"
    ))
  )
})

test_that("glycan-centric GSEA accepts custom grouping columns without trait", {
  local_mocked_bindings(
    compareCluster = function(formula, data, ...) data,
    .package = "clusterProfiler"
  )
  df <- .mock_gc_gsea_dea_df()
  expected <- suppressMessages(enrich_gc_gsea_kegg(df, rank_by = "log2fc"))
  df$data <- df$trait
  df$trait <- NULL
  expect_equal(
    suppressMessages(enrich_gc_gsea_kegg(df, by = "data", rank_by = "log2fc")),
    expected
  )
  by_protein <- suppressMessages(enrich_gc_gsea_kegg(
    df,
    by = "protein",
    rank_by = "log2fc"
  ))
  expect_identical(by_protein$gene, by_protein$trait)
  df$motif <- df$data
  expect_equal(
    suppressMessages(enrich_gc_gsea_kegg(df, rank_by = "log2fc")),
    expected
  )
})

test_that("glycan-centric GSEA rejects invalid by before enrichment", {
  df <- .mock_gc_gsea_dea_df()
  expect_snapshot(error = TRUE, enrich_gc_gsea_kegg(df, by = "missing"))
  df$trait <- as.list(df$trait)
  expect_snapshot(error = TRUE, enrich_gc_gsea_kegg(df, by = "trait"))
})

test_that("enrich_gc_gsea_go forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
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
      formula = formula,
      data = data,
      fun = fun,
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
    sentinel
  }

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_go(
        .mock_gc_gsea_dea_df(),
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, clusterProfiler::gseGO)
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
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c("P01308", "P04637", "P00533", "P42345"),
      score = c(2, -3, 2.5, 0.5),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea_go builds trait-specific ranked lists from glystats results", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(formula, data, fun, ...) {
    captured <<- list(
      formula = formula,
      data = data,
      fun = fun,
      dots = list(...)
    )
    sentinel
  }

  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:6),
    trait = c("trait_A", "trait_A", "trait_A", "trait_B", "trait_B", "trait_B"),
    protein = c("P01308", "P01308", "P04637", "P42345", "P00533", "P00533"),
    p_val = c(0.01, 0.04, 0.001, 0.2, 0.03, 0.001),
    p_adj = c(0.02, 0.5, 0.001, 0.3, 0.04, 0.002),
    log2fc = c(2, -1, -3, 0.5, 1.5, 2.5),
    estimate = rep(1, 6)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_go(
        dea_res,
        rank_by = "signed_log10p",
        aggr = "max",
        orgdb = "org.Hs.eg.db"
      )
    )
  )

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, clusterProfiler::gseGO)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c("P01308", "P04637", "P00533", "P42345"),
      score = c(1.69897, -3, 2.69897, 0.5228787),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    ),
    tolerance = 1e-6
  )
})

test_that("enrich_gc_gsea_kegg forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
    organism = NULL,
    keyType = NULL,
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
      formula = formula,
      data = data,
      fun = fun,
      organism = organism,
      keyType = keyType,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_kegg(
        .mock_gc_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        organism = "mmu",
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, clusterProfiler::gseKEGG)
  expect_identical(captured$organism, "mmu")
  expect_identical(captured$keyType, "uniprot")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
})

test_that("enrich_gc_gsea_reactome forwards args to compareCluster through all layers", {
  skip_if_not_installed("ReactomePA")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
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
      formula = formula,
      data = data,
      fun = fun,
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
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
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
      enrich_gc_gsea_reactome(
        .mock_gc_gsea_dea_df(),
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, ReactomePA::gsePathway)
  expect_identical(captured$organism, "human")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c(
        "ENTREZ_P01308",
        "ENTREZ_P04637",
        "ENTREZ_P00533",
        "ENTREZ_P42345"
      ),
      score = c(2, -3, 2.5, 0.5),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea_wp forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
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
      formula = formula,
      data = data,
      fun = fun,
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
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
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
      enrich_gc_gsea_wp(
        .mock_gc_gsea_dea_df(),
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, clusterProfiler::gseWP)
  expect_identical(captured$organism, "Homo sapiens")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c(
        "ENTREZ_P01308",
        "ENTREZ_P04637",
        "ENTREZ_P00533",
        "ENTREZ_P42345"
      ),
      score = c(2, -3, 2.5, 0.5),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea_do forwards args to compareCluster through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
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
      formula = formula,
      data = data,
      fun = fun,
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
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
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
      enrich_gc_gsea_do(
        .mock_gc_gsea_dea_df(),
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, DOSE::gseDO)
  expect_identical(captured$ont, "HDO")
  expect_identical(captured$organism, "hsa")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c(
        "ENTREZ_P01308",
        "ENTREZ_P04637",
        "ENTREZ_P00533",
        "ENTREZ_P42345"
      ),
      score = c(2, -3, 2.5, 0.5),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea_ncg forwards args to compareCluster through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    formula,
    data,
    fun,
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
      formula = formula,
      data = data,
      fun = fun,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      exponent = exponent,
      eps = eps,
      seed = seed,
      dots = list(...)
    )
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
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
      enrich_gc_gsea_ncg(
        .mock_gc_gsea_dea_df(),
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

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, DOSE::gseNCG)
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
  expect_identical(captured$exponent, 2)
  expect_identical(captured$eps, 1e-8)
  expect_identical(captured$seed, FALSE)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c(
        "ENTREZ_P01308",
        "ENTREZ_P04637",
        "ENTREZ_P00533",
        "ENTREZ_P42345"
      ),
      score = c(2, -3, 2.5, 0.5),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea collapses duplicate Entrez IDs after conversion", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(formula, data, fun, ...) {
    captured <<- list(
      formula = formula,
      data = data,
      fun = fun,
      dots = list(...)
    )
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
  local_mocked_bindings(
    .reactome_orgdb = function(organism) paste0("MOCK_REACTOME_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb, drop_na = TRUE) {
      expect_identical(drop_na, FALSE)
      expect_identical(orgdb, "MOCK_REACTOME_human")
      dplyr::recode(
        uniprot,
        P01308 = "ENTREZ_DUP",
        P04637 = "ENTREZ_DUP",
        P00533 = "ENTREZ_DUP",
        P42345 = "ENTREZ_UNIQUE"
      )
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_reactome(
        .mock_gc_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        organism = "human"
      )
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c("ENTREZ_DUP", "ENTREZ_DUP", "ENTREZ_UNIQUE"),
      score = c(2, 2.5, 0.5),
      trait = c("trait_A", "trait_B", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea drops non-finite ranking scores before compareCluster", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(formula, data, fun, ...) {
    captured <<- list(
      formula = formula,
      data = data,
      fun = fun,
      dots = list(...)
    )
    sentinel
  }

  dea_res <- tibble::tibble(
    protein = c("P1", "P2", "P3", "P4"),
    site = c("s1", "s1", "s1", "s1"),
    trait = c("trait_A", "trait_A", "trait_B", "trait_B"),
    p_val = c(0, 0.01, NA, 0.001),
    log2fc = c(1, 2, -1, -3)
  )

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_go(
        dea_res,
        rank_by = "signed_log10p",
        aggr = "max",
        orgdb = "org.Hs.eg.db"
      )
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c("P2", "P4"),
      score = c(2, -3),
      trait = c("trait_A", "trait_B")
    )
  )
})

test_that("enrich_gc_gsea_go errors on data.frame with missing columns", {
  dea_res_missing_trait <- tibble::tibble(
    protein = c("P01308", "P04637"),
    site = c("s1", "s2"),
    p_val = c(0.001, 0.002),
    log2fc = c(2.5, -1.2)
  )
  expect_snapshot(
    error = TRUE,
    suppressMessages(enrich_gc_gsea_go(dea_res_missing_trait))
  )
})

test_that("enrich_gc_gsea_go errors when glycan trait columns are missing", {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:3),
    protein = c("P01308", "P04637", "P42345"),
    p_val = c(0.001, 0.002, 0.003),
    p_adj = c(0.001, 0.002, 0.003),
    log2fc = c(2.5, -1.2, 0.7),
    estimate = rep(1, 3)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  expect_error(
    enrich_gc_gsea_go(dea_res),
    "Cannot determine glycan traits"
  )
})

test_that("enrich_gc_gsea_go returns NULL when no terms are enriched", {
  local_mocked_bindings(
    compareCluster = function(...) NULL,
    .package = "clusterProfiler"
  )
  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_gc_gsea_go(
        .mock_gc_gsea_dea_df(),
        rank_by = "log2fc",
        aggr = "max",
        orgdb = "org.Hs.eg.db"
      )
    )
  )

  expect_null(result)
})
