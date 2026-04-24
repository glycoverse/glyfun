skip_if_not_installed("clusterProfiler")

# Build a minimal glystats_ttest_res for contract/unit tests
.mock_dea_res <- function() {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:10),
    trait = rep(c("trait_A", "trait_B"), each = 5),
    protein = c(
      "P01308",
      "P04637",
      "P42345",
      "P00533",
      "P42336",
      "P01116",
      "P04049",
      "P22607",
      "P28482",
      "P24941"
    ),
    p_val = rep(0.001, 10),
    p_adj = rep(0.001, 10),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2, -2.5, -3.0, -1.5, -2.0, -2.2),
    estimate = rep(1, 10)
  )
  raw_result <- list()
  structure(
    list(tidy_result = tidy_result, raw_result = raw_result),
    class = c("glystats_ttest_res", "glystats_res")
  )
}

.mock_gc_dea_df <- function() {
  tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336", "P01116"),
    site = rep("site1", 6),
    trait = c("trait_A", "trait_A", "trait_B", "trait_B", "trait_B", "trait_B"),
    p_val = c(0.001, 0.2, 0.001, 0.001, 0.5, 0.001),
    log2fc = c(2.5, 3.0, 0.3, -2.0, 1.8, -0.2)
  )
}

test_that("enrich_gc_ora_go forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    keyType = NULL,
    OrgDb = NULL,
    ont = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      keyType = keyType,
      OrgDb = OrgDb,
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
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
    enrich_gc_ora_go(
      .mock_gc_dea_df(),
      orgdb = "org.Hs.eg.db",
      ont = "BP",
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 3,
      max_gs_size = 99
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "P01308", trait_B = "P00533")
  )
  expect_identical(captured$fun, clusterProfiler::enrichGO)
  expect_equal(captured$universe, c("P01308", "P04637", "P42345"))
  expect_identical(captured$keyType, "UNIPROT")
  expect_identical(captured$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$ont, "BP")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 3)
  expect_identical(captured$maxGSSize, 99)
})

test_that(".gc_ora.data.frame groups proteins by trait after filtering", {
  sentinel <- list(source = "gc_ora_impl")
  captured <- NULL
  mock_gc_ora_impl <- function(
    dea_res,
    enrich_fun,
    dea_p_cutoff,
    dea_log2fc_cutoff,
    ...,
    pro_list_fun = NULL,
    uniprot_to_entrez = FALSE
  ) {
    captured <<- list(
      protein_list = pro_list_fun(dea_res),
      dea_p_cutoff = dea_p_cutoff,
      dea_log2fc_cutoff = dea_log2fc_cutoff,
      uniprot_to_entrez = uniprot_to_entrez,
      dots = list(...)
    )
    sentinel
  }
  local_mocked_bindings(.gc_ora_impl = mock_gc_ora_impl, .package = "glyfun")

  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = c("trait_A", "trait_A", "trait_B", "trait_B", "trait_B"),
    p_val = c(0.001, 0.2, 0.001, 0.001, 0.5),
    log2fc = c(2.5, 3.0, 0.3, -2.0, 1.8)
  )

  result <- glyfun:::.gc_ora(
    dea_res = dea_res,
    enrich_fun = function(...) list(source = "mock_enrich"),
    dea_p_cutoff = 0.05,
    dea_log2fc_cutoff = c(-1, 1)
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$protein_list,
    list(trait_A = "P01308", trait_B = "P00533")
  )
  expect_identical(captured$dea_p_cutoff, 0.05)
  expect_equal(captured$dea_log2fc_cutoff, c(-1, 1))
  expect_false(captured$uniprot_to_entrez)
})

test_that(".gc_ora.glystats_res groups proteins by trait after filtering", {
  sentinel <- list(source = "gc_ora_impl")
  captured <- NULL
  mock_gc_ora_impl <- function(
    dea_res,
    enrich_fun,
    dea_p_cutoff,
    dea_log2fc_cutoff,
    ...,
    pro_list_fun = NULL,
    uniprot_to_entrez = FALSE
  ) {
    captured <<- list(
      protein_list = pro_list_fun(dea_res),
      dea_p_cutoff = dea_p_cutoff,
      dea_log2fc_cutoff = dea_log2fc_cutoff,
      uniprot_to_entrez = uniprot_to_entrez,
      dots = list(...)
    )
    sentinel
  }
  local_mocked_bindings(.gc_ora_impl = mock_gc_ora_impl, .package = "glyfun")

  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:6),
    trait = c("trait_A", "trait_A", "trait_A", "trait_B", "trait_B", "trait_B"),
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336", "P01116"),
    p_val = rep(0.001, 6),
    p_adj = c(0.001, 0.2, 0.001, 0.001, 0.001, 0.2),
    log2fc = c(2.5, 3.0, 0.3, -2.0, -0.2, -3.0),
    estimate = rep(1, 6)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  result <- glyfun:::.gc_ora(
    dea_res = dea_res,
    enrich_fun = function(...) list(source = "mock_enrich"),
    dea_p_cutoff = 0.05,
    dea_log2fc_cutoff = c(-1, 1)
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$protein_list,
    list(trait_A = "P01308", trait_B = "P00533")
  )
  expect_identical(captured$dea_p_cutoff, 0.05)
  expect_equal(captured$dea_log2fc_cutoff, c(-1, 1))
  expect_false(captured$uniprot_to_entrez)
})

test_that("enrich_gc_ora_kegg forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    keyType = NULL,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      keyType = keyType,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      dots = list(...)
    )
    sentinel
  }
  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    enrich_gc_ora_kegg(
      .mock_gc_dea_df(),
      organism = "mmu",
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "P01308", trait_B = "P00533")
  )
  expect_identical(captured$fun, clusterProfiler::enrichKEGG)
  expect_equal(captured$universe, c("P01308", "P04637", "P42345"))
  expect_identical(captured$keyType, "uniprot")
  expect_identical(captured$organism, "mmu")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

test_that("enrich_gc_ora_reactome forwards args to compareCluster through all layers", {
  skip_if_not_installed("ReactomePA")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
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
    .uniprot_to_entrez_prolist = function(pro_list, orgdb) {
      expect_identical(orgdb, "MOCK_REACTOME_human")
      purrr::map(pro_list, ~ paste0("ENTREZ_", .x))
    },
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_REACTOME_human")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_gc_ora_reactome(
      .mock_gc_dea_df(),
      organism = "human",
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "ENTREZ_P01308", trait_B = "ENTREZ_P00533")
  )
  expect_identical(captured$fun, ReactomePA::enrichPathway)
  expect_equal(
    captured$universe,
    c("ENTREZ_P01308", "ENTREZ_P04637", "ENTREZ_P42345")
  )
  expect_identical(captured$organism, "human")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

test_that("enrich_gc_ora_wp forwards args to compareCluster through all layers", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
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
    .uniprot_to_entrez_prolist = function(pro_list, orgdb) {
      expect_identical(orgdb, "MOCK_WP_Homo sapiens")
      purrr::map(pro_list, ~ paste0("ENTREZ_", .x))
    },
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_WP_Homo sapiens")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_gc_ora_wp(
      .mock_gc_dea_df(),
      organism = "Homo sapiens",
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "ENTREZ_P01308", trait_B = "ENTREZ_P00533")
  )
  expect_identical(captured$fun, clusterProfiler::enrichWP)
  expect_equal(
    captured$universe,
    c("ENTREZ_P01308", "ENTREZ_P04637", "ENTREZ_P42345")
  )
  expect_identical(captured$organism, "Homo sapiens")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

test_that("enrich_gc_ora_do forwards args to compareCluster through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    ont = NULL,
    organism = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      ont = ont,
      organism = organism,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
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
    .uniprot_to_entrez_prolist = function(pro_list, orgdb) {
      expect_identical(orgdb, "MOCK_DO_hsa")
      purrr::map(pro_list, ~ paste0("ENTREZ_", .x))
    },
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_DO_hsa")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_gc_ora_do(
      .mock_gc_dea_df(),
      ont = "HDO",
      organism = "hsa",
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "ENTREZ_P01308", trait_B = "ENTREZ_P00533")
  )
  expect_identical(captured$fun, DOSE::enrichDO)
  expect_equal(
    captured$universe,
    c("ENTREZ_P01308", "ENTREZ_P04637", "ENTREZ_P42345")
  )
  expect_identical(captured$ont, "HDO")
  expect_identical(captured$organism, "hsa")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

test_that("enrich_gc_ora_ncg forwards args to compareCluster through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(
    geneCluster,
    fun,
    universe = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      geneCluster = geneCluster,
      fun = fun,
      universe = universe,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
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
    .uniprot_to_entrez_prolist = function(pro_list, orgdb) {
      expect_identical(orgdb, "MOCK_org.Hs.eg.db")
      purrr::map(pro_list, ~ paste0("ENTREZ_", .x))
    },
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_org.Hs.eg.db")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_gc_ora_ncg(
      .mock_gc_dea_df(),
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(trait_A = "ENTREZ_P01308", trait_B = "ENTREZ_P00533")
  )
  expect_identical(captured$fun, DOSE::enrichNCG)
  expect_equal(
    captured$universe,
    c("ENTREZ_P01308", "ENTREZ_P04637", "ENTREZ_P42345")
  )
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

# Error handling tests ----

test_that("enrich_gc_ora_go errors on data.frame with missing columns", {
  # Missing 'trait' column
  dea_res_missing_trait <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_gc_ora_go(dea_res_missing_trait),
    "must have all expected columns"
  )

  # Missing 'protein' column
  dea_res_missing_protein <- tibble::tibble(
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_gc_ora_go(dea_res_missing_protein),
    "must have all expected columns"
  )
})

test_that("enrich_gc_ora_go errors on invalid dea_p_cutoff", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep(c("trait_A", "trait_B"), c(3, 2)),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )

  # Negative p_cutoff
  expect_error(
    enrich_gc_ora_go(dea_res, dea_p_cutoff = -0.1),
    "Assertion on 'p_cutoff' failed"
  )

  # p_cutoff > 1
  expect_error(
    enrich_gc_ora_go(dea_res, dea_p_cutoff = 1.5),
    "Assertion on 'p_cutoff' failed"
  )

  # Non-numeric p_cutoff
  expect_error(
    enrich_gc_ora_go(dea_res, dea_p_cutoff = "invalid"),
    "Assertion on 'p_cutoff' failed"
  )

  # NA value
  expect_error(
    enrich_gc_ora_go(dea_res, dea_p_cutoff = NA),
    "Assertion on 'p_cutoff' failed"
  )

  # NaN value
  expect_error(
    enrich_gc_ora_go(dea_res, dea_p_cutoff = NaN),
    "Assertion on 'p_cutoff' failed"
  )
})

test_that("enrich_gc_ora_go errors on invalid dea_log2fc_cutoff", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep(c("trait_A", "trait_B"), c(3, 2)),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )

  # Wrong length (1 element)
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = c(-1)),
    "must have exactly 2 elements"
  )

  # Wrong length (3 elements)
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = c(-1, 0, 1)),
    "must have exactly 2 elements"
  )

  # First element positive
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = c(0.5, 1)),
    "must be 0 or negative"
  )

  # Second element negative
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = c(-1, -0.5)),
    "must be 0 or positive"
  )

  # Non-numeric
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = "invalid"),
    "Assertion on 'log2fc_cutoff' failed"
  )

  # NULL value
  expect_error(
    enrich_gc_ora_go(dea_res, dea_log2fc_cutoff = NULL),
    "Assertion on 'log2fc_cutoff' failed"
  )
})

test_that("enrich_gc_ora_go returns NULL when no terms are enriched", {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:3),
    trait = rep("trait_A", 3),
    protein = c("P01308", "P04637", "P42345"),
    p_val = rep(0.001, 3),
    p_adj = rep(0.001, 3),
    log2fc = c(2.5, 3.0, 1.5),
    estimate = rep(1, 3)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  # Mock downstream dependencies to keep this test fully offline.
  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  mock_compare_cluster <- function(...) NULL
  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )

  result <- suppressMessages(
    suppressWarnings(enrich_gc_ora_go(
      dea_res,
      orgdb = "org.Hs.eg.db",
      ont = "MF"
    ))
  )

  expect_null(result)
})

test_that(".gc_ora_impl converts protein list and universe when uniprot_to_entrez is TRUE", {
  sentinel <- list(source = "compareCluster")
  captured <- NULL
  mock_compare_cluster <- function(geneCluster, universe = NULL, ...) {
    captured <<- list(
      geneCluster = geneCluster,
      universe = universe,
      dots = list(...)
    )
    sentinel
  }
  mock_uniprot_to_entrez_prolist <- function(pro_list, orgdb) {
    expect_identical(orgdb, "MOCK_ORGDB")
    purrr::map(pro_list, ~ paste0("ENTREZ_", .x))
  }
  mock_uniprot_to_entrez <- function(uniprot, orgdb) {
    expect_identical(orgdb, "MOCK_ORGDB")
    paste0("ENTREZ_", uniprot)
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
  local_mocked_bindings(
    .uniprot_to_entrez_prolist = mock_uniprot_to_entrez_prolist,
    .uniprot_to_entrez = mock_uniprot_to_entrez,
    .package = "glyfun"
  )

  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637"),
    site = c("s1", "s2"),
    trait = c("t1", "t2"),
    p_val = c(0.001, 0.001),
    log2fc = c(2, -2)
  )

  result <- suppressMessages(
    glyfun:::.gc_ora_impl(
      dea_res = dea_res,
      enrich_fun = function(...) list(source = "mock_enrich"),
      universe = c("U1", "U2"),
      bitr_orgdb = "MOCK_ORGDB",
      uniprot_to_entrez = TRUE,
      pro_list_fun = function(dea_res) {
        expect_equal(nrow(dea_res), 2)
        list(trait_A = c("P01308", "P04637"), trait_B = "P00533")
      }
    )
  )

  expect_identical(result, sentinel)
  expect_equal(
    captured$geneCluster,
    list(
      trait_A = c("ENTREZ_P01308", "ENTREZ_P04637"),
      trait_B = "ENTREZ_P00533"
    )
  )
  expect_equal(captured$universe, c("ENTREZ_U1", "ENTREZ_U2"))
})

test_that(".gc_ora_impl errors when uniprot_to_entrez is TRUE but bitr_orgdb is missing", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637"),
    site = c("s1", "s2"),
    trait = c("t1", "t2"),
    p_val = c(0.001, 0.001),
    log2fc = c(2, -2)
  )

  expect_error(
    glyfun:::.gc_ora_impl(
      dea_res = dea_res,
      enrich_fun = function(...) list(source = "mock_enrich"),
      uniprot_to_entrez = TRUE,
      pro_list_fun = function(dea_res) {
        list(trait_A = dea_res$protein)
      }
    ),
    "bitr_orgdb"
  )
})

test_that("enrich_gc_ora_go errors when glycan trait columns are missing", {
  # Create a glystats_res without glycan_structure, glycan_composition, trait, or motif columns
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:5),
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    p_val = rep(0.001, 5),
    p_adj = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2),
    estimate = rep(1, 5)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  expect_error(
    enrich_gc_ora_go(dea_res),
    "Cannot determine glycan traits"
  )
})
