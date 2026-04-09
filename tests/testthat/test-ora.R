skip_if_not_installed("clusterProfiler")

# Build a minimal glystats_ttest_res for contract/unit tests
.mock_dea_res <- function() {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:10),
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

.mock_dea_df <- function() {
  tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336", "P01116"),
    site = rep("site1", 6),
    trait = rep("trait1", 6),
    p_val = c(0.001, 0.001, 0.2, 0.001, 0.001, 0.9),
    log2fc = c(2.5, -3.0, 1.5, 0.8, -0.5, -2.2)
  )
}

test_that("enrich_ora_go forwards args to enrichGO through all layers", {
  # We capture arguments in the parent test scope and return a fixed sentinel
  # so this test can independently verify both dimensions:
  # (1) all wrapper/dispatch layers forward parameters correctly, and
  # (2) the top-level API returns the downstream object unchanged.
  sentinel <- list(source = "enrichGO")
  captured <- NULL
  mock_enrich_go <- function(
    gene,
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
      gene = gene,
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
  local_mocked_bindings(enrichGO = mock_enrich_go, .package = "clusterProfiler")

  result <- suppressMessages(
    enrich_ora_go(
      .mock_dea_df(),
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
  expect_equal(captured$gene, c("P01308", "P04637"))
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

test_that(".ora.data.frame filters proteins by p_val and log2fc", {
  sentinel <- list(source = "ora_impl")
  captured <- NULL
  mock_ora_impl <- function(
    dea_res,
    enrich_fun,
    result_class,
    dea_p_cutoff,
    dea_log2fc_cutoff,
    ...,
    pro_fun = NULL,
    uniprot_to_entrez = FALSE
  ) {
    captured <<- list(
      proteins = pro_fun(dea_res),
      dea_p_cutoff = dea_p_cutoff,
      dea_log2fc_cutoff = dea_log2fc_cutoff,
      uniprot_to_entrez = uniprot_to_entrez,
      dots = list(...)
    )
    sentinel
  }
  local_mocked_bindings(.ora_impl = mock_ora_impl, .package = "glyfun")

  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = c(0.001, 0.001, 0.2, 0.001, 0.001),
    log2fc = c(2.5, -3.0, 1.5, 0.8, -0.5)
  )

  result <- glyfun:::.ora(
    dea_res = dea_res,
    enrich_fun = function(...) list(source = "mock_enrich"),
    result_class = "unused",
    dea_p_cutoff = 0.05,
    dea_log2fc_cutoff = c(-1, 1)
  )

  expect_identical(result, sentinel)
  expect_equal(captured$proteins, c("P01308", "P04637"))
  expect_identical(captured$dea_p_cutoff, 0.05)
  expect_equal(captured$dea_log2fc_cutoff, c(-1, 1))
  expect_false(captured$uniprot_to_entrez)
})

test_that(".ora.glystats_res filters proteins by p_adj and log2fc", {
  sentinel <- list(source = "ora_impl")
  captured <- NULL
  mock_ora_impl <- function(
    dea_res,
    enrich_fun,
    result_class,
    dea_p_cutoff,
    dea_log2fc_cutoff,
    ...,
    pro_fun = NULL,
    uniprot_to_entrez = FALSE
  ) {
    captured <<- list(
      proteins = pro_fun(dea_res),
      dea_p_cutoff = dea_p_cutoff,
      dea_log2fc_cutoff = dea_log2fc_cutoff,
      uniprot_to_entrez = uniprot_to_entrez,
      dots = list(...)
    )
    sentinel
  }
  local_mocked_bindings(.ora_impl = mock_ora_impl, .package = "glyfun")

  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:5),
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    p_val = rep(0.001, 5),
    p_adj = c(0.001, 0.2, 0.001, 0.001, 0.01),
    log2fc = c(2.5, 3.0, 0.3, -2.0, 0.2),
    estimate = rep(1, 5)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
  )

  result <- glyfun:::.ora(
    dea_res = dea_res,
    enrich_fun = function(...) list(source = "mock_enrich"),
    result_class = "unused",
    dea_p_cutoff = 0.05,
    dea_log2fc_cutoff = c(-1, 1)
  )

  expect_identical(result, sentinel)
  expect_equal(captured$proteins, c("P01308", "P00533"))
  expect_identical(captured$dea_p_cutoff, 0.05)
  expect_equal(captured$dea_log2fc_cutoff, c(-1, 1))
  expect_false(captured$uniprot_to_entrez)
})

test_that("enrich_ora_kegg forwards args to enrichKEGG through all layers", {
  sentinel <- list(source = "enrichKEGG")
  captured <- NULL
  mock_enrich_kegg <- function(
    gene,
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
      gene = gene,
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
  local_mocked_bindings(enrichKEGG = mock_enrich_kegg, .package = "clusterProfiler")

  result <- suppressMessages(
    enrich_ora_kegg(
      .mock_dea_df(),
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
  expect_equal(captured$gene, c("P01308", "P04637"))
  expect_equal(captured$universe, c("P01308", "P04637", "P42345"))
  expect_identical(captured$keyType, "uniprot")
  expect_identical(captured$organism, "mmu")
  expect_identical(captured$pAdjustMethod, "BY")
  expect_identical(captured$pvalueCutoff, 0.01)
  expect_identical(captured$qvalueCutoff, 0.1)
  expect_identical(captured$minGSSize, 5)
  expect_identical(captured$maxGSSize, 50)
})

test_that("enrich_ora_reactome forwards args to enrichPathway through all layers", {
  skip_if_not_installed("ReactomePA")

  sentinel <- list(source = "enrichPathway")
  captured <- NULL
  mock_enrich_pathway <- function(
    gene,
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
      gene = gene,
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

  local_mocked_bindings(enrichPathway = mock_enrich_pathway, .package = "ReactomePA")
  local_mocked_bindings(
    .reactome_orgdb = function(organism) paste0("MOCK_REACTOME_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_REACTOME_human")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_ora_reactome(
      .mock_dea_df(),
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
  expect_equal(captured$gene, c("ENTREZ_P01308", "ENTREZ_P04637"))
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

test_that("enrich_ora_wp forwards args to enrichWP through all layers", {
  sentinel <- list(source = "enrichWP")
  captured <- NULL
  mock_enrich_wp <- function(
    gene,
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
      gene = gene,
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

  local_mocked_bindings(enrichWP = mock_enrich_wp, .package = "clusterProfiler")
  local_mocked_bindings(
    .wp_orgdb = function(organism) paste0("MOCK_WP_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_WP_Homo sapiens")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_ora_wp(
      .mock_dea_df(),
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
  expect_equal(captured$gene, c("ENTREZ_P01308", "ENTREZ_P04637"))
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

test_that("enrich_ora_do forwards args to enrichDO through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "enrichDO")
  captured <- NULL
  mock_enrich_do <- function(
    gene,
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
      gene = gene,
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

  local_mocked_bindings(enrichDO = mock_enrich_do, .package = "DOSE")
  local_mocked_bindings(
    .do_orgdb = function(organism) paste0("MOCK_DO_", organism),
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_DO_hsa")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_ora_do(
      .mock_dea_df(),
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
  expect_equal(captured$gene, c("ENTREZ_P01308", "ENTREZ_P04637"))
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

test_that("enrich_ora_ncg forwards args to enrichNCG through all layers", {
  skip_if_not_installed("DOSE")

  sentinel <- list(source = "enrichNCG")
  captured <- NULL
  mock_enrich_ncg <- function(
    gene,
    universe = NULL,
    pAdjustMethod = NULL,
    pvalueCutoff = NULL,
    qvalueCutoff = NULL,
    minGSSize = NULL,
    maxGSSize = NULL,
    ...
  ) {
    captured <<- list(
      gene = gene,
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

  local_mocked_bindings(enrichNCG = mock_enrich_ncg, .package = "DOSE")
  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .uniprot_to_entrez = function(uniprot, orgdb) {
      expect_identical(orgdb, "MOCK_org.Hs.eg.db")
      paste0("ENTREZ_", uniprot)
    },
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_ora_ncg(
      .mock_dea_df(),
      universe = c("P01308", "P04637", "P42345"),
      p_adj_method = "BY",
      p_cutoff = 0.01,
      q_cutoff = 0.1,
      min_gs_size = 5,
      max_gs_size = 50
    )
  )

  expect_identical(result, sentinel)
  expect_equal(captured$gene, c("ENTREZ_P01308", "ENTREZ_P04637"))
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

test_that("enrich_ora_go errors on data.frame with missing columns", {
  # Missing 'protein' column
  dea_res_missing_protein <- tibble::tibble(
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_protein),
    "must have all expected columns"
  )

  # Missing 'site' column
  dea_res_missing_site <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_site),
    "must have all expected columns"
  )

  # Missing 'trait' column
  dea_res_missing_trait <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_trait),
    "must have all expected columns"
  )

  # Missing 'p_val' column
  dea_res_missing_pval <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_pval),
    "must have all expected columns"
  )

  # Missing 'log2fc' column
  dea_res_missing_log2fc <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_log2fc),
    "must have all expected columns"
  )

  # Missing multiple columns
  dea_res_missing_multiple <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    p_val = rep(0.001, 5)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_multiple),
    "must have all expected columns"
  )
})

test_that("enrich_ora_go errors on invalid dea_p_cutoff", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )

  # Negative p_cutoff
  expect_error(
    enrich_ora_go(dea_res, dea_p_cutoff = -0.1),
    "Assertion on 'p_cutoff' failed"
  )

  # p_cutoff > 1
  expect_error(
    enrich_ora_go(dea_res, dea_p_cutoff = 1.5),
    "Assertion on 'p_cutoff' failed"
  )

  # Non-numeric p_cutoff
  expect_error(
    enrich_ora_go(dea_res, dea_p_cutoff = "invalid"),
    "Assertion on 'p_cutoff' failed"
  )

  # Multiple values
  expect_error(
    enrich_ora_go(dea_res, dea_p_cutoff = c(0.05, 0.01)),
    "Assertion on 'p_cutoff' failed"
  )
})

test_that("enrich_ora_go errors on invalid dea_log2fc_cutoff", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )

  # Wrong length
  expect_error(
    enrich_ora_go(dea_res, dea_log2fc_cutoff = c(-1)),
    "must have exactly 2 elements"
  )

  expect_error(
    enrich_ora_go(dea_res, dea_log2fc_cutoff = c(-1, 1, 2)),
    "must have exactly 2 elements"
  )

  # First element positive
  expect_error(
    enrich_ora_go(dea_res, dea_log2fc_cutoff = c(1, 2)),
    "must be 0 or negative"
  )

  # Second element negative
  expect_error(
    enrich_ora_go(dea_res, dea_log2fc_cutoff = c(-2, -1)),
    "must be 0 or positive"
  )

  # Non-numeric
  expect_error(
    enrich_ora_go(dea_res, dea_log2fc_cutoff = c("a", "b")),
    "Assertion on 'log2fc_cutoff' failed"
  )
})

test_that("enrich_ora_go errors on unsupported glystats_res class", {
  dea_res <- list(tidy_res = tibble::tibble(), raw_res = list())
  class(dea_res) <- c("glystats_invalid_res", "glystats_res")

  expect_error(
    enrich_ora_go(dea_res),
    "Unsupported input class for `dea_res`"
  )
})

test_that("enrich_ora_go errors on glystats_res with invalid experiment type", {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:5),
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    p_val = rep(0.001, 5),
    p_adj = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2),
    estimate = rep(1, 5)
  )
  dea_res <- structure(
    list(
      tidy_result = tidy_result,
      raw_result = list(),
      meta_data = list(exp_type = "invalid_type")
    ),
    class = c("glystats_ttest_res", "glystats_res")
  )

  expect_error(
    enrich_ora_go(dea_res),
    "must be of"
  )
})

test_that("enrich_ora_go errors on glystats_res without protein column", {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:5),
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
    enrich_ora_go(dea_res),
    "A protein column must be in `dea_res`"
  )
})

test_that("enrich_ora_go errors on multi-group glystats_limma_res", {
  # Create a glystats_limma_res mock with multi-group contrasts
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:10),
    protein = rep(c("P01308", "P04637"), 5),
    p_val = rep(0.001, 10),
    p_adj = rep(0.001, 10),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2, -2.5, -3.0, -1.5, -2.0, -2.2),
    estimate = rep(1, 10),
    ref_group = rep(c("A", "B"), each = 5),
    test_group = rep(c("C", "D"), each = 5)
  )
  # Use tidy_result (not tidy_res) to match glystats convention
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_limma_res", "glystats_res")
  )

  expect_error(
    enrich_ora_go(dea_res),
    "does not support multi-group"
  )
})

test_that("enrich_ora_go returns NULL when no terms are enriched", {
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:3),
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
  mock_enrich_go <- function(...) NULL
  local_mocked_bindings(enrichGO = mock_enrich_go, .package = "clusterProfiler")

  result <- suppressMessages(
    suppressWarnings(enrich_ora_go(dea_res, orgdb = "org.Hs.eg.db", ont = "MF"))
  )

  expect_null(result)
})

test_that(".ora_impl converts proteins and universe when uniprot_to_entrez is TRUE", {
  captured <- NULL
  mock_enrich <- function(gene, universe = NULL, ...) {
    captured <<- list(gene = gene, universe = universe, dots = list(...))
    list(source = "mock_enrich")
  }
  mock_uniprot_to_entrez <- function(uniprot, orgdb) {
    expect_identical(orgdb, "MOCK_ORGDB")
    paste0("ENTREZ_", uniprot)
  }
  local_mocked_bindings(
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

  result <- glyfun:::.ora_impl(
    dea_res = dea_res,
    enrich_fun = mock_enrich,
    result_class = "unused",
    universe = c("U1", "U2"),
    bitr_orgdb = "MOCK_ORGDB",
    uniprot_to_entrez = TRUE,
    pro_fun = function(dea_res) {
      expect_equal(nrow(dea_res), 2)
      c("P01308", "P04637")
    }
  )

  expect_equal(captured$gene, c("ENTREZ_P01308", "ENTREZ_P04637"))
  expect_equal(captured$universe, c("ENTREZ_U1", "ENTREZ_U2"))
  expect_identical(result$source, "mock_enrich")
})

test_that(".ora_impl errors when uniprot_to_entrez is TRUE but bitr_orgdb is missing", {
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637"),
    site = c("s1", "s2"),
    trait = c("t1", "t2"),
    p_val = c(0.001, 0.001),
    log2fc = c(2, -2)
  )

  expect_error(
    glyfun:::.ora_impl(
      dea_res = dea_res,
      enrich_fun = function(...) list(source = "mock_enrich"),
      result_class = "unused",
      uniprot_to_entrez = TRUE,
      pro_fun = function(dea_res) {
        dea_res$protein
      }
    ),
    "bitr_orgdb"
  )
})
