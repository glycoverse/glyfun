skip_if_not_installed("clusterProfiler")

.mock_unified_dea_df <- function() {
  tibble::tibble(
    protein = c("P1", "P2", "P3", "P4"),
    site = rep("site1", 4),
    trait = c("trait_A", "trait_A", "trait_B", "trait_B"),
    p_val = c(0.001, 0.2, 0.001, 0.001),
    log2fc = c(2, 0.5, -2, 0.2)
  )
}

test_that("enrich_glyco routes protein-centric ORA to the selected database", {
  captured <- NULL
  sentinel <- list(source = "enrichGO")
  mock_enrich_go <- function(gene, universe = NULL, ...) {
    captured <<- list(gene = gene, universe = universe, dots = list(...))
    sentinel
  }

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(enrichGO = mock_enrich_go, .package = "clusterProfiler")

  result <- enrich_glyco(
    .mock_unified_dea_df(),
    method = "ora",
    database = "GO",
    scope = "protein",
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    universe = c("P1", "P2", "P3")
  )

  expect_identical(result, sentinel)
  expect_identical(captured$gene, c("P1", "P3"))
  expect_identical(captured$universe, c("P1", "P2", "P3"))
  expect_identical(captured$dots$keyType, "UNIPROT")
  expect_identical(captured$dots$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$dots$ont, "BP")
})

test_that("enrich_glyco routes glycan-trait ORA through compareCluster", {
  captured <- NULL
  sentinel <- list(source = "compareCluster")
  mock_compare_cluster <- function(geneCluster, fun, ...) {
    captured <<- list(geneCluster = geneCluster, fun = fun, dots = list(...))
    sentinel
  }

  local_mocked_bindings(
    compareCluster = mock_compare_cluster,
    .package = "clusterProfiler"
  )
  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )

  result <- suppressMessages(
    enrich_glyco(
      .mock_unified_dea_df(),
      method = "ora",
      database = "GO",
      scope = "glycan_trait",
      orgdb = "org.Hs.eg.db",
      ont = "BP"
    )
  )

  expect_identical(result, sentinel)
  expect_identical(captured$fun, clusterProfiler::enrichGO)
  expect_equal(captured$geneCluster, list(trait_A = "P1", trait_B = "P3"))
  expect_identical(captured$dots$keyType, "UNIPROT")
  expect_identical(captured$dots$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$dots$ont, "BP")
})

test_that("enrich_glyco routes protein-centric GSEA to the selected database", {
  captured <- NULL
  sentinel <- tibble::tibble(ID = "GO:1")
  mock_gse_go <- function(geneList, ...) {
    captured <<- list(geneList = geneList, dots = list(...))
    sentinel
  }

  local_mocked_bindings(
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )
  local_mocked_bindings(gseGO = mock_gse_go, .package = "clusterProfiler")

  result <- suppressWarnings(
    enrich_glyco(
      .mock_unified_dea_df(),
      method = "gsea",
      database = "GO",
      scope = "protein",
      rank_by = "log2fc",
      aggr = "max",
      orgdb = "org.Hs.eg.db",
      ont = "BP"
    )
  )

  expect_identical(result, sentinel)
  expect_equal(unname(captured$geneList), c(2, 0.5, 0.2, -2))
  expect_identical(names(captured$geneList), c("P1", "P2", "P4", "P3"))
  expect_identical(captured$dots$keyType, "UNIPROT")
  expect_identical(captured$dots$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$dots$ont, "BP")
})

test_that("enrich_glyco routes glycan-trait GSEA through compareCluster", {
  captured <- NULL
  sentinel <- list(source = "compareCluster")
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
    .prepare_orgdb = function(orgdb) paste0("MOCK_", orgdb),
    .package = "glyfun"
  )

  result <- suppressMessages(
    suppressWarnings(
      enrich_glyco(
        .mock_unified_dea_df(),
        method = "gsea",
        database = "GO",
        scope = "glycan_trait",
        rank_by = "log2fc",
        aggr = "max",
        orgdb = "org.Hs.eg.db",
        ont = "BP"
      )
    )
  )

  expect_identical(result, sentinel)
  expect_identical(deparse(captured$formula), "gene | score ~ trait")
  expect_identical(captured$fun, clusterProfiler::gseGO)
  expect_equal(
    captured$data,
    tibble::tibble(
      gene = c("P1", "P2", "P4", "P3"),
      score = c(2, 0.5, 0.2, -2),
      trait = c("trait_A", "trait_A", "trait_B", "trait_B")
    )
  )
  expect_identical(captured$dots$keyType, "UNIPROT")
  expect_identical(captured$dots$OrgDb, "MOCK_org.Hs.eg.db")
  expect_identical(captured$dots$ont, "BP")
})
