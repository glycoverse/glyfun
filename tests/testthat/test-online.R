skip_if_not_installed("clusterProfiler")

.mock_online_ora_dea_res <- function() {
  tibble::tibble(
    protein = c(
      "P01308",
      "P04637",
      "P42345",
      "P00533",
      "P42336",
      "P01116",
      "P04049",
      "P22607"
    ),
    site = rep("site1", 8),
    trait = rep(c("trait_A", "trait_B"), each = 4),
    p_val = rep(0.001, 8),
    log2FC = c(2.5, 3.0, 1.5, 2.0, -2.5, -3.0, -1.5, -2.0)
  )
}

.skip_online_integration <- function() {
  skip_if_online_tests_disabled()
  skip_if_offline()
  skip_on_cran()
  skip_on_ci()
}

.expect_s4_or_null_with_external_skip <- function(expr, expected_s4) {
  result <- tryCatch(
    suppressMessages(suppressWarnings(force(expr))),
    error = identity
  )

  if (inherits(result, "error")) {
    msg <- tolower(conditionMessage(result))
    external_error_patterns <- c(
      "cannot open url",
      "timeout",
      "timed out",
      "connection",
      "http",
      "https",
      "ssl",
      "download",
      "service unavailable",
      "internal server error",
      "bad gateway",
      "network"
    )
    if (any(vapply(external_error_patterns, grepl, logical(1), x = msg, fixed = TRUE))) {
      skip(paste("External service/network issue:", conditionMessage(result)))
    }
    stop(result)
  }

  expect_true(is.null(result) || methods::is(result, expected_s4))
}

test_that("online integration smoke: enrich_ora_go returns enrichResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_go(dea_res, orgdb = "org.Hs.eg.db", ont = "MF", p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_ora_kegg returns enrichResult or NULL", {
  .skip_online_integration()
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_kegg(dea_res, organism = "hsa", p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_ora_reactome returns enrichResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("ReactomePA")
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_reactome(dea_res, organism = "human", p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_ora_wp returns enrichResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_wp(dea_res, organism = "Homo sapiens", p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_ora_do returns enrichResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("DOSE")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_no_hdo()
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_do(dea_res, ont = "HDO", organism = "hsa", p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_ora_ncg returns enrichResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("DOSE")
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_ora_ncg(dea_res, p_cutoff = 0.05),
    "enrichResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_go returns compareClusterResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_go(dea_res, orgdb = "org.Hs.eg.db", ont = "MF", p_cutoff = 0.05),
    "compareClusterResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_kegg returns compareClusterResult or NULL", {
  .skip_online_integration()
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_kegg(dea_res, organism = "hsa", p_cutoff = 0.05),
    "compareClusterResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_reactome returns compareClusterResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("ReactomePA")
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_reactome(dea_res, organism = "human", p_cutoff = 0.05),
    "compareClusterResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_wp returns compareClusterResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_wp(dea_res, organism = "Homo sapiens", p_cutoff = 0.05),
    "compareClusterResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_do returns compareClusterResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("DOSE")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_no_hdo()
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_do(dea_res, ont = "HDO", organism = "hsa", p_cutoff = 0.05),
    "compareClusterResult"
  )
})

test_that("online integration smoke: enrich_gc_ora_ncg returns compareClusterResult or NULL", {
  .skip_online_integration()
  skip_if_not_installed("DOSE")
  skip_if_not_installed("org.Hs.eg.db")
  dea_res <- .mock_online_ora_dea_res()
  .expect_s4_or_null_with_external_skip(
    enrich_gc_ora_ncg(dea_res, p_cutoff = 0.05),
    "compareClusterResult"
  )
})
