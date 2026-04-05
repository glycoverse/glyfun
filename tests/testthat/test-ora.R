skip_if_not_installed("clusterProfiler")
skip_if_not_installed("org.Hs.eg.db")

# Build a minimal glystats_ttest_res with real UniProt IDs for integration tests
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

test_that("enrich_ora_go works with tibble input on happy path (integration)", {
  # Create a tibble with required columns
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )

  result <- suppressMessages(
    enrich_ora_go(dea_res, orgdb = "org.Hs.eg.db", ont = "MF", p_cutoff = 0.05)
  )

  expect_s3_class(
    result,
    c("glyfun_ora_go_res", "glyfun_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "id",
    "description",
    "gene_ratio",
    "bg_ratio",
    "rich_factor",
    "fold_enrichment",
    "z_score",
    "p_val",
    "p_adj",
    "q_val",
    "gene_id",
    "count"
  )
  expect_true(all(expected_cols %in% colnames(result$tidy_result)))
})

test_that("enrich_ora_go returns correct structure on happy path (integration)", {
  dea_res <- .mock_dea_res()

  result <- suppressMessages(
    enrich_ora_go(dea_res, orgdb = "org.Hs.eg.db", ont = "MF", p_cutoff = 0.05)
  )

  expect_s3_class(
    result,
    c("glyfun_ora_go_res", "glyfun_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "id",
    "description",
    "gene_ratio",
    "bg_ratio",
    "rich_factor",
    "fold_enrichment",
    "z_score",
    "p_val",
    "p_adj",
    "q_val",
    "gene_id",
    "count"
  )
  expect_true(all(expected_cols %in% colnames(result$tidy_result)))
})

test_that("enrich_ora_kegg returns correct structure on happy path (integration)", {
  skip_if_offline()

  dea_res <- .mock_dea_res()

  result <- try(
    suppressMessages(
      enrich_ora_kegg(dea_res, organism = "hsa", p_cutoff = 0.05)
    ),
    silent = TRUE
  )

  skip_if(
    inherits(result, "try-error"),
    "KEGG enrichment failed, likely due to network issues"
  )

  expect_s3_class(
    result,
    c("glyfun_ora_kegg_res", "glyfun_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "id",
    "description",
    "gene_ratio",
    "bg_ratio",
    "rich_factor",
    "fold_enrichment",
    "z_score",
    "p_val",
    "p_adj",
    "q_val",
    "gene_id",
    "count"
  )
  expect_true(all(expected_cols %in% colnames(result$tidy_result)))
})

# Error handling tests ----

test_that("enrich_ora_go errors on data.frame with missing columns", {
  # Missing 'protein' column
  dea_res_missing_protein <- tibble::tibble(
    site = rep("site1", 5),
    trait = rep("trait1", 5),
    p_val = rep(0.001, 5),
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
  )
  expect_error(
    enrich_ora_go(dea_res_missing_pval),
    "must have all expected columns"
  )

  # Missing 'log2FC' column
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
