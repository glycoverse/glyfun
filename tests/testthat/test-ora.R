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
