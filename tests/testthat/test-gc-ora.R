skip_if_not_installed("clusterProfiler")
skip_if_not_installed("org.Hs.eg.db")

# Build a minimal glystats_ttest_res with real UniProt IDs for integration tests
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

test_that("enrich_gc_ora_go works with tibble input on happy path (integration)", {
  # Create a tibble with required columns
  dea_res <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336", "P01116"),
    site = rep("site1", 6),
    trait = rep(c("trait_A", "trait_B"), each = 3),
    p_val = rep(0.001, 6),
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2, -2.5)
  )

  result <- suppressMessages(
    enrich_gc_ora_go(
      dea_res,
      orgdb = "org.Hs.eg.db",
      ont = "MF",
      p_cutoff = 0.05
    )
  )

  expect_s3_class(
    result,
    c("glyfun_gc_ora_go_res", "glyfun_gc_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "trait",
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

test_that("enrich_gc_ora_go returns correct structure on happy path (integration)", {
  dea_res <- .mock_dea_res()

  result <- suppressMessages(
    enrich_gc_ora_go(
      dea_res,
      orgdb = "org.Hs.eg.db",
      ont = "MF",
      p_cutoff = 0.05
    )
  )

  expect_s3_class(
    result,
    c("glyfun_gc_ora_go_res", "glyfun_gc_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "trait",
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

test_that("enrich_gc_ora_kegg returns correct structure on happy path (integration)", {
  skip_if_offline()

  dea_res <- .mock_dea_res()

  suppressMessages(
    result <- enrich_gc_ora_kegg(dea_res, organism = "hsa", p_cutoff = 0.05)
  )

  expect_s3_class(
    result,
    c("glyfun_ora_kegg_res", "glyfun_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "trait",
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

test_that("enrich_gc_ora_reactome returns correct structure on happy path (integration)", {
  skip_if_offline()
  skip_if_not_installed("ReactomePA")
  skip_if_not_installed("org.Hs.eg.db")

  dea_res <- .mock_dea_res()

  suppressMessages(
    result <- enrich_gc_ora_reactome(
      dea_res,
      organism = "human",
      p_cutoff = 0.05
    )
  )

  expect_s3_class(
    result,
    c("glyfun_gc_ora_reactome_res", "glyfun_gc_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "trait",
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

test_that("enrich_gc_ora_wp returns correct structure on happy path (integration)", {
  skip_if_offline()
  skip_if_not_installed("org.Hs.eg.db")

  dea_res <- .mock_dea_res()

  suppressMessages(
    result <- enrich_gc_ora_wp(
      dea_res,
      organism = "Homo sapiens",
      p_cutoff = 0.05
    )
  )

  expect_s3_class(
    result,
    c("glyfun_gc_ora_wp_res", "glyfun_gc_ora_res", "glyfun_res")
  )
  expect_named(result, c("tidy_result", "raw_result"))
  expect_true(tibble::is_tibble(result$tidy_result))

  expected_cols <- c(
    "trait",
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

test_that("enrich_gc_ora_go errors on data.frame with missing columns", {
  # Missing 'trait' column
  dea_res_missing_trait <- tibble::tibble(
    protein = c("P01308", "P04637", "P42345", "P00533", "P42336"),
    site = rep("site1", 5),
    p_val = rep(0.001, 5),
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
    log2FC = c(2.5, 3.0, 1.5, 2.0, 2.2)
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
  # Use fake proteins that won't match any GO terms
  tidy_result <- tibble::tibble(
    variable = paste0("var", 1:3),
    trait = rep("trait_A", 3),
    protein = c("FAKE1", "FAKE2", "FAKE3"),
    p_val = rep(0.001, 3),
    p_adj = rep(0.001, 3),
    log2fc = c(2.5, 3.0, 1.5),
    estimate = rep(1, 3)
  )
  dea_res <- structure(
    list(tidy_result = tidy_result, raw_result = list()),
    class = c("glystats_ttest_res", "glystats_res")
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
