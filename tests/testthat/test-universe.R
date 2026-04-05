# Helper to create a minimal glyexp_experiment for error tests
.mock_glyexp_experiment <- function(
  protein_col = TRUE,
  exp_type = "glycoproteomics"
) {
  expr_mat <- matrix(runif(10), nrow = 5, ncol = 2)
  colnames(expr_mat) <- c("S1", "S2")
  rownames(expr_mat) <- c("V1", "V2", "V3", "V4", "V5")

  var_info <- tibble::tibble(
    variable = c("V1", "V2", "V3", "V4", "V5"),
    glycan_composition = glyrepr::glycan_composition(
      c(Hex = 1, HexNAc = 1),
      c(Hex = 1, HexNAc = 1),
      c(Hex = 1, HexNAc = 1),
      c(Hex = 1, HexNAc = 1),
      c(Hex = 1, HexNAc = 1)
    ),
    protein_site = c(1L, 2L, 3L, 4L, 5L)
  )
  if (protein_col) {
    var_info$protein <- c("P01308", "P04637", "P42345", "P00533", "P42336")
  }

  sample_info <- tibble::tibble(sample = c("S1", "S2"))
  meta_data <- list(exp_type = exp_type, glycan_type = "N")

  structure(
    list(
      expr_mat = expr_mat,
      var_info = var_info,
      sample_info = sample_info,
      meta_data = meta_data
    ),
    class = "glyexp_experiment"
  )
}


.mock_glystats_res <- function(
  protein_col = TRUE,
  exp_type = "glycoproteomics"
) {
  tidy_result <- tibble::tibble(
    variable = c("var1", "var2", "var3", "var4", "var5"),
    p_val = rep(0.001, 5),
    p_adj = rep(0.001, 5),
    log2fc = c(2.5, 3.0, 1.5, 2.0, 2.2),
    estimate = rep(1, 5)
  )
  if (protein_col) {
    tidy_result$protein <- c("P01308", "P04637", "P42345", "P00533", "P42336")
  }

  raw_result <- list()
  meta_data <- list(exp_type = exp_type, glycan_type = "N")

  structure(
    list(
      tidy_result = tidy_result,
      raw_result = raw_result,
      meta_data = meta_data
    ),
    class = c("glystats_ttest_res", "glystats_res")
  )
}


# Happy path tests ----

test_that("detected_universe returns proteins from glyexp_experiment", {
  exp <- .mock_glyexp_experiment()
  result <- detected_universe(exp)
  expect_equal(result, c("P01308", "P04637", "P42345", "P00533", "P42336"))
})


test_that("detected_universe returns proteins from glystats_res", {
  res <- .mock_glystats_res()
  result <- detected_universe(res)
  expect_equal(result, c("P01308", "P04637", "P42345", "P00533", "P42336"))
})


test_that("detected_universe works with traitproteomics glystats_res", {
  res <- .mock_glystats_res(exp_type = "traitproteomics")
  result <- detected_universe(res)
  expect_equal(result, c("P01308", "P04637", "P42345", "P00533", "P42336"))
})


# Error handling tests ----

test_that("detected_universe errors on glyexp_experiment without protein column", {
  exp <- .mock_glyexp_experiment(protein_col = FALSE)
  expect_error(
    detected_universe(exp),
    "There must be a"
  )
})


test_that("detected_universe errors on glystats_res without protein column", {
  res <- .mock_glystats_res(protein_col = FALSE)
  expect_error(
    detected_universe(res),
    "There must be a"
  )
})


test_that("detected_universe errors on experiment with invalid exp_type", {
  exp <- .mock_glyexp_experiment(exp_type = "glycomics")
  expect_error(
    detected_universe(exp),
    "Experiment type must be"
  )
})


test_that("detected_universe errors on glystats_res with invalid exp_type", {
  res <- .mock_glystats_res(exp_type = "glycomics")
  expect_error(
    detected_universe(res),
    "Experiment type must be"
  )
})


test_that("detected_universe errors on experiment missing exp_type in meta_data", {
  exp <- .mock_glyexp_experiment()
  exp$meta_data$exp_type <- NULL
  expect_error(
    detected_universe(exp),
    "There is no"
  )
})


test_that("detected_universe errors on glystats_res with NULL meta_data", {
  res <- .mock_glystats_res()
  res$meta_data <- NULL
  expect_error(
    detected_universe(res),
    "Cannot decide"
  )
})


test_that("detected_universe errors on glystats_res missing exp_type in meta_data", {
  res <- .mock_glystats_res()
  res$meta_data$exp_type <- NULL
  expect_error(
    detected_universe(res),
    "There is no"
  )
})
