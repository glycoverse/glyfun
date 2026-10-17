test_that("automatic trait selection preserves candidate priority", {
  candidates <- c("glycan_structure", "glycan_composition", "trait", "motif")
  for (candidate in candidates) {
    df <- tibble::tibble(value = c("A", "B"))
    names(df) <- candidate
    expect_no_message(expect_identical(.infer_gc_trait_col(df), candidate))
  }

  df <- tibble::tibble(
    motif = "M",
    trait = "T",
    glycan_composition = "C",
    glycan_structure = "S"
  )
  expect_snapshot({
    .infer_gc_trait_col(df)
    .infer_gc_trait_col(df, prefer_trait = TRUE)
  })
  expect_identical(
    suppressMessages(.infer_gc_trait_col(df[c(
      "motif",
      "trait",
      "glycan_composition"
    )])),
    "glycan_composition"
  )
  expect_no_message(expect_identical(
    .infer_gc_trait_col(df, by = "motif"),
    "motif"
  ))
})

test_that("explicit trait selection accepts custom columns and factor labels", {
  df <- tibble::tibble(`custom trait` = factor(c("B", "A")), number = c(1, 2))
  expect_identical(.infer_gc_trait_col(df, by = "custom trait"), "custom trait")
  expect_identical(.infer_gc_trait_col(df, by = "number"), "number")
})

test_that("trait selection rejects invalid names and grouping columns clearly", {
  df <- tibble::tibble(trait = c("A", "B"))
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = "missing"))
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = character()))
  expect_snapshot(
    error = TRUE,
    .infer_gc_trait_col(df, by = c("trait", "motif"))
  )
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = NA_character_))
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = ""))
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = 1))
  expect_snapshot(
    error = TRUE,
    .infer_gc_trait_col(tibble::tibble(protein = "P1"))
  )

  df$trait <- list("A", "B")
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = "trait"))
  df$trait <- matrix(c("A", "B"), ncol = 1)
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = "trait"))
  df <- data.frame(trait = "A", other = "B")
  names(df) <- c("trait", "trait")
  expect_snapshot(error = TRUE, .infer_gc_trait_col(df, by = "trait"))
})
