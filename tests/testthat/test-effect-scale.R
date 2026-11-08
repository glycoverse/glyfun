test_that("GSEA only requires quantities consumed by the selected ranking", {
  df <- data.frame(protein = c("P1", "P2"), p_val = c(.01, .1))
  expect_equal(.prepare_pro_list(df, "log10p", "median"), c(P1 = 2, P2 = 1))
  for (method in c("log2fc", "abs_log2fc", "signed_log10p", "log2fc_log10p")) {
    expect_error(.prepare_pro_list(df, method, "median"), "log2fc")
    df$log2fc <- NA_real_
    expect_error(.prepare_pro_list(df, method, "median"), "no finite")
    df$log2fc <- NULL
  }
  df$log2fc <- c(-2, 3)
  expect_equal(.calcu_rank_scores(df, "log2fc"), c(-2, 3))
  expect_equal(.calcu_rank_scores(df, "abs_log2fc"), c(2, 3))
  expect_equal(.calcu_rank_scores(df, "signed_log10p"), c(-2, 1))
  expect_equal(.calcu_rank_scores(df, "log2fc_log10p"), c(-4, 3))
  df$p_adj <- c(.1, .5)
  expect_equal(.calcu_rank_scores(df, "log10p"), -log10(df$p_adj))
})

test_that("invalid and empty enrichment quantities are explicit", {
  df <- data.frame(
    protein = c("P1", "P2"),
    log2fc = c(NA_real_, 2),
    p_val = c(.01, .1)
  )
  expect_warning(.check_ora_quantities(df), "non-finite")
  expect_warning(.prepare_pro_list(df, "log2fc", "median"), "non-finite")
  df$log2fc <- c(Inf, 2)
  expect_warning(
    expect_equal(
      .prepare_pro_list(df, "signed_log10p", "median"),
      c(P2 = 1)
    ),
    "non-finite ranking"
  )
  df$log2fc[] <- NA_real_
  expect_error(.check_ora_quantities(df), "no finite")
  df$p_val <- c(-.1, .1)
  expect_error(.calcu_rank_scores(df, "log10p"), "probabilities")
  df$p_val[] <- 0
  expect_error(.prepare_pro_list(df, "log10p", "median"), "no finite ranking")
  df$p_val[] <- NA_real_
  expect_error(.calcu_rank_scores(df, "log10p"), "no finite")
})

test_that("ORA does not send empty valid selections to enrichment backends", {
  df <- data.frame(
    protein = "P1",
    site = "s",
    trait = "t",
    p_val = .5,
    log2fc = .1
  )
  backend <- function(...) stop("backend must not run")
  expect_warning(expect_null(.ora(df, backend)), "No proteins pass")
  expect_warning(expect_null(.gc_ora(df, backend)), "No proteins pass")
})

test_that("current glystats results preserve scale contracts at enrichment entry", {
  skip_if_not_installed("glystats", "0.14.0")
  skip_if_not_installed("SummarizedExperiment")
  values <- rbind(c(1, 2, 3, 7, 8, 9), c(8, 9, 10, 2, 3, 4))
  dimnames(values) <- list(c("v1", "v2"), paste0("s", 1:6))
  for (scale in c("raw", "log", "log2", "log10", "logx")) {
    x <- switch(
      scale,
      raw = values,
      log = log(values),
      log2 = log2(values),
      log10 = log10(values),
      logx = log2(values) - 2
    )
    exp <- SummarizedExperiment::SummarizedExperiment(
      assays = list(abundance = x),
      rowData = S4Vectors::DataFrame(
        protein = c("P1", "P2"),
        site = c("s", "s"),
        trait = c("t", "t")
      ),
      colData = S4Vectors::DataFrame(group = rep(c("A", "B"), each = 3)),
      metadata = list(scale = scale, exp_type = "glycoproteomics")
    )
    res <- glystats::gly_ttest(exp)
    tidy <- glystats::get_tidy_result(res)
    expect_silent(.check_dea_res(res))
    expect_equal(.calcu_rank_scores(tidy, "log10p"), -log10(tidy$p_adj))
    if (scale == "raw") {
      limma_res <- glystats::gly_limma(exp, log_transform = FALSE)
      expect_error(.check_ora_quantities(limma_res), "no finite")
      expect_true(all(is.finite(.calcu_rank_scores(
        glystats::get_tidy_result(limma_res),
        "log10p"
      ))))
    }
    if (scale == "logx") {
      expect_true(all(is.na(tidy$log2fc)))
      expect_error(.check_ora_quantities(res), "no finite")
      expect_error(
        .prepare_pro_list(tidy, "signed_log10p", "median"),
        "no finite"
      )
    } else {
      expect_silent(.check_ora_quantities(res))
      expect_equal(.calcu_rank_scores(tidy, "log2fc"), tidy$log2fc)
    }
  }
})

test_that("GSEA entry points accept missing log2FC for significance ranking", {
  df <- data.frame(
    protein = c("P1", "P2"),
    site = "s",
    trait = "t",
    p_val = c(.01, .1)
  )
  captured <- NULL
  backend <- function(geneList, ...) {
    captured <<- geneList
    data.frame(ID = "term")
  }
  .gsea(df, backend, rank_by = "log10p", aggr = "median")
  expect_equal(captured, c(P1 = 2, P2 = 1))
  expect_equal(
    .prepare_gc_gsea_df(df, "trait", "log10p", "median")$score,
    c(2, 1)
  )
  df$log2fc <- NA_real_
  expect_error(.gc_ora(df, backend), "no finite")
})
