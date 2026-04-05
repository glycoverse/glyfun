# TODO: Add a `check_parameters()` function.

#' Check argument `dea_res`
#' @noRd
.check_dea_res <- function(dea_res) {
  UseMethod(".check_dea_res")
}

.check_dea_res.glystats_res <- function(dea_res) {
  # Check type
  basic_class <- class(dea_res)[[1]]
  supported_classes <- c(
    "glystats_limma_res",
    "glystats_ttest_res",
    "glystats_wilcox_res"
  )
  if (!basic_class %in% supported_classes) {
    cli::cli_abort(c(
      "Unsupported input class for {.arg dea_res}.",
      "i" = "Expected: {.cls {supported_classes}}",
      "x" = "Got: {.cls {basic_class}}"
    ))
  }

  # Check meta-data
  exp_type <- dea_res$meta_data$exp_type
  if (!is.null(exp_type)) {
    if (!exp_type %in% c("glycoproteomics", "traitproteomics")) {
      cli::cli_abort(c(
        "The experiment for glystats DEA functions must be of {.val glycoproteomics} or {.val traitproteomics} type.",
        "x" = "Got {.val {exp_type}}"
      ))
    }
  }

  # Check columns
  if (!"protein" %in% colnames(dea_res$tidy_res)) {
    cli::cli_abort(c(
      "A {.field protein} column must be in {.arg dea_res}",
      "i" = "Did you mistakenly set {.arg add_info} to `FALSE` when calling the {.pkg glystats} function?"
    ))
  }

  # For glystats_limma_res, check number of groups.
  if (inherits(dea_res, "glystats_limma_res")) {
    tidy_res <- glystats::get_tidy_result(dea_res)
    contrasts <- paste0(tidy_res$ref_group, "-", tidy_res$test_group)
    if (length(unique(contrasts)) > 1) {
      cli::cli_abort(
        "{.pkg glyfun} functions does not support multi-group {.fn glystats::gly_limma} results."
      )
    }
  }
}

.check_dea_res.data.frame <- function(dea_res) {
  expected_cols <- c("protein", "site", "trait", "p_val", "log2FC")
  missing_cols <- setdiff(expected_cols, colnames(dea_res))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg dea_res} must have all expected columns.",
      "i" = "Expected columns: {.field {expected_cols}}",
      "x" = "Missing columns: {.field {missing_cols}}"
    ))
  }
}

.check_p_cutoff_arg <- function(p_cutoff) {
  checkmate::assert_number(p_cutoff, lower = 0, upper = 1)
}

.check_log2fc_cutoff_arg <- function(log2fc_cutoff) {
  checkmate::assert_numeric(log2fc_cutoff)
  if (length(log2fc_cutoff) != 2) {
    cli::cli_abort(
      "{.arg log2fc_cutoff} must have exactly 2 elements.",
      "x" = "Got {.val {length(log2fc_cutoff)}}"
    )
  }
  if (log2fc_cutoff[[1]] > 0) {
    cli::cli_abort(
      "The first element of {.arg log2fc_cutoff} must be 0 or negative."
    )
  }
  if (log2fc_cutoff[[2]] < 0) {
    cli::cli_abort(
      "The second element of {.arg log2fc_cutoff} must be 0 or positive."
    )
  }
}

#' Convert UniProt IDs to Entrez IDs
#' @param uniprot The UniProt IDs.
#' @param orgdb An OrgDb object.
#' @returns The transformed Entrez IDs.
#' @noRd
.uniprot_to_entrez <- function(uniprot, orgdb) {
  suppressWarnings(suppressMessages(
    entrez_ids <- clusterProfiler::bitr(
      uniprot,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = orgdb
    )$ENTREZID
  ))
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  n_failed <- length(uniprot) - length(entrez_ids)
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(uniprot) * 100, 1)
    cli::cli_alert_warning(
      "{.val {n_failed}} of {.val {length(uniprot)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs."
    )
  }
  entrez_ids
}
