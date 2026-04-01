#' Process `by` argument for glystats input
#'
#' If `by` is NULL, return a valid value.
#' Otherwise, check if `by` exists.
#'
#' @param dea_res DEA result.
#' @param by `by` paramerter.
#'
#' @returns A valid `by` value.
#' @noRd
.process_by_arg_glystats <- function(dea_res, by) {
  checkmate::assert_string(by, null.ok = TRUE)

  tidy_res <- glystats::get_tidy_result(dea_res)

  # 1. `by` is NULL
  if (is.null(by)) {
    by <- dplyr::case_when(
      "glycan_structure" %in% colnames(tidy_res) ~ "glycan_structure",
      "glycan_composition" %in% colnames(tidy_res) ~ "glycan_composition",
      "trait" %in% colnames(tidy_res) ~ "trait",
      "motif" %in% colnames(tidy_res) ~ "motif",
      .unmatched = "error"
    )
    return(by)
  }

  # 2. `by` has a value
  available_cols <- .available_by_columns(tidy_res)
  if (!by %in% available_cols) {
    if (by %in% c("variable", "protein")) {
      cli::cli_abort(c(
        "{.arg by} cannot be {.val variable} or {.val protein}.",
        "i" = "Available columns: {.val {available_cols}}"
      ))
    }
    msg1 <- "{.arg by} must be a column in the {.field var_info} of the {.fn experiment} before calling the {.pkg glystats} function."
    msg2 <- "{.val {by}} is not in {.var dea_res}"
    if (length(available_cols) == 0) {
      cli::cli_abort(c(msg1, "x" = msg2, "i" = "Did you mistakenly set {.arg add_info} to `FALSE` when calling the {.pkg glystats} function?"))
    } else {
      cli::cli_abort(c(msg1, "x" = msg2, "i" = "Available columns: {.val {available_cols}}"))
    }
  }
  by
}

#' Process `by` argument for dataframe input
#'
#' If `by` is NULL, defaults to "trait".
#' Otherwise, it must be "trait".
#'
#' @param dea_res DEA result. A data.frame.
#' @param by `by` parameter.
#'
#' @returns A valid `by` value.
#' @noRd
.process_by_arg_df <- function(dea_res, by) {
  if (is.null(by)) {
    return("trait")
  } else {
    if (by != "trait") {
      cli::cli_abort("{.arg by} must be {.val trait} when {.arg dea_res} is a data.frame.")
    }
    return("trait")
  }
}

#' Guess available columns for `by`
#' @param res The tidy glystats result.
#' @returns A character vector of available `by` values.
#' @noRd
.available_by_columns <- function(res) {
  UseMethod(".available_by_columns")
}

.available_by_columns.glystats_limma_res <- function(res) {
  .available_by_columns_impl(res, "log2fc")
}

.available_by_columns.glystats_ttest_res <- function(res) {
  .available_by_columns_impl(res, "estimate")
}

.available_by_columns.glystats_wilcox_res <- function(res) {
  .available_by_columns_impl(res, "statistic")
}

.available_by_columns_impl <- function(res, col_before) {
  all_cols <- colnames(res)
  idx <- which(all_cols == col_before)
  cols <- all_cols[1:(idx - 1)]
  setdiff(cols, c("variable", "protein"))
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
    cli::cli_abort("The first element of {.arg log2fc_cutoff} must be 0 or negative.")
  }
  if (log2fc_cutoff[[2]] < 0) {
    cli::cli_abort("The second element of {.arg log2fc_cutoff} must be 0 or positive.")
  }
}

#' Rename parameters and call `clusterProfiler::compareCluster()`
#' @param ... Parameters to be passed.
#' @noRd
.call_compare_cluster <- function(proteins, ...) {
  dots <- rlang::list2(...)
  param_name_mapping <- c(
    "orgdb" = "OrgDb",
    "keytype" = "keyType",
    "p_adj_method" = "pAdjustMethod",
    "p_cutoff" = "pvalueCutoff",
    "q_cutoff" = "qvalueCutoff"
  )
  dots <- .rename_list(dots, param_name_mapping)
  rlang::exec(clusterProfiler::compareCluster, proteins, !!!dots)
}

#' Rename a list with a mapping
#'
#' If a name in `x` is in the names of `mapping`, it will be replaced.
#' Otherwise, the old names are used.
#'
#' @param x The list to be renamed.
#' @param mapping A character vector of name mapping, with the format `c("old_name" = "new_name")`.
#' @returns The renamed list.
#' @noRd
.rename_list <- function(x, mapping) {
  old_names <- names(x)
  new_names <- mapping[old_names]
  new_names <- ifelse(is.na(new_names), old_names, new_names)
  rlang::set_names(x, new_names)
}