#' Helper function to prepare the `universe` parameter
#'
#' This function extracts all detected proteins in a [glyexp::experiment()] or a `glystats` result.
#' It can be readily passed to the `universe` parameter of all `glyfun` functions.
#'
#' @param x A [glyexp::experiment()] or a `glystats` result.
#' @returns A character vector of protein UniProt IDs.
#' @examples
#' library(glyexp)
#' universe <- detected_universe(real_experiment)
#' length(universe)
#' universe[1:5]
#' @export
detected_universe <- function(x) {
  UseMethod("detected_universe")
}

#' @export
detected_universe.glyexp_experiment <- function(x) {
  .check_universe_meta_data(x$meta_data)
  proteins <- x$var_info$protein
  if (is.null(proteins)) {
    cli::cli_abort(c(
      "There must be a {.field protein} column in the {.field var_info} of the experiment.",
      "i" = "Did you accidentally remove the {.field protein} column?"
    ))
  }
  proteins
}

#' @export
detected_universe.glystats_res <- function(x) {
  .check_universe_meta_data(x$meta_data)
  proteins <- x$tidy_result$protein
  if (is.null(proteins)) {
    cli::cli_abort(c(
      "There must be a {.field protein} column in the {.field var_info} of the experiment.",
      "i" = "Did you accidentally set `add_info = FALSE` when calling the {.pkg glystats} function?"
    ))
  }
  proteins
}

#' Check the meta-data of input for `detected_universe()`
#' @param meta The meta data.
#' @param input_type "exp" for experiment, "res" for glystats result.
#' @noRd
.check_universe_meta_data <- function(meta, input_type = c("exp", "res")) {
  # If `meta`` is NULL, `input_type`` cannot be "exp",
  # because a `meta_data` field in enforced by an `experiment()` object.
  if (is.null(meta) && input_type == "res") {
    # The only reason for this case is that the user is using
    # a glystats version before 0.8.0.
    cli::cli_abort(c(
      "Cannot decide experiment type because a {.field meta_data} field is missing.",
      "i" = "You might need to update {.pkg glystats} to the lastest version.",
      "i" = "Current {.pkg glystats} version: {.val {.get_glystats_version()}}"
    ))
  }

  exp_type <- meta$exp_type
  if (is.null(exp_type)) {
    cli::cli_abort(c(
      "There is no {.field exp_type} in the meta data.",
      "i" = "This might be a bug. Please create a issue on GitHub."
    ))
  }

  if (!exp_type %in% c("glycoproteomics", "traitproteomics")) {
    cli::cli_abort(c(
      "Experiment type must be one of {.val glycoproteomics} or {.val traitproteomics}",
      "x" = "Got: {.val {exp_type}}"
    ))
  }
}

.get_glystats_version <- function() {
  tryCatch(
    utils::packageVersion("glystats"),
    error = function(e) NULL
  )
}