hdo_available <- function() {
  if (
    !requireNamespace("DOSE", quietly = TRUE) ||
      !requireNamespace("GOSemSim", quietly = TRUE)
  ) {
    return(FALSE)
  }

  isTRUE(tryCatch(
    {
      suppressMessages(suppressWarnings(
        GOSemSim:::get_onto_data("HDO", table = "term", output = "data.frame")
      ))
      TRUE
    },
    error = function(e) FALSE
  ))
}

skip_if_no_hdo <- function() {
  if (!hdo_available()) {
    testthat::skip("HDO ontology DB unavailable in this environment.")
  }
}

is_r_devel <- function() {
  grepl("Under development", R.version.string, fixed = TRUE)
}
