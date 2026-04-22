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
  expected_cols <- c("protein", "site", "trait", "p_val", "log2fc")
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

#' Infer the glycan trait column from a tidy DEA result
#' @param tidy_dea_res A tidy DEA result.
#' @returns The column name representing glycan traits.
#' @noRd
.infer_gc_trait_col <- function(tidy_dea_res) {
  by <- dplyr::case_when(
    "glycan_structure" %in% colnames(tidy_dea_res) ~ "glycan_structure",
    "glycan_composition" %in% colnames(tidy_dea_res) ~ "glycan_composition",
    "trait" %in% colnames(tidy_dea_res) ~ "trait",
    "motif" %in% colnames(tidy_dea_res) ~ "motif"
  )
  if (is.na(by)) {
    required_cols <- c(
      "glycan_structure",
      "glycan_composition",
      "trait",
      "motif"
    )
    cli::cli_abort(c(
      "Cannot determine glycan traits.",
      "i" = "At least one of these columns is needed: {.field {required_cols}}",
      "i" = "Did you accidentally set {.arg add_info} to `FALSE` when performing DEA with {.pkg glystats}?"
    ))
  }
  by
}

#' Prepare the `orgdb` parameter
#'
#' If `orgdb` is an `OrgDb` object, return it directly.
#' If it is a string, check it the package is installed, then return the `OrgDb` object.
#'
#' @param orgdb An `OrgDb` object or a string.
#' @returns An `OrgDb` object.
#' @noRd
.prepare_orgdb <- function(orgdb) {
  if (inherits(orgdb, "OrgDb")) {
    return(orgdb)
  }

  # Deal with character input
  available <- c(
    "org.Ag.eg.db", # Anopheles
    "org.At.tair.db", # Arabidopsis
    "org.Bt.eg.db", # Bovine
    "org.Ce.eg.db", # Worm
    "org.Cf.eg.db", # Canine
    "org.Dm.eg.db", # Fly
    "org.Dr.eg.db", # Zebrafish
    "org.EcK12.eg.db", # E coli strain K12
    "org.EcSakai.eg.db", # E coli strain Sakai
    "org.Gg.eg.db", # Chicken
    "org.Hbacteriophora.eg.db", # Heterorhabditis bacteriophora
    "org.Hs.eg.db", # Human
    "org.Mm.eg.db", # Mouse
    "org.Mmu.eg.db", # Rhesus
    "org.Mxanthus.db", # Myxococcus xanthus DK 1622
    "org.Pf.plasmo.db", # Malaria
    "org.Pt.eg.db", # Chimp
    "org.Rn.eg.db", # Rat
    "org.Sc.sgd.db", # Yeast
    "org.Ss.eg.db", # Pig
    "org.Xl.eg.db" # Xenopus
  )
  if (!orgdb %in% available) {
    cli::cli_abort(c(
      "Unsupported OrgDb type: {.val {orgdb}}",
      "i" = "Available OrgDb: {.val {available}}"
    ))
  }
  rlang::check_installed(orgdb)
  getExportedValue(orgdb, orgdb)
}

#' Convert UniProt IDs to Entrez IDs
#' @param uniprot The UniProt IDs.
#' @param orgdb An OrgDb object.
#' @returns The transformed Entrez IDs.
#' @noRd
.uniprot_to_entrez <- function(uniprot, orgdb, drop_na = TRUE) {
  unique_uniprot <- unique(uniprot)
  suppressWarnings(
    suppressMessages(
      mapping <- clusterProfiler::bitr(
        unique_uniprot,
        fromType = "UNIPROT",
        toType = "ENTREZID",
        OrgDb = orgdb
      )
    )
  )
  lookup <- rlang::set_names(mapping$ENTREZID, mapping$UNIPROT)
  entrez_ids <- unname(lookup[uniprot])
  n_failed <- sum(is.na(entrez_ids))
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(uniprot) * 100, 1)
    cli::cli_alert_warning(
      "{.val {n_failed}} of {.val {length(uniprot)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs."
    )
  }
  if (drop_na) {
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  }
  entrez_ids
}

#' Convert UniProt IDs to Entrez IDs for a named list of proteins
#'
#' This is more efficient than calling `.uniprot_to_entrez()` on each list
#' element separately, as it performs the mapping in a single database query.
#'
#' @param pro_list A named list where values are character vectors of UniProt IDs.
#' @param orgdb An OrgDb object.
#' @returns A named list with the same structure as `pro_list`, containing Entrez IDs.
#' @noRd
.uniprot_to_entrez_prolist <- function(pro_list, orgdb) {
  # Get all unique proteins
  all_proteins <- unique(unlist(pro_list, use.names = FALSE))

  # Build mapping in a single bitr call
  mapping <- suppressWarnings(
    suppressMessages(
      clusterProfiler::bitr(
        all_proteins,
        fromType = "UNIPROT",
        toType = "ENTREZID",
        OrgDb = orgdb
      )
    )
  )

  # Report conversion failures
  n_failed <- length(all_proteins) - length(unique(mapping$UNIPROT))
  if (n_failed > 0) {
    pct_failed <- round(n_failed / length(all_proteins) * 100, 1)
    cli::cli_alert_warning(
      "{.val {n_failed}} of {.val {length(all_proteins)}} ({.val {pct_failed}}%) proteins failed to map to Entrez IDs."
    )
  }

  # Create lookup table: uniprot -> entrez
  lookup <- rlang::set_names(mapping$ENTREZID, mapping$UNIPROT)

  # Map each list element
  purrr::map(pro_list, function(proteins) {
    entrez_ids <- lookup[proteins]
    entrez_ids[!is.na(entrez_ids)]
  })
}

.reactome_orgdb <- function(organism) {
  checkmate::assert_choice(
    organism,
    c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
  )
  if (organism == "celegans") {
    cli::cli_abort("Celegans is not supported yet.")
  }
  orgdb_name <- dplyr::recode_values(
    organism,
    "human" ~ "org.Hs.eg.db",
    "rat" ~ "org.Rn.eg.db",
    "mouse" ~ "org.Mm.eg.db",
    "yeast" ~ "org.Sc.sgd.db",
    "zebrafish" ~ "org.Dr.eg.db",
    "fly" ~ "org.Dm.eg.db"
  )
  .prepare_orgdb(orgdb_name)
}

.wp_orgdb <- function(organism) {
  checkmate::assert_choice(organism, clusterProfiler::get_wp_organisms())
  # WikiPathways uses latin names, e.g. Homo sapiens.

  unsupported <- c(
    "Equus caballus",
    "Populus trichocarpa",
    "Solanum lycopersicum",
    "Zea mays"
  )
  if (organism %in% c(unsupported)) {
    cli::cli_abort("{organism} is not supported yet.")
  }

  orgdb_name <- dplyr::recode_values(
    organism,
    "Anopheles gambiae" ~ "org.Ag.eg.db",
    "Arabidopsis thaliana" ~ "org.At.tair.db",
    "Bos taurus" ~ "org.Bt.eg.db",
    "Caenorhabditis elegans" ~ "org.Ce.eg.db",
    "Canis familiaris" ~ "org.Cf.eg.db",
    "Drosophila melanogaster" ~ "org.Dm.eg.db",
    "Danio rerio" ~ "org.Dr.eg.db",
    "Gallus gallus" ~ "org.Gg.eg.db",
    "Homo sapiens" ~ "org.Hs.eg.db",
    "Mus musculus" ~ "org.Mm.eg.db",
    "Pan troglodytes" ~ "org.Pt.eg.db",
    "Rattus norvegicus" ~ "org.Rn.eg.db",
    "Saccharomyces cerevisiae" ~ "org.Sc.sgd.db",
    "Sus scrofa" ~ "org.Ss.eg.db"
  )
  .prepare_orgdb(orgdb_name)
}

.do_orgdb <- function(organism) {
  checkmate::assert_choice(organism, c("hsa", "mmu"))
  orgdb_name <- switch(organism, hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db")
  .prepare_orgdb(orgdb_name)
}

#' Check gene set size arguments
#' @noRd
.check_gs_size_args <- function(min_gs_size, max_gs_size) {
  checkmate::assert_int(min_gs_size, lower = 1)
  checkmate::assert_int(max_gs_size, lower = 1)
  if (min_gs_size > max_gs_size) {
    cli::cli_abort(
      "{.arg min_gs_size} must be less than or equal to {.arg max_gs_size}."
    )
  }
}
