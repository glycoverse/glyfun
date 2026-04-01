#' Glycan-Centric GO Over Representation Analysis
#'
#' @description
#' This function first groups the proteins in `dea_res` according to `by`,
#' then performs GO ORA analysis with `clusterProfiler::compareCluster()` and
#' `clusterProfiler::enrichGO()`
#'
#' @details
#' # What is glycan-centric enrichment?
#'
#' In traditional glycoproteomics data analysis,
#' we usually perform differential expression analysis (DEA) on glycoforms,
#' extract proteins that have dysregulated glycosylation,
#' then perform functional enrichment (e.g. GO) on these proteins.
#' This is what all the enrichment functions in glystats do (e.g. `glystats::gly_enrich_go()`).
#'
#' `glyfun` functions differ in that they link specific glycan traits with functional annotations.
#' Instead of answering the question
#' "Which functions are enriched in dysregulated glycoproteins?",
#' `glyfun` answers questions like
#' “Which functions are enriched in proteins with dysregulated core-fucosylation?”
#' Higher specificity, deeper insights. By focusing on distinct glycan motifs,
#' glyfun helps you pinpoint the functional relevance of specific glycosylation changes.
#'
#' # Common usage pattern
#'
#' A common pattern of using this function is:
#'
#' ```r
#' # 1. Use `glydet` to calculate derived traits or motif quantification.
#' trait_exp <- derive_traits(exp)  # or `quantify_motifs()`
#'
#' # 2. Perform differential analysis with `glystats`.
#' dea_res <- gly_ttest(trait_exp)
#'
#' # 3. Use this function.
#' go_res <- gc_ora_go(dea_res)
#' ```
#'
#' @param dea_res Differential analysis result. Can be one of:
#'   - Result from [glystats::gly_limma()] (two groups), [glystats::gly_ttest()], or [glystats::gly_wilcox()],
#'     called on an [glyexp::experiment()] of "traitproteomics" type.
#'   - A tibble with the following columns:
#'     - `protein`: Uniprot ID of proteins
#'     - `trait`: A glycosylation trait (e.g. "TFc" for proportion of core-fucosylated glycans)
#'     - `p_val`: p-values, preferably adjusted p-values
#'     - `log2FC`: log2 of fold change
#' @param by A column to group the proteins by.
#'   - If `dea_res` is a [glyexp::experiment()]: the column name in `var_info` of the experiment.
#'   - If `dea_res` is a tibble: the column name in the tibble (defaults to "trait").
#' @param p_cutoff P-value cutoff for statistical significance. Defaults to 0.05.
#'   For `glystats` result input, adjusted p-values are used.
#' @param log2fc_cutoff Log2 fold change cutoff statistical significance.
#'   A length-2 numeric vector, being negative and positive boundaries, respectively.
#'   For example, `c(-1, 1)` means "log2FC < -1 or log2FC > 1", and `c(-Inf, 1)` means "log2FC > 1".
#'   Defaults to `c(-1, 1)`.
#'
#' @return A list with two elements:
#'  - `tidy_result`: A tibble with enrichment results containing the following columns:
#'    - `trait`: Glycan trait
#'    - `id`: Term ID
#'    - `description`: Term description
#'    - `gene_ratio`: Ratio of genes in the term to total genes in the input
#'    - `bg_ratio`: Ratio of genes in the term to total genes in the background
#'    - `p_val`: Raw p-value from hypergeometric test
#'    - `p_adj`: Adjusted p-value
#'    - `q_value`: Q-value (FDR)
#'    - `gene_id`: Gene IDs in the term (separated by "/")
#'    - `count`: Number of genes in the term
#'  - `raw_result`: The raw clusterProfiler clusterProfResult object
#' The list has classes `glystats_go_ora_res` and `glystats_res`.
#'
#' @seealso [clusterProfiler::compareCluster()], [clusterProfiler::enrichGO()]
#' @export
gc_ora_go <- function(dea_res, by = NULL, p_cutoff = 0.05, log2fc_cutoff = c(-1, 1)) {
  UseMethod("gc_ora_go")
}