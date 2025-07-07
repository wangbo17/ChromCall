#' Test region counts in a chromcall sample
#'
#' @description
#' Conduct Poisson-based statistical tests on read counts across genomic regions
#' in a chromcall [RangedSummarizedExperiment] object. Adds assays for logFC,
#' background lambda, p-values, adjusted p-values, and enrichment score.
#'
#' @param x A [RangedSummarizedExperiment] object containing counts and metadata.
#' @param trim_blacklist Logical; if `TRUE`, removes rows where `blacklist == TRUE` before testing.
#'
#' @return A modified [RangedSummarizedExperiment] object with additional assays:
#' `logFC`, `lambda_t`, `p_value`, `p_adj`, `score`, and a rowData column `modification_factor`.
#'
#' @examples
#' # See ?build_chromcall_sample for how to construct a sample
#' sample <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall"),
#'     H3K4me3 = system.file("extdata", "h3k4me3_sampleA.bam", package = "chromcall"),
#'     Control = system.file("extdata", "control_sampleA.bam", package = "chromcall")
#'   ),
#'   control_name = "Control",
#'   genome_file = system.file("extdata", "genome.txt", package = "chromcall"),
#'   region_file = system.file("extdata", "example.bed", package = "chromcall"),
#'   window_size = 10000,
#'   blacklist_file = system.file("extdata", "blacklist.bed", package = "chromcall"),
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#' )
#' result <- test_region_counts(sample)
#'
#' # Explore output
#' SummarizedExperiment::assayNames(result)       # Assays: counts, logFC, lambda_t, p_value, p_adj
#' head(SummarizedExperiment::assay(result, "counts"))
#' head(SummarizedExperiment::assay(result, "lambda_t"))
#' head(SummarizedExperiment::assay(result, "p_adj"))
#' head(SummarizedExperiment::assay(result, "score"))
#' head(SummarizedExperiment::rowData(result)$expression)
#' SummarizedExperiment::colData(result)          # Experiment metadata (e.g. lambda_g)
#' head(SummarizedExperiment::rowData(result))    # Region-level info (e.g. modification_factor)
#'
#' @export
test_region_counts <- function(x, trim_blacklist = TRUE) {
  if (identical(trim_blacklist, TRUE)) {
    x <- x[!SummarizedExperiment::rowData(x)$blacklist]
    SummarizedExperiment::rowData(x)$blacklist <- NULL
  }

  pseudo <- 0.01
  counts <- SummarizedExperiment::assays(x)$counts
  control_col <- S4Vectors::metadata(x)$control_name
  control_counts <- counts[, control_col]

  # Log fold change: sample vs control
  logfc_matrix <- log2((counts + pseudo) / (control_counts + pseudo))

  # Modification factor: control sample observed / expected (global lambda)
  modfac <- counts[, control_col] / SummarizedExperiment::colData(x)[control_col, "lambda_g"]
  modfac[modfac < 1] <- 1
  SummarizedExperiment::rowData(x)$modification_factor <- modfac

  # Adjusted expected lambda per sample
  lambda_g <- SummarizedExperiment::colData(x)$lambda_g
  testlambda_matrix <- do.call(cbind, lapply(lambda_g, function(l) modfac * l))

  # p-value matrix from Poisson test (observed vs expected)
  pval_matrix <- do.call(cbind, lapply(seq_len(ncol(x)), function(sample_i) {
    obs <- counts[, sample_i]
    exp <- testlambda_matrix[, sample_i]
    vapply(seq_along(obs), function(j) {
      stats::poisson.test(obs[j], exp[j], alternative = "greater")$p.value
    }, numeric(1))
  }))
  padj_matrix <- apply(pval_matrix, 2, stats::p.adjust, method = "fdr")

  # Log2 enrichment: observed vs adjusted expected
  enrichment_matrix <- log2((counts + pseudo) / (testlambda_matrix + pseudo))

  # Add all to assays
  SummarizedExperiment::assays(x, withDimnames = FALSE)$logFC         <- logfc_matrix
  SummarizedExperiment::assays(x, withDimnames = FALSE)$lambda_t      <- testlambda_matrix
  SummarizedExperiment::assays(x, withDimnames = FALSE)$p_value       <- pval_matrix
  SummarizedExperiment::assays(x, withDimnames = FALSE)$p_adj         <- padj_matrix
  SummarizedExperiment::assays(x, withDimnames = FALSE)$score <- enrichment_matrix

  return(x)
}
