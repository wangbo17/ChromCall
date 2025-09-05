#' Test region counts in a chromcall sample
#'
#' @description
#' Conduct Poisson-based statistical tests on read counts across genomic regions
#' in a chromcall [RangedSummarizedExperiment] object. Adds assays for logFC,
#' background lambda, p-values, adjusted p-values, enrichment score, and a
#' Poisson Z-score.
#'
#' @param x A [RangedSummarizedExperiment] object containing counts and metadata.
#' @param trim_blacklist Logical; if `TRUE`, removes rows where `blacklist == TRUE` before testing.
#'
#' @return A modified [RangedSummarizedExperiment] object with additional assays:
#' `logFC`, `lambda_t`, `p_value`, `p_adj`, `score`, `z_pois`,
#' and a rowData column `modification_factor`.
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
#' head(SummarizedExperiment::assay(result, "z_pois"))
#' head(SummarizedExperiment::rowData(result)$expression)
#' SummarizedExperiment::colData(result)          # Experiment metadata (e.g. lambda_g)
#' head(SummarizedExperiment::rowData(result))    # Region-level info (e.g. modification_factor)
#'
#' @export
test_region_counts <- function(x, trim_blacklist = TRUE) {
  if (isTRUE(trim_blacklist) && "blacklist" %in% colnames(SummarizedExperiment::rowData(x))) {
    keep <- !SummarizedExperiment::rowData(x)$blacklist
    keep[is.na(keep)] <- TRUE
    x <- x[keep]
    SummarizedExperiment::rowData(x)$blacklist <- NULL
  }

  pseudo  <- 0.01
  counts  <- SummarizedExperiment::assays(x)$counts
  rn      <- rownames(counts)
  cn      <- colnames(counts)

  control_col <- S4Vectors::metadata(x)$control_name
  if (!(control_col %in% cn)) {
    stop("control_name '", control_col, "' not found in colnames(x).")
  }
  control_counts <- counts[, control_col]

  # logFC vs control
  logfc_matrix <- log2((counts + pseudo) / (control_counts + pseudo))
  dimnames(logfc_matrix) <- list(rn, cn)

  # region-specific modulation factor from control
  lambda_ctrl <- as.numeric(SummarizedExperiment::colData(x)[control_col, "lambda_g"])
  if (!is.finite(lambda_ctrl) || lambda_ctrl <= 0) {
    stop("Control lambda_g must be positive and finite.")
  }
  modfac <- control_counts / lambda_ctrl
  modfac[modfac < 1] <- 1
  SummarizedExperiment::rowData(x)$modification_factor <- modfac

  # expected lambda per experiment j: E_ij = modfac_i * lambda_gj
  lambda_g <- as.numeric(SummarizedExperiment::colData(x)$lambda_g)
  testlambda_matrix <- do.call(cbind, lapply(lambda_g, function(l) modfac * l))
  dimnames(testlambda_matrix) <- list(rn, cn)

  # exact one-sided Poisson test p-values
  pval_matrix <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) {
    obs <- counts[, i]
    exp <- testlambda_matrix[, i]
    vapply(seq_along(obs), function(j) {
      stats::poisson.test(obs[j], exp[j], alternative = "greater")$p.value
    }, numeric(1))
  }))
  dimnames(pval_matrix) <- list(rn, cn)

  padj_matrix <- apply(pval_matrix, 2, stats::p.adjust, method = "fdr")
  if (is.null(dimnames(padj_matrix))) dimnames(padj_matrix) <- list(rn, cn)

  # effect size: log2(O/E) with small pseudo
  enrichment_matrix <- log2((counts + pseudo) / (testlambda_matrix + pseudo))
  dimnames(enrichment_matrix) <- list(rn, cn)

  # Poisson Z-score
  z_pois_matrix <- (counts - testlambda_matrix) / sqrt(testlambda_matrix)
  dimnames(z_pois_matrix) <- list(rn, cn)

  # write outputs
  SummarizedExperiment::assays(x)$logFC    <- logfc_matrix
  SummarizedExperiment::assays(x)$lambda_t <- testlambda_matrix
  SummarizedExperiment::assays(x)$p_value  <- pval_matrix
  SummarizedExperiment::assays(x)$p_adj    <- padj_matrix
  SummarizedExperiment::assays(x)$score    <- enrichment_matrix
  SummarizedExperiment::assays(x)$z_pois   <- z_pois_matrix

  return(x)
}
