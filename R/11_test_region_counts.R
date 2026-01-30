#' Test region counts in a chromcall sample
#'
#' @description
#' Conduct Poisson-based statistical tests on read counts across genomic regions
#' in a chromcall [RangedSummarizedExperiment] object. Adds assays for logFC,
#' background lambda, p-values, adjusted p-values, enrichment score, and a
#' Poisson Z-score.
#'
#' In control mode, region-specific background weighting is estimated from the control
#' track (`modification_factor`). In no-control mode, background weighting is set to 1
#' everywhere (i.e. no input-based correction).
#'
#' @param x A [RangedSummarizedExperiment] object containing counts and metadata.
#' @param trim_blacklist Logical; if `TRUE`, removes rows where `blacklist == TRUE` before testing.
#'
#' @return A modified [RangedSummarizedExperiment] object with additional assays:
#' `logFC`, `lambda_t`, `p_value`, `p_adj`, `score`, `z_pois`,
#' and a rowData column `modification_factor`.
#'
#' @examples
#' ## With control (ChIP/CUT&RUN-style)
#' sample <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall"),
#'     H3K4me3  = system.file("extdata", "h3k4me3_sampleA.bam", package = "chromcall"),
#'     Control  = system.file("extdata", "control_sampleA.bam", package = "chromcall")
#'   ),
#'   control_name    = "Control",
#'   genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
#'   region_file     = system.file("extdata", "example.bed", package = "chromcall"),
#'   window_size     = 10000,
#'   blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall"),
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#' )
#'
#' # Run region-level testing
#' result <- test_region_counts(sample)
#'
#' # Explore outputs
#' result
#' SummarizedExperiment::assayNames(result)
#' SummarizedExperiment::colData(result)
#' head(SummarizedExperiment::assay(result, "counts"))
#' head(SummarizedExperiment::assay(result, "lambda_t"))
#' head(SummarizedExperiment::assay(result, "p_adj"))
#' head(SummarizedExperiment::assay(result, "score"))
#' head(SummarizedExperiment::assay(result, "z_pois"))
#' head(SummarizedExperiment::rowData(result)$modification_factor)
#' head(SummarizedExperiment::rowData(result)$expression)
#'
#' ## No-control mode (ATAC-seq-style)
#' # Note: we use a package example BAM here for demonstration only.
#' atac_sample <- build_chromcall_sample(
#'   sample_name     = "ATAC_A",
#'   experiments     = list(ATAC = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall")),
#'   control_name    = NULL,
#'   genome_file     = system.file("extdata", "genome.txt", package = "chromcall"),
#'   region_file     = system.file("extdata", "example.bed", package = "chromcall"),
#'   window_size     = 10000,
#'   blacklist_file  = system.file("extdata", "blacklist.bed", package = "chromcall"),
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#' )
#'
#' atac_res <- test_region_counts(atac_sample)
#'
#' # Explore outputs (no-control)
#' atac_res
#' SummarizedExperiment::assayNames(atac_res)
#' SummarizedExperiment::colData(atac_res)
#' head(SummarizedExperiment::assay(atac_res, "counts"))
#' head(SummarizedExperiment::assay(atac_res, "lambda_t"))
#' head(SummarizedExperiment::assay(atac_res, "p_adj"))
#' head(SummarizedExperiment::assay(atac_res, "score"))
#' head(SummarizedExperiment::assay(atac_res, "z_pois"))
#'
#' # Note: in no-control mode, logFC is NA because there is no matched control track.
#' head(SummarizedExperiment::assay(atac_res, "logFC"))
#'
#' @export
test_region_counts <- function(x, trim_blacklist = TRUE) {
  if (isTRUE(trim_blacklist) && "blacklist" %in% colnames(SummarizedExperiment::rowData(x))) {
    keep <- !SummarizedExperiment::rowData(x)$blacklist
    keep[is.na(keep)] <- TRUE
    x <- x[keep]
    SummarizedExperiment::rowData(x)$blacklist <- NULL
  }

  pseudo <- 0.01
  eps    <- 1e-8

  counts <- SummarizedExperiment::assays(x)$counts
  rn     <- rownames(counts)
  cn     <- colnames(counts)

  ctrl <- S4Vectors::metadata(x)$control_name
  no_control <- isTRUE(S4Vectors::metadata(x)$no_control) ||
    is.null(ctrl) || is.na(ctrl) || !nzchar(ctrl) || !(ctrl %in% cn)

  # ----------------------------
  # 1) modification_factor (weighting)
  # ----------------------------
  if (!no_control) {
    lambda_ctrl <- as.numeric(SummarizedExperiment::colData(x)[ctrl, "lambda_g"])
    if (!is.finite(lambda_ctrl) || lambda_ctrl <= 0) {
      stop("Control lambda_g must be positive and finite.")
    }
    modfac <- (counts[, ctrl] / lambda_ctrl)
    modfac[modfac < 1] <- 1
  } else {
    modfac <- rep(1, nrow(counts))
  }
  SummarizedExperiment::rowData(x)$modification_factor <- modfac

  # ----------------------------
  # 2) expected lambda per experiment
  #    Control:    E_ij = modfac_i * lambda_gj
  #    No control: E_ij = 1       * lambda_gj
  # ----------------------------
  lambda_g <- as.numeric(SummarizedExperiment::colData(x)$lambda_g)
  testlambda_matrix <- do.call(cbind, lapply(lambda_g, function(l) modfac * l))
  dimnames(testlambda_matrix) <- list(rn, cn)

  # Avoid zeros in denominator / sqrt
  testlambda_safe <- pmax(testlambda_matrix, eps)

  # ----------------------------
  # 3) logFC
  # ----------------------------
  if (!no_control) {
    control_counts <- counts[, ctrl]
    logfc_matrix <- log2((counts + pseudo) / (control_counts + pseudo))
  } else {
    # No-control mode: logFC is undefined (no matched control track)
    logfc_matrix <- matrix(
      NA_real_,
      nrow = nrow(counts),
      ncol = ncol(counts),
      dimnames = list(rn, cn)
    )
  }

  # ----------------------------
  # 4) one-sided Poisson p-values
  # ----------------------------
  pval_matrix <- do.call(cbind, lapply(seq_len(ncol(counts)), function(i) {
    obs <- counts[, i]
    exp <- testlambda_safe[, i]

    vapply(seq_along(obs), function(j) {
      stats::poisson.test(obs[j], exp[j], alternative = "greater")$p.value
    }, numeric(1))
  }))
  dimnames(pval_matrix) <- list(rn, cn)

  padj_matrix <- apply(pval_matrix, 2, stats::p.adjust, method = "fdr")
  if (is.null(dimnames(padj_matrix))) dimnames(padj_matrix) <- list(rn, cn)

  # ----------------------------
  # 5) effect size + z-score
  # ----------------------------
  enrichment_matrix <- log2((counts + pseudo) / (testlambda_safe + pseudo))
  dimnames(enrichment_matrix) <- list(rn, cn)

  z_pois_matrix <- (counts - testlambda_safe) / sqrt(testlambda_safe)
  dimnames(z_pois_matrix) <- list(rn, cn)

  # ----------------------------
  # 6) write outputs
  # ----------------------------
  SummarizedExperiment::assays(x)$logFC    <- logfc_matrix
  SummarizedExperiment::assays(x)$lambda_t <- testlambda_matrix
  SummarizedExperiment::assays(x)$p_value  <- pval_matrix
  SummarizedExperiment::assays(x)$p_adj    <- padj_matrix
  SummarizedExperiment::assays(x)$score    <- enrichment_matrix
  SummarizedExperiment::assays(x)$z_pois   <- z_pois_matrix

  return(x)
}
