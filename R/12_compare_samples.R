#' Compare two chromcall samples
#'
#' @description
#' Compare the specified experiments in two chromcall sample [`RangedSummarizedExperiment`] objects.
#' This function extracts adjusted p-values, assigns significance based on a threshold, computes
#' enrichment scores (observed / expected counts), and calculates delta scores between samples.
#' If both samples include gene expression values (`expression` column in `rowData`), it also adds
#' `<sample>_expression` and `log2FC_expression` to the output.
#'
#' @param x A chromcall sample [`RangedSummarizedExperiment`] object
#' @param y A chromcall sample [`RangedSummarizedExperiment`] object
#' @param experiments Optional character vector specifying which experiments to compare.
#' If not specified, compares all shared experiments.
#' @param threshold The adjusted p-value significance threshold (default: 0.05)
#'
#' @return A [`GRanges`] object containing per-region comparison results, including:
#' - Per-sample adjusted p-values, class (significance), and enrichment scores (counts / expected)
#' - Delta scores for each experiment (score_y - score_x)
#' - Optional: expression values and `log2FC_expression` if available in both samples
#'
#' @examples
#' # Load two chromcall sample objects with expression
#' sampleA <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall"),
#'     H3K4me3  = system.file("extdata", "h3k4me3_sampleA.bam", package = "chromcall"),
#'     Control  = system.file("extdata", "control_sampleA.bam", package = "chromcall")
#'   ),
#'   control_name = "Control",
#'   genome_file = system.file("extdata", "genome.txt", package = "chromcall"),
#'   region_file = system.file("extdata", "example.bed", package = "chromcall"),
#'   window_size = 10000,
#'   blacklist_file = system.file("extdata", "blacklist.bed", package = "chromcall"),
#'   expression_file = system.file("extdata", "expression_sampleA.bed", package = "chromcall")
#' )
#'
#' sampleB <- build_chromcall_sample(
#'   sample_name = "sampleB",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleB.bam", package = "chromcall"),
#'     H3K4me3  = system.file("extdata", "h3k4me3_sampleB.bam", package = "chromcall"),
#'     Control  = system.file("extdata", "control_sampleB.bam", package = "chromcall")
#'   ),
#'   control_name = "Control",
#'   genome_file = system.file("extdata", "genome.txt", package = "chromcall"),
#'   region_file = system.file("extdata", "example.bed", package = "chromcall"),
#'   window_size = 10000,
#'   blacklist_file = system.file("extdata", "blacklist.bed", package = "chromcall"),
#'   expression_file = system.file("extdata", "expression_sampleB.bed", package = "chromcall")
#' )
#'
#' # Run region-level differential testing
#' resultA <- test_region_counts(sampleA)
#' resultB <- test_region_counts(sampleB)
#'
#' # Compare samples (returns a GRanges object)
#' comparison <- compare_samples(resultA, resultB, threshold = 0.05)
#' comparison
#'
#' @export
compare_samples <- function(x, y, experiments = NULL, threshold = 0.05) {
  stopifnot(
    "Range mismatch: x and y must have identical rowRanges" =
      all(SummarizedExperiment::rowRanges(x) == SummarizedExperiment::rowRanges(y))
  )

  if (is.null(experiments)) {
    experiments <- intersect(get_experiment_names(x), get_experiment_names(y))
  }

  for (exp in experiments) {
    if (!(exp %in% get_experiment_names(x) && exp %in% get_experiment_names(y))) {
      stop(sprintf("Experiment '%s' not found in both samples", exp))
    }
  }

  sample_names <- c(S4Vectors::metadata(x)$sample_name, S4Vectors::metadata(y)$sample_name)
  output <- SummarizedExperiment::rowRanges(x)

  if ("modification_factor" %in% names(GenomicRanges::mcols(output))) {
    GenomicRanges::mcols(output)$modification_factor <- NULL
  }

  if ("expression" %in% names(GenomicRanges::mcols(output))) {
    GenomicRanges::mcols(output)$expression <- NULL
  }

  for (i in 1:2) {
    s <- if (i == 1) x else y
    s_name <- sample_names[i]

    if ("expression" %in% names(SummarizedExperiment::rowData(s))) {
      GenomicRanges::mcols(output)[[paste0(s_name, "_expression")]] <-
        SummarizedExperiment::rowData(s)$expression
    }

    for (exp in experiments) {
      padj <- SummarizedExperiment::assay(s, "p_adj")[, exp]
      enrichment <- SummarizedExperiment::assay(s, "score")[, exp]

      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_padj")]] <- padj
      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_class")]] <- as.numeric(padj <= threshold)
      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_score")]] <- enrichment
    }
  }

  expr1_col <- paste0(sample_names[1], "_expression")
  expr2_col <- paste0(sample_names[2], "_expression")
  if (expr1_col %in% names(GenomicRanges::mcols(output)) &&
      expr2_col %in% names(GenomicRanges::mcols(output))) {
    log2fc <- log2((GenomicRanges::mcols(output)[[expr2_col]] + 1e-3) /
                     (GenomicRanges::mcols(output)[[expr1_col]] + 1e-3))
    GenomicRanges::mcols(output)[["log2FC_expression"]] <- log2fc
  }

  for (exp in experiments) {
    score_x <- GenomicRanges::mcols(output)[[paste0(sample_names[1], "_", exp, "_score")]]
    score_y <- GenomicRanges::mcols(output)[[paste0(sample_names[2], "_", exp, "_score")]]
    GenomicRanges::mcols(output)[[paste0(exp, "_DeltaScore")]] <- score_y - score_x
  }

  return(output)
}
