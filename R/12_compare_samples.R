#' Compare two chromcall samples
#'
#' @description
#' Compare the specified experiments in two chromcall sample [`RangedSummarizedExperiment`] objects.
#' This function extracts adjusted p-values, assigns significance based on a threshold, and
#' collects enrichment scores (log2(O/E), stored as `score`) and z-scores (`z_pois`).
#' It also calculates per-experiment delta metrics between samples.
#'
#' If both samples include gene expression values (`expression` column in `rowData`), the output
#' additionally includes `<sample>_expression` and `log2FC_expression`.
#'
#' For interpretation, it is recommended that both samples are processed under the same mode
#' (control vs no-control).
#'
#' Note: This function expects that both `x` and `y` have already been processed by
#' [test_region_counts()], so that assays `p_adj`, `score`, and `z_pois` exist.
#'
#' @param x A chromcall sample [`RangedSummarizedExperiment`] object (after [test_region_counts()]).
#' @param y A chromcall sample [`RangedSummarizedExperiment`] object (after [test_region_counts()]).
#' @param experiments Optional character vector specifying which experiments to compare.
#'   If not specified, compares all shared experiments (excluding control if present).
#' @param threshold The adjusted p-value significance threshold (default: 0.05).
#'
#' @return A [`GRanges`] object containing per-region comparison results, including:
#' - Per-sample adjusted p-values (`_padj`), significance class (`_class`),
#'   enrichment scores (`_enrichment_score`), and z-scores (`_z_score`)
#' - Delta metrics for each experiment: `<exp>_DeltaEnrichment` and `<exp>_DeltaZscore`
#' - Optional: expression values and `log2FC_expression` if available in both samples
#'
#' @examples
#' ## Build and test two samples (control-mode)
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' region_file <- system.file("extdata", "example.bed", package = "chromcall")
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#' expr_file <- system.file("extdata", "expression_tss.bed", package = "chromcall")
#'
#' sampleA <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall"),
#'     H3K4me3  = system.file("extdata", "h3k4me3_sampleA.bam", package = "chromcall"),
#'     Control  = system.file("extdata", "control_sampleA.bam", package = "chromcall")
#'   ),
#'   control_name    = "Control",
#'   genome_file     = genome_file,
#'   region_file     = region_file,
#'   window_size     = 10000,
#'   blacklist_file  = blacklist_file,
#'   expression_file = expr_file
#' )
#'
#' sampleB <- build_chromcall_sample(
#'   sample_name = "sampleB",
#'   experiments = list(
#'     H3K27me3 = system.file("extdata", "h3k27me3_sampleB.bam", package = "chromcall"),
#'     H3K4me3  = system.file("extdata", "h3k4me3_sampleB.bam", package = "chromcall"),
#'     Control  = system.file("extdata", "control_sampleB.bam", package = "chromcall")
#'   ),
#'   control_name    = "Control",
#'   genome_file     = genome_file,
#'   region_file     = region_file,
#'   window_size     = 10000,
#'   blacklist_file  = blacklist_file,
#'   expression_file = expr_file
#' )
#'
#' resultA <- test_region_counts(sampleA)
#' resultB <- test_region_counts(sampleB)
#'
#' comparison <- compare_samples(resultA, resultB, threshold = 0.05)
#' comparison
#'
#' @export
compare_samples <- function(x, y, experiments = NULL, threshold = 0.05) {

  # ---- sanity checks: required assays ----
  required_assays <- c("p_adj", "score", "z_pois")
  ax <- SummarizedExperiment::assayNames(x)
  ay <- SummarizedExperiment::assayNames(y)
  if (!all(required_assays %in% ax) || !all(required_assays %in% ay)) {
    stop(
      "compare_samples() requires both inputs to have assays: ",
      paste(required_assays, collapse = ", "),
      ". Please run test_region_counts() on both samples first."
    )
  }

  # ---- range match check (coordinate-level) ----
  rx <- SummarizedExperiment::rowRanges(x)
  ry <- SummarizedExperiment::rowRanges(y)

  same_ranges <- length(rx) == length(ry) &&
    all(as.character(GenomicRanges::seqnames(rx)) == as.character(GenomicRanges::seqnames(ry))) &&
    all(GenomicRanges::start(rx) == GenomicRanges::start(ry)) &&
    all(GenomicRanges::end(rx) == GenomicRanges::end(ry)) &&
    all(as.character(GenomicRanges::strand(rx)) == as.character(GenomicRanges::strand(ry)))

  stopifnot("Range mismatch: x and y must have identical genomic coordinates" = same_ranges)

  # ---- experiments to compare ----
  x_exp_default <- get_experiment_names(x)  # excludes control if present
  y_exp_default <- get_experiment_names(y)

  if (is.null(experiments)) {
    experiments <- intersect(x_exp_default, y_exp_default)
  } else {
    # if user explicitly provides experiments, validate against all colnames (allow Control if asked)
    if (!all(experiments %in% colnames(x)) || !all(experiments %in% colnames(y))) {
      stop(
        "Some experiments are not present in both samples. ",
        "x colnames: ", paste(colnames(x), collapse = ", "),
        "; y colnames: ", paste(colnames(y), collapse = ", ")
      )
    }
  }

  if (length(experiments) == 0) {
    stop("No shared experiments to compare between x and y.")
  }

  # ---- safe sample names (avoid NA/dup collisions) ----
  s1 <- S4Vectors::metadata(x)$sample_name
  s2 <- S4Vectors::metadata(y)$sample_name
  if (is.null(s1) || is.na(s1) || !nzchar(s1)) s1 <- "sample1"
  if (is.null(s2) || is.na(s2) || !nzchar(s2)) s2 <- "sample2"
  if (identical(s1, s2)) {
    s1 <- paste0(s1, "_1")
    s2 <- paste0(s2, "_2")
  }
  sample_names <- c(s1, s2)

  output <- rx

  # clean any existing columns to avoid conflicts
  if ("modification_factor" %in% names(GenomicRanges::mcols(output))) {
    GenomicRanges::mcols(output)$modification_factor <- NULL
  }
  if ("expression" %in% names(GenomicRanges::mcols(output))) {
    GenomicRanges::mcols(output)$expression <- NULL
  }

  # ---- optional expression columns ----
  for (i in 1:2) {
    s <- if (i == 1) x else y
    s_name <- sample_names[i]
    if ("expression" %in% names(SummarizedExperiment::rowData(s))) {
      GenomicRanges::mcols(output)[[paste0(s_name, "_expression")]] <-
        SummarizedExperiment::rowData(s)$expression
    }
  }

  # ---- per-sample per-experiment stats ----
  for (i in 1:2) {
    s <- if (i == 1) x else y
    s_name <- sample_names[i]

    for (exp in experiments) {
      padj   <- SummarizedExperiment::assay(s, "p_adj")[, exp]
      enrich <- SummarizedExperiment::assay(s, "score")[, exp]
      z_vals <- SummarizedExperiment::assay(s, "z_pois")[, exp]

      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_padj")]]   <- padj
      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_class")]]  <- as.numeric(padj <= threshold)
      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_enrichment_score")]] <- enrich
      GenomicRanges::mcols(output)[[paste0(s_name, "_", exp, "_z_score")]]          <- z_vals
    }
  }

  # ---- expression log2FC if both exist ----
  expr1_col <- paste0(sample_names[1], "_expression")
  expr2_col <- paste0(sample_names[2], "_expression")
  if (expr1_col %in% names(GenomicRanges::mcols(output)) &&
      expr2_col %in% names(GenomicRanges::mcols(output))) {
    log2fc <- log2((GenomicRanges::mcols(output)[[expr2_col]] + 1e-3) /
                     (GenomicRanges::mcols(output)[[expr1_col]] + 1e-3))
    GenomicRanges::mcols(output)[["log2FC_expression"]] <- log2fc
  }

  # ---- delta metrics ----
  for (exp in experiments) {
    es_x <- GenomicRanges::mcols(output)[[paste0(sample_names[1], "_", exp, "_enrichment_score")]]
    es_y <- GenomicRanges::mcols(output)[[paste0(sample_names[2], "_", exp, "_enrichment_score")]]
    GenomicRanges::mcols(output)[[paste0(exp, "_DeltaEnrichment")]] <- es_y - es_x

    z_x <- GenomicRanges::mcols(output)[[paste0(sample_names[1], "_", exp, "_z_score")]]
    z_y <- GenomicRanges::mcols(output)[[paste0(sample_names[2], "_", exp, "_z_score")]]
    GenomicRanges::mcols(output)[[paste0(exp, "_DeltaZscore")]] <- z_y - z_x
  }

  return(output)
}
