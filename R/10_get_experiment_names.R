#' Return Experiment Names from a Chromcall Sample
#'
#' @description
#' Extract experiment names from a chromcall sample object produced by
#' [build_chromcall_sample()]. By default, the control experiment (if present)
#' is excluded.
#'
#' In *no-control mode* (i.e., `control_name` is `NULL`), this function returns
#' all experiment names.
#'
#' @param x A [SummarizedExperiment::RangedSummarizedExperiment] object,
#'   typically produced by [build_chromcall_sample()].
#' @param include_control Logical; if `TRUE`, include the control experiment in
#'   the result. Default is `FALSE`.
#'
#' @return A character vector of experiment names.
#'
#' @examples
#' ## With control (ChIP/CUT&RUN-style)
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' region_file <- system.file("extdata", "example.bed", package = "chromcall")
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#'
#' bam_dir <- system.file("extdata", package = "chromcall")
#' experiments <- list(
#'   H3K27me3 = file.path(bam_dir, "h3k27me3_sampleA.bam"),
#'   H3K4me3  = file.path(bam_dir, "h3k4me3_sampleA.bam"),
#'   Control  = file.path(bam_dir, "control_sampleA.bam")
#' )
#'
#' sample <- build_chromcall_sample(
#'   sample_name     = "sampleA",
#'   experiments     = experiments,
#'   control_name    = "Control",
#'   genome_file     = genome_file,
#'   region_file     = region_file,
#'   window_size     = 10000,
#'   blacklist_file  = blacklist_file,
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#' )
#'
#' # Explore sample
#' sample
#' SummarizedExperiment::colData(sample)
#'
#' # Get experiment names (excluding control)
#' get_experiment_names(sample)
#'
#' # Get all names including control
#' get_experiment_names(sample, include_control = TRUE)
#'
#' ## No-control mode (ATAC-seq-style)
#' # Note: we use a package example BAM here for demonstration only.
#' atac_sample <- build_chromcall_sample(
#'   sample_name     = "ATAC_A",
#'   experiments     = list(ATAC = file.path(bam_dir, "h3k27me3_sampleA.bam")),
#'   control_name    = NULL,
#'   genome_file     = genome_file,
#'   region_file     = region_file,
#'   window_size     = 10000,
#'   blacklist_file  = blacklist_file,
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#' )
#'
#' # Explore no-control sample
#' atac_sample
#' SummarizedExperiment::colData(atac_sample)
#'
#' # Get experiment names (no control present)
#' get_experiment_names(atac_sample)
#'
#' @export
get_experiment_names <- function(x, include_control = FALSE) {
  # Prefer the explicit experiment name column if present (more robust)
  cd <- SummarizedExperiment::colData(x)

  if ("name" %in% colnames(cd)) {
    nms <- as.character(cd$name)
    # fall back if name is somehow empty
    if (all(is.na(nms)) || all(!nzchar(nms))) {
      nms <- colnames(x)
    }
  } else {
    nms <- colnames(x)
  }

  ctrl <- S4Vectors::metadata(x)$control_name
  no_control <- isTRUE(S4Vectors::metadata(x)$no_control) ||
    is.null(ctrl) || is.na(ctrl) || !nzchar(ctrl)

  if (no_control) {
    return(nms)
  }

  if (!isTRUE(include_control)) {
    nms <- nms[nms != ctrl]
  }

  nms
}
