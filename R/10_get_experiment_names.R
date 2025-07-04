#' Return Experiment Names from a Chromcall Sample
#'
#' @description
#' Extracts the names of experimental samples from a chromcall sample object,
#' optionally excluding the control sample.
#'
#' @param x A [RangedSummarizedExperiment] object, typically produced by [build_chromcall_sample()].
#' @param include_control Logical; if `TRUE`, include the control sample in the result. Default is `FALSE`.
#'
#' @return A character vector of experiment names.
#'
#' @examples
#' # Load example files
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' region_file <- system.file("extdata", "example.bed", package = "chromcall")
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#'
#' # Define BAM files for a sample
#' bam_dir <- system.file("extdata", package = "chromcall")
#' experiments <- list(
#'   H3K27me3 = file.path(bam_dir, "h3k27me3_sampleA.bam"),
#'   H3K4me3  = file.path(bam_dir, "h3k4me3_sampleA.bam"),
#'   Control  = file.path(bam_dir, "control_sampleA.bam")
#' )
#'
#' # Build sample object
#' sample <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = experiments,
#'   control_name = "Control",
#'   genome_file = genome_file,
#'   region_file = region_file,
#'   window_size = 10000,
#'   blacklist_file = blacklist_file
#' )
#'
#' # Get experiment names (excluding control)
#' get_experiment_names(sample)
#'
#' # Get all names including control
#' get_experiment_names(sample, include_control = TRUE)
#'
#' @export
get_experiment_names <- function(x, include_control = FALSE) {
  names <- colnames(x)
  if (!isTRUE(include_control)) {
    names <- setdiff(names, S4Vectors::metadata(x)$control_name)
  }
  names
}
