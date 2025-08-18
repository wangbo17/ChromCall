#' Create a chromcall sample from input files
#'
#' @description
#' Constructs a single-sample [SummarizedExperiment] object by loading multiple ChIP/CUT&RUN experiments
#' with a shared genome, promoter windows, and blacklist information. Computes per-experiment background lambda.
#'
#' @param sample_name The name of the sample (string).
#' @param experiments A named list of BAM file paths for each experiment (e.g., H3K27me3, H3K4me3, control).
#' @param paired Logical. If TRUE (default), input BAM is assumed to be paired-end
#'   and reads are loaded as fragments. If FALSE, input BAM is treated as single-end.
#' @param control_name The name (matching `names(experiments)`) of the experiment to treat as control.
#' @param genome_file Path to the genome info text file.
#' @param region_file Path to the promoter BED file.
#' @param window_size Genome tile size in base pairs.
#' @param blacklist_file (Optional) Path to BED file of blacklist regions, or NULL to disable.
#' @param expression_file (Optional) Path to a BED-format TSS expression file to merge expression into regions.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] object containing region counts and lambda estimates.
#'
#' @examples
#' # Define input files
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' region_file <- system.file("extdata", "example.bed", package = "chromcall")
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#'
#' # Define BAM files for sample A
#' bam_dir <- system.file("extdata", package = "chromcall")
#' experiments <- list(
#'   H3K27me3 = file.path(bam_dir, "h3k27me3_sampleA.bam"),
#'   H3K4me3  = file.path(bam_dir, "h3k4me3_sampleA.bam"),
#'   Control  = file.path(bam_dir, "control_sampleA.bam")
#' )
#'
#' # Create sample object
#' sample <- build_chromcall_sample(
#'   sample_name = "sampleA",
#'   experiments = experiments,
#'   control_name = "Control",
#'   genome_file = genome_file,
#'   region_file = region_file,
#'   window_size = 10000,
#'   blacklist_file = blacklist_file,
#'   expression_file = system.file("extdata", "expression_tss.bed", package = "chromcall")
#'   )
#'
#' # Explore output
#' sample
#' SummarizedExperiment::assayNames(sample)
#' SummarizedExperiment::colData(sample)
#' head(SummarizedExperiment::assay(sample, "counts"))
#' head(SummarizedExperiment::rowData(sample)$expression)
#'
#' @export
build_chromcall_sample <- function(
    sample_name,
    experiments,
    paired = TRUE,
    control_name,
    genome_file,
    region_file,
    window_size,
    blacklist_file = NULL,
    expression_file = NULL
) {
  if (!control_name %in% names(experiments)) {
    stop("The specified control_name '", control_name, "' is not in names(experiments).")
  }

  genome <- load_genome(genome_file)

  blacklist <- if (!is.null(blacklist_file)) {
    load_bedfile(blacklist_file, genome = genome, reduce = TRUE)
  } else {
    NULL
  }

  genome_tiles <- tile_genome(genome, window_size = window_size, blacklist = blacklist)
  regions <- load_windows(region_file, genome = genome, blacklist = blacklist)

  if (!is.null(expression_file)) {
    exprs <- load_expression(expression_file, genome)

    overlaps <- GenomicRanges::findOverlaps(
      query = GenomicRanges::resize(regions, width = 1, fix = "center"),
      subject = exprs
    )

    S4Vectors::mcols(regions)$expression <- 0

    S4Vectors::mcols(regions)[
      S4Vectors::queryHits(overlaps), "expression"
    ] <- exprs[S4Vectors::subjectHits(overlaps)]$expression
  }

  experiment_counts <- lapply(
    experiments,
    load_experiment_counts,
    genome_tiles = genome_tiles,
    regions = regions,
    paired = paired
  )

  experiment_metadata <- S4Vectors::DataFrame(
    name = names(experiments),
    control = names(experiments) == control_name,
    lambda_g = vapply(experiment_counts, function(x) x$lambda_g, numeric(1))
  )

  count_matrix <- do.call(cbind, lapply(experiment_counts, function(x) x$counts))

  SummarizedExperiment::SummarizedExperiment(
    rowRanges = regions,
    assays = S4Vectors::SimpleList(counts = count_matrix),
    colData = experiment_metadata,
    metadata = list(
      sample_name = sample_name,
      control_name = control_name,
      genome_tiles = genome_tiles
    )
  )
}
