#' Create a chromcall sample from input files
#'
#' @description
#' Constructs a single-sample [SummarizedExperiment::SummarizedExperiment] object by loading
#' one or more experiments (e.g., histone marks, CUT&RUN, or ATAC-seq) with a shared genome,
#' promoter windows, and optional blacklist information. For each experiment, it counts reads
#' in the target regions and estimates a genome-wide background lambda using tiled genome windows.
#'
#' By default, ChromCall expects a matched input/control track (specified via `control_name`).
#' If `control_name` is `NULL`, ChromCall will run in *no-control mode*, which is useful for data
#' types such as ATAC-seq. In no-control mode, downstream testing should assume background
#' weighting is constant (effectively 1), i.e. no input-based correction is applied.
#'
#' @param sample_name The name of the sample (string).
#' @param experiments A named list of BAM file paths for each experiment
#'   (e.g., \code{list(H3K27me3="...", Control="...")} or for ATAC-seq \code{list(ATAC="...")}).
#' @param paired Logical. If TRUE (default), input BAM is assumed to be paired-end and reads are
#'   loaded as fragments. If FALSE, input BAM is treated as single-end.
#' @param control_name The name (matching \code{names(experiments)}) of the experiment to treat
#'   as control. Set to `NULL` to enable no-control mode.
#' @param genome_file Path to the genome info text file.
#' @param region_file Path to the BED file defining target regions (e.g., promoters).
#' @param window_size Genome tile size in base pairs.
#' @param blacklist_file Optional path to a BED file of blacklist regions, or NULL to disable.
#' @param expression_file Optional path to a BED-format TSS expression file to merge expression into regions.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] object containing region counts and lambda estimates.
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
#' # Explore output (with control)
#' sample
#' SummarizedExperiment::assayNames(sample)
#' SummarizedExperiment::colData(sample)
#' head(SummarizedExperiment::assay(sample, "counts"))
#' head(SummarizedExperiment::rowData(sample)$expression)
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
#' # Explore output (no-control)
#' atac_sample
#' SummarizedExperiment::assayNames(atac_sample)
#' SummarizedExperiment::colData(atac_sample)
#' head(SummarizedExperiment::assay(atac_sample, "counts"))
#'
#' # After region-level testing (no-control)
#' atac_res <- test_region_counts(atac_sample)
#' SummarizedExperiment::assayNames(atac_res)
#' head(SummarizedExperiment::assay(atac_res, "score"))
#'
#' @export
build_chromcall_sample <- function(
    sample_name,
    experiments,
    paired = TRUE,
    control_name = NULL,
    genome_file,
    region_file,
    window_size,
    blacklist_file = NULL,
    expression_file = NULL
) {
  stopifnot(
    "experiments must be a named list" = is.list(experiments) && !is.null(names(experiments)),
    "experiments must have at least one entry" = length(experiments) >= 1
  )

  # Basic experiment name/path checks
  if (any(!nzchar(names(experiments)))) {
    stop("All experiments must have non-empty names.")
  }
  if (any(!nzchar(as.character(experiments)))) {
    stop("All experiments must have non-empty BAM paths.")
  }
  missing_bams <- !file.exists(as.character(experiments))
  if (any(missing_bams)) {
    stop("BAM file(s) not found: ", paste(as.character(experiments)[missing_bams], collapse = ", "))
  }

  # Decide mode
  if (!is.null(control_name) && !is.character(control_name)) {
    stop("control_name must be a character string or NULL.")
  }
  no_control <- is.null(control_name) || is.na(control_name) || !nzchar(control_name)

  # If control specified, validate it exists
  if (!no_control && !(control_name %in% names(experiments))) {
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
      query   = GenomicRanges::resize(regions, width = 1, fix = "center"),
      subject = exprs
    )

    S4Vectors::mcols(regions)$expression <- 0
    S4Vectors::mcols(regions)[S4Vectors::queryHits(overlaps), "expression"] <-
      exprs[S4Vectors::subjectHits(overlaps)]$expression
  }

  experiment_counts <- lapply(
    experiments,
    load_experiment_counts,
    genome_tiles = genome_tiles,
    regions      = regions,
    paired       = paired
  )

  count_matrix <- do.call(cbind, lapply(experiment_counts, function(x) x$counts))
  colnames(count_matrix) <- names(experiments)

  experiment_metadata <- S4Vectors::DataFrame(
    name     = names(experiments),
    control  = if (no_control) rep(FALSE, length(experiments)) else names(experiments) == control_name,
    lambda_g = vapply(experiment_counts, function(x) x$lambda_g, numeric(1))
  )
  rownames(experiment_metadata) <- names(experiments)

  SummarizedExperiment::SummarizedExperiment(
    rowRanges = regions,
    assays    = S4Vectors::SimpleList(counts = count_matrix),
    colData   = experiment_metadata,
    metadata  = list(
      sample_name  = sample_name,
      control_name = if (no_control) NULL else control_name,
      no_control   = no_control,
      genome_tiles = genome_tiles
    )
  )
}
