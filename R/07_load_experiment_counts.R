#' Load Experiment Counts from a BAM File
#'
#' @description
#' Load sequencing reads from a BAM file, count overlaps with specified genome tiles and target regions,
#' and compute a genomic lambda value representing the background signal.
#'
#' @param file Path to the BAM file.
#' @param paired Logical. If TRUE, input BAM is assumed to be paired-end
#'   and reads are loaded as fragments. If FALSE, input BAM is treated as single-end.
#' @param genome_tiles A [GenomicRanges::GRanges] object defining genome-wide tiling regions.
#' @param regions A [GenomicRanges::GRanges] object defining target regions (e.g., promoters). Can include a 'blacklist' column.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{counts}{An integer vector of read counts for each region in \code{regions}.}
#'   \item{lambda_g}{A numeric value representing the estimated genomic background lambda.}
#' }
#'
#' @examples
#' # Load example BAM file and genome info
#' bam_file <- system.file("extdata", "example.bam", package = "chromcall")
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' genome <- load_genome(genome_file)
#' # Load blacklist
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#' blacklist <- load_bedfile(blacklist_file, genome = genome)
#'
#' # Create genome tiles (e.g., 10kb)
#' tiles <- tile_genome(genome, window_size = 10000)
#' tiles_with_bl <- tile_genome(genome, window_size = 10000, blacklist = blacklist)
#'
#' # Load promoter regions and annotate blacklist overlaps
#' promoter_file <- system.file("extdata", "example.bed", package = "chromcall")
#' promoters <- load_windows(promoter_file, genome = genome, blacklist = blacklist)
#'
#' # Run experiment counts and lambda estimation
#' res <- load_experiment_counts(
#'   file = bam_file,
#'   paired = TRUE,
#'   genome_tiles = tiles_with_bl,
#'   regions = promoters
#' )
#'
#' # Inspect results
#' res$counts     # Read counts per region
#' res$lambda_g   # Estimated lambda background
#'
#' @export
load_experiment_counts <- function(file, paired, genome_tiles, regions) {

  reads <- load_read_locations(file, paired = paired)

  tile_counts_vec <- GenomicRanges::countOverlaps(genome_tiles, reads, ignore.strand = TRUE)
  region_counts_vec <- GenomicRanges::countOverlaps(regions, reads, ignore.strand = TRUE)

  tile_has_bl <- "blacklist" %in% names(S4Vectors::mcols(genome_tiles))
  region_has_bl <- "blacklist" %in% names(S4Vectors::mcols(regions))
  if (tile_has_bl != region_has_bl) {
    stop(
      "Inconsistent blacklist columns: ",
      "genome_tiles and regions must either both include a 'blacklist' column or neither."
    )
  }

  genome_tiles_tmp <- genome_tiles
  S4Vectors::mcols(genome_tiles_tmp)$counts <- tile_counts_vec
  lambda_g <- calculate_lambda(
    genome_tiles_tmp,
    counts_col    = "counts",
    blacklist_col = if (tile_has_bl) "blacklist" else NULL,
    rm_zero       = FALSE
  )

  # optional safety: ensure lambda_g is finite and non-negative
  if (!is.finite(lambda_g) || lambda_g < 0) {
    stop("Computed lambda_g is not finite or is negative. Please check input BAM and genome_tiles.")
  }

  list(
    counts   = as.integer(region_counts_vec),
    lambda_g = lambda_g
  )
}
