#' Load Experiment Counts from a BAM File
#'
#' @description
#' Load sequencing reads from a BAM file, count overlaps with specified genome tiles and target regions,
#' and compute a genomic lambda value representing the background signal.
#'
#' @param file Path to the BAM file.
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
#'
#' # Create genome tiles (e.g., 10kb)
#' tiles <- tile_genome(genome, window_size = 10000)
#' tiles_with_bl <- tile_genome(genome, window_size = 10000, blacklist = blacklist)
#'
#' # Load blacklist
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#' blacklist <- load_bedfile(blacklist_file, genome = genome)
#'
#' # Load promoter regions and annotate blacklist overlaps
#' promoter_file <- system.file("extdata", "example.bed", package = "chromcall")
#' promoters <- load_windows(promoter_file, genome = genome, blacklist = blacklist)
#'
#' # Run experiment counts and lambda estimation
#' res <- load_experiment_counts(
#'   file = bam_file,
#'   genome_tiles = tiles_with_bl,
#'   regions = promoters
#' )
#'
#' # Inspect results
#' res$counts     # Read counts per region
#' res$lambda_g   # Estimated lambda background
#'
#' @export
load_experiment_counts <- function(file, genome_tiles, regions) {

  count_reads <- function(reads, targets) {
    targets$counts <- GenomicRanges::countOverlaps(targets, reads)
    targets
  }

  reads <- load_read_locations(file)

  tile_counts <- count_reads(reads, genome_tiles)
  region_counts <- count_reads(reads, regions)

  tile_has_bl <- "blacklist" %in% names(S4Vectors::mcols(tile_counts))
  region_has_bl <- "blacklist" %in% names(S4Vectors::mcols(region_counts))

  if (tile_has_bl != region_has_bl) {
    stop(
      "Inconsistent blacklist columns: ",
      "genome_tiles and regions must either both include a 'blacklist' column or neither."
    )
  }

  lambda_g <- calculate_lambda(
    tile_counts,
    counts_col = "counts",
    blacklist_col = if (tile_has_bl) "blacklist" else NULL,
    rm_zero = FALSE
  )

  list(
    counts = S4Vectors::mcols(region_counts)$counts,
    lambda_g = lambda_g
  )
}
