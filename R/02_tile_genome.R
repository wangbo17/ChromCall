#' Generate Fixed-Width Windows by Tiling a Genome
#'
#' @description
#' Generate fixed-size, non-overlapping windows (tiles) across a genome defined by a
#' [GenomeInfoDb::Seqinfo] object. If a blacklist is provided, each tile is annotated
#' with a logical `blacklist` column indicating whether it overlaps any blacklisted region.
#'
#' @param genome A [GenomeInfoDb::Seqinfo] object defining the genome to tile.
#' @param window_size Integer. Width of each tile (e.g., 10000 for 10kb windows).
#' @param blacklist Optional [GenomicRanges::GRanges] object. If provided, overlapping tiles
#' will be annotated with `blacklist = TRUE`.
#'
#' @return A [GenomicRanges::GRanges] object representing tiled windows. Includes a
#' `blacklist` column if `blacklist` is provided.
#'
#' @examples
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' genome <- load_genome(genome_file)
#'
#' # Generate 10kb windows without blacklist
#' tiles <- tile_genome(genome, window_size = 10000)
#' tiles
#'
#' # Generate 10kb windows with blacklist
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#' blacklist <- load_bedfile(blacklist_file, genome = genome)
#' tiles_with_bl <- tile_genome(genome, window_size = 10000, blacklist = blacklist)
#' tiles_with_bl
#'
#' @export
tile_genome <- function(genome, window_size, blacklist = NULL) {

  window_size <- as.integer(window_size)
  stopifnot(
    "window_size must not be NA" = !is.na(window_size),
    "window_size must be positive" = window_size > 0
  )

  tiles <- GenomicRanges::tileGenome(
    seqlengths = genome,
    tilewidth = window_size,
    cut.last.tile.in.chrom = TRUE
  )

  if (!is.null(blacklist) && inherits(blacklist, "GRanges")) {
    tiles$blacklist <- IRanges::overlapsAny(tiles, blacklist)
  }

  return(tiles)
}
