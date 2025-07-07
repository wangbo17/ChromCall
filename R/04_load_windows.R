#' Load Promoter Window Regions from a BED File
#'
#' @description
#' Load genomic regions (e.g., promoter windows) from a BED-format file and associate
#' them with a given [GenomeInfoDb::Seqinfo] genome object. If a blacklist is provided,
#' an additional logical column `blacklist` is added to indicate whether each region
#' overlaps with any blacklisted region.
#'
#' @param file Path to a BED-format file specifying genomic windows (e.g., promoters).
#' @param genome A [GenomeInfoDb::Seqinfo] object used to define the genome context.
#' @param blacklist Optional [GenomicRanges::GRanges] object. If provided, each region
#' will be annotated with a logical flag `blacklist = TRUE` if it overlaps any blacklist region.
#'
#' @return A [GenomicRanges::GRanges] object. If `blacklist` is provided, it includes an
#' additional metadata column `blacklist`.
#'
#' @examples
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' genome <- load_genome(genome_file)
#' blacklist_file <- system.file("extdata", "blacklist.bed", package = "chromcall")
#' blacklist <- load_bedfile(blacklist_file, genome = genome)
#' promoter_file <- system.file("extdata", "example.bed", package = "chromcall")
#' promoters <- load_windows(promoter_file, genome, blacklist = blacklist)
#' promoters
#'
#' @export
load_windows <- function(file, genome, blacklist = NULL) {
  w <- load_bedfile(file, genome = genome)

  if (!is.null(blacklist) && inherits(blacklist, "GRanges")) {
    w$blacklist <- IRanges::overlapsAny(w, blacklist)
  }

  return(w)
}
