#' Load Reads from a BAM File as Central Points
#'
#' @description
#' Load aligned sequencing reads from a BAM file and convert each read into a
#' single central point, returned as a [GenomicRanges::GRanges] object.
#' This is commonly used in footprinting or single-base signal profiling,
#' where only the center of each fragment is needed.
#'
#' @param file Path to a BAM file. Must be coordinate-sorted and indexed (BAI file required).
#'
#' @return A [GenomicRanges::GRanges] object, with each read resized to a single base
#' at its center.
#'
#' @examples
#' bam_file <- system.file("extdata", "example.bam", package = "chromcall")
#' read_points <- load_read_locations(bam_file)
#' read_points
#'
#' @export
load_read_locations <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  galp <- GenomicAlignments::readGAlignmentPairs(file)
  gr <- GenomicRanges::granges(galp)

  GenomicRanges::resize(
    gr,
    width = 1,
    fix = "center"
  )
}
