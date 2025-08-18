#' Load Reads from a BAM File as Central Points
#'
#' @description
#' Load aligned sequencing reads from a BAM file and convert each read into a
#' single central point, returned as a [GenomicRanges::GRanges] object.
#' This is commonly used in footprinting or single-base signal profiling,
#' where only the center of each fragment is needed.
#'
#' @param file Path to a BAM file. Must be coordinate-sorted and indexed (BAI file required).
#' @param paired Logical. If TRUE (default), input BAM is assumed to be paired-end
#'   and reads are loaded as fragments. If FALSE, input BAM is treated as single-end.
#'
#' @return A [GenomicRanges::GRanges] object, with each read resized to a single base
#' at its center.
#'
#' @examples
#' bam_file <- system.file("extdata", "example.bam", package = "chromcall")
#' read_points <- load_read_locations(bam_file, paired = TRUE)
#' read_points
#'
#' @export
load_read_locations <- function(file, paired) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  if (paired) {
    gal <- GenomicAlignments::readGAlignmentPairs(file)
    gr <- GenomicRanges::granges(gal)
  } else {
    gal <- GenomicAlignments::readGAlignments(file)
    gr <- GenomicRanges::granges(gal)
  }

  GenomicRanges::resize(
    gr,
    width = 1,
    fix = "center"
  )
}
