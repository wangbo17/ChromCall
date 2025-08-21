#' Load a BED File as a GRanges Object
#'
#' @description
#' Load a BED-format file into a [GenomicRanges::GRanges] object. The file must contain
#' at least three columns: chromosome (`chr`), start, and end positions. Only the first
#' three columns are used; others are ignored. If a [GenomeInfoDb::Seqinfo] object is
#' provided via `genome`, only chromosomes present in that object are retained.
#'
#' @param file Path to a BED-format file with at least three tab-delimited columns.
#' @param genome Optional [GenomeInfoDb::Seqinfo] object used to annotate and filter chromosomes.
#' @param reduce Logical; if `TRUE`, the resulting [GenomicRanges::GRanges] object will be reduced
#' using [GenomicRanges::reduce()].
#'
#' @return A [GenomicRanges::GRanges] object.
#'
#' @examples
#' bed_file <- system.file("extdata", "example.bed", package = "chromcall")
#' genome <- load_genome(system.file("extdata", "genome.txt", package = "chromcall"))
#' bed_regions <- load_bedfile(bed_file, genome = genome, reduce = TRUE)
#' bed_regions
#'
#' @export
load_bedfile <- function(file, genome = NULL, reduce = FALSE) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  bed_coltypes <- sapply(
    utils::read.table(
      file,
      sep = "\t",
      header = FALSE,
      nrows = 1
    ),
    class
  )

  bed <- utils::read.table(
    file,
    sep = "\t",
    header = FALSE,
    colClasses = bed_coltypes,
    comment.char = "#"
  )

  if (ncol(bed) < 3) {
    stop("BED file must contain at least 3 columns: chr, start, and end.")
  }

  bed <- bed[, 1:3, drop = FALSE]
  colnames(bed) <- c("chr", "start", "end")

  if (!is.null(genome)) {
    bed <- bed[bed$chr %in% GenomeInfoDb::seqnames(genome), , drop = FALSE]
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = bed,
    seqinfo = genome,
    starts.in.df.are.0based = F
  )

  if (isTRUE(reduce)) {
    gr <- GenomicRanges::reduce(gr)
  }

  return(gr)
}
