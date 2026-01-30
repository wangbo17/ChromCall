#' Load Expression Data from TSS BED-like File
#'
#' @description
#' Load a tab-delimited file containing transcription start site (TSS) coordinates and expression values,
#' and convert it into a [GenomicRanges::GRanges] object with an "expression" metadata column.
#'
#' **Coordinate convention:** This function assumes the input `tss` column is **1-based**.
#' This matches the current ChromCall window files used in the project.
#'
#' @param file Path to a tab-delimited file with at least 3 columns: chromosome, TSS, and expression.
#' @param genome Optional [GenomeInfoDb::Seqinfo] object to restrict chromosomes and assign seqinfo.
#'
#' @return A [GenomicRanges::GRanges] object with 1-bp ranges and a metadata column `expression`.
#'
#' @examples
#' expr_file <- system.file("extdata", "expression_tss.bed", package = "chromcall")
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' genome <- load_genome(genome_file)
#' expr_gr <- load_expression(expr_file, genome = genome)
#' expr_gr
#'
#' @export
load_expression <- function(file, genome = NULL) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  bed <- utils::read.table(
    file,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    colClasses = c("character", "numeric", "numeric")
  )

  stopifnot(
    "Input must have at least 3 columns: chr, tss, expression" = ncol(bed) >= 3
  )

  bed <- bed[, 1:3, drop = FALSE]
  colnames(bed) <- c("chr", "tss", "expression")

  if (!is.null(genome)) {
    bed <- bed[bed$chr %in% GenomeInfoDb::seqnames(genome), , drop = FALSE]
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    bed,
    start.field = "tss",
    end.field = "tss",
    keep.extra.columns = TRUE,
    seqinfo = genome,
    starts.in.df.are.0based = FALSE
  )

  gr
}
