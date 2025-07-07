#' Load Expression Data from TSS BED File
#'
#' @description
#' Load a BED-format file containing transcription start site (TSS) coordinates and expression values,
#' and convert it into a [GenomicRanges::GRanges] object with an "expression" metadata column.
#'
#' @param file Path to a BED-format file with at least 3 columns: chromosome, TSS, and expression.
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
  bed <- utils::read.table(
    file,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    colClasses = c("character", "numeric", "numeric")
  )

  stopifnot(
    "BED file must have at least 3 columns: chr, tss, expression" = ncol(bed) >= 3
  )

  bed <- bed[, 1:3]
  colnames(bed) <- c("chr", "tss", "expression")

  if (!is.null(genome)) {
    bed <- bed[bed$chr %in% GenomeInfoDb::seqnames(genome), , drop = FALSE]
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    bed,
    start.field = "tss",
    end.field = "tss",
    keep.extra.columns = TRUE,
    seqinfo = genome
  )

  return(gr)
}
