#' Load Genome Information from a 4-Column Text File
#'
#' @description
#' Load genome sequence information from a tab-delimited text file into a
#' [GenomeInfoDb::Seqinfo] object. The file must contain four columns in the following order:
#' chromosome name (`chr`), chromosome length (`length`), circularity (`isCircular`),
#' and genome build (`genome`).
#'
#' @param file Path to a tab-delimited text file with four columns: chromosome name
#' (character), length (numeric), circularity (logical), and genome build (character).
#'
#' @return A [GenomeInfoDb::Seqinfo] object.
#'
#' @examples
#' genome_file <- system.file("extdata", "genome.txt", package = "chromcall")
#' genome <- load_genome(genome_file)
#' genome
#'
#' @export
load_genome <- function(file) {
  genome <- utils::read.table(
    file,
    sep = "\t",
    header = FALSE,
    colClasses = c("character", "numeric", "logical", "character"),
    col.names = c("chr", "length", "isCircular", "genome")
  )

  GenomeInfoDb::Seqinfo(
    seqnames = genome$chr,
    seqlengths = genome$length,
    isCircular = genome$isCircular,
    genome = genome$genome
  )
}
