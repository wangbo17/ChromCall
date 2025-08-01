% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_load_bedfile.R
\name{load_bedfile}
\alias{load_bedfile}
\title{Load a BED File as a GRanges Object}
\usage{
load_bedfile(file, genome = NULL, reduce = FALSE)
}
\arguments{
\item{file}{Path to a BED-format file with at least three tab-delimited columns.}

\item{genome}{Optional \link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb::Seqinfo} object used to annotate and filter chromosomes.}

\item{reduce}{Logical; if \code{TRUE}, the resulting \link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges} object will be reduced
using \code{\link[GenomicRanges:inter-range-methods]{GenomicRanges::reduce()}}.}
}
\value{
A \link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges} object.
}
\description{
Load a BED-format file into a \link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges} object. The file must contain
at least three columns: chromosome (\code{chr}), start, and end positions. Only the first
three columns are used; others are ignored. If a \link[GenomeInfoDb:Seqinfo-class]{GenomeInfoDb::Seqinfo} object is
provided via \code{genome}, only chromosomes present in that object are retained.
}
\examples{
bed_file <- system.file("extdata", "example.bed", package = "chromcall")
genome <- load_genome(system.file("extdata", "genome.txt", package = "chromcall"))
bed_regions <- load_bedfile(bed_file, genome = genome, reduce = TRUE)
bed_regions

}
