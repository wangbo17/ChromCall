% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/12_compare_samples.R
\name{compare_samples}
\alias{compare_samples}
\title{Compare two chromcall samples}
\usage{
compare_samples(x, y, experiments = NULL, threshold = 0.05)
}
\arguments{
\item{x}{A chromcall sample \code{\link{RangedSummarizedExperiment}} object}

\item{y}{A chromcall sample \code{\link{RangedSummarizedExperiment}} object}

\item{experiments}{Optional character vector specifying which experiments to compare.
If not specified, compares all shared experiments.}

\item{threshold}{The adjusted p-value significance threshold (default: 0.05)}
}
\value{
A \code{\link{GRanges}} object containing per-region comparison results, including:
\itemize{
\item Per-sample adjusted p-values, class (significance), and enrichment scores (counts / expected)
\item Delta scores for each experiment (score_y - score_x)
\item Optional: expression values and \code{log2FC_expression} if available in both samples
}
}
\description{
Compare the specified experiments in two chromcall sample \code{\link{RangedSummarizedExperiment}} objects.
This function extracts adjusted p-values, assigns significance based on a threshold, computes
enrichment scores (observed / expected counts), and calculates delta scores between samples.
If both samples include gene expression values (\code{expression} column in \code{rowData}), it also adds
\verb{<sample>_expression} and \code{log2FC_expression} to the output.
}
\examples{
# Load two chromcall sample objects with expression
sampleA <- build_chromcall_sample(
  sample_name = "sampleA",
  experiments = list(
    H3K27me3 = system.file("extdata", "h3k27me3_sampleA.bam", package = "chromcall"),
    H3K4me3  = system.file("extdata", "h3k4me3_sampleA.bam", package = "chromcall"),
    Control  = system.file("extdata", "control_sampleA.bam", package = "chromcall")
  ),
  control_name = "Control",
  genome_file = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file = system.file("extdata", "example.bed", package = "chromcall"),
  window_size = 10000,
  blacklist_file = system.file("extdata", "blacklist.bed", package = "chromcall"),
  expression_file = system.file("extdata", "expression_sampleA.bed", package = "chromcall")
)

sampleB <- build_chromcall_sample(
  sample_name = "sampleB",
  experiments = list(
    H3K27me3 = system.file("extdata", "h3k27me3_sampleB.bam", package = "chromcall"),
    H3K4me3  = system.file("extdata", "h3k4me3_sampleB.bam", package = "chromcall"),
    Control  = system.file("extdata", "control_sampleB.bam", package = "chromcall")
  ),
  control_name = "Control",
  genome_file = system.file("extdata", "genome.txt", package = "chromcall"),
  region_file = system.file("extdata", "example.bed", package = "chromcall"),
  window_size = 10000,
  blacklist_file = system.file("extdata", "blacklist.bed", package = "chromcall"),
  expression_file = system.file("extdata", "expression_sampleB.bed", package = "chromcall")
)

# Run region-level differential testing
resultA <- test_region_counts(sampleA)
resultB <- test_region_counts(sampleB)

# Compare samples (returns a GRanges object)
comparison <- compare_samples(resultA, resultB, threshold = 0.05)
comparison

}
