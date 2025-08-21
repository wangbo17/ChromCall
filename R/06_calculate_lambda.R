#' Calculate Genomic Lambda Value
#'
#' @description
#' Computes the genomic lambda value from a [GenomicRanges::GRanges] object by summarizing a specified counts column.
#' Optionally excludes blacklisted regions and zero-count entries before summarization.
#'
#' @param x A [GenomicRanges::GRanges] object containing count values and optional metadata.
#' @param counts_col Character. Column name containing count values. Default is `"counts"`.
#' @param blacklist_col Character (optional). Logical column indicating blacklisted regions.
#' @param rm_zero Logical. Whether to exclude zero-count regions before summarization. Default is TRUE.
#' @param summary_fun Function. Summary function to apply (e.g., [mean()], [median()]). Default is [mean()].
#'
#' @return A single numeric value representing the background lambda level.
#'
#' @examples
#' # GRanges object with counts and blacklist
#' gr <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = c(1, 1001, 2001), width = 1000),
#'   counts = c(8, 0, 30),
#'   blacklist = c(FALSE, FALSE, TRUE)
#' )
#'
#' # No blacklist filtering: 'blacklist_col' not specified
#' calculate_lambda(gr)
#'
#' # With blacklist filtering: exclude third region
#' calculate_lambda(gr, blacklist_col = "blacklist")
#'
#' # With blacklist and zero filtering: exclude third (blacklist) and second (zero)
#' calculate_lambda(gr, blacklist_col = "blacklist", rm_zero = TRUE)
#'
#' # Use median instead of mean (no blacklist)
#' calculate_lambda(gr, summary_fun = median)
#'
#' @export
calculate_lambda <- function(
    x,
    counts_col = "counts",
    blacklist_col,
    rm_zero = TRUE,
    summary_fun = mean
) {

  counts <- S4Vectors::mcols(x)[[counts_col]]

  if (!missing(blacklist_col)) {
    counts <- counts[!S4Vectors::mcols(x)[[blacklist_col]]]
  }

  if (isTRUE(rm_zero)) {
    counts <- counts[counts != 0]
  }

  summary_fun(counts)
}
