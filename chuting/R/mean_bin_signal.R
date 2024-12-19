#' Compute Mean Signal per Bin
#'
#' Given a set of genomic bins (GRanges) and a coverage object (list of Rle),
#' this function computes the mean coverage value for each bin.
#'
#' @param bins A GRanges object representing genomic bins.
#' @param cov A named list of coverage objects (Rle) keyed by chromosome.
#' @return A numeric vector of the same length as \code{bins}, where each element
#'   is the mean coverage over that bin.
#'
#' @importFrom GenomicRanges seqnames trim
#' @importFrom IRanges Views viewMeans
#' @export
mean_bin_signal <- function(bins, cov) {
  bins_by_chr <- split(bins, seqnames(bins))
  result <- numeric(length(bins))
  start_idx <- 1
  for (chr in names(bins_by_chr)) {
    chr_bins <- bins_by_chr[[chr]]
    if (!chr %in% names(cov)) {
      # Chromosome not found in coverage
      end_idx <- start_idx + length(chr_bins) - 1
      result[start_idx:end_idx] <- NA
      start_idx <- end_idx + 1
      next
    }
    chr_cov <- cov[[chr]]
    chr_bins <- trim(chr_bins)
    if (length(chr_bins) == 0) {
      # No bins after trimming
      end_idx <- start_idx + length(chr_bins) - 1
      result[start_idx:end_idx] <- NA
      start_idx <- end_idx + 1
      next
    }
    v <- Views(chr_cov, start=start(chr_bins), end=end(chr_bins))
    means <- viewMeans(v)
    if (length(means) == 0) means <- rep(NA, length(chr_bins))
    end_idx <- start_idx + length(chr_bins) - 1
    result[start_idx:end_idx] <- means
    start_idx <- end_idx + 1
  }
  result
}
