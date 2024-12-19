#' Find the Best Bin with a Fixed Pseudocount
#'
#' This function identifies which bin has the strongest (absolute) correlation
#' with gene expression levels, given coverage data and a fixed pseudocount.
#'
#' @param hist_cov A named list of coverage Rle objects by chromosome.
#' @param bins_list A list of GRanges objects, one per gene, representing gene bins.
#' @param expression A numeric vector of gene expression values, one per gene.
#' @return A list with:
#' \itemize{
#'   \item \code{pseudocount}: The fixed pseudocount (0.1).
#'   \item \code{best_bin}: The index of the bin with the strongest absolute correlation.
#'   \item \code{correlations}: A numeric vector of correlations (one per bin).
#' }
#'
#' @details
#' We first compute mean signals for each gene's bins using \code{\link{mean_bin_signal}}.
#' Then we add a pseudocount (0.1), log2-transform, and correlate each bin's coverage
#' profile with the expression vector.
#'
#' @importFrom stats cor
#' @export
find_best_bin_fixed_pc <- function(hist_cov, bins_list, expression) {
  pc <- 0.1
  signals_list <- lapply(bins_list, mean_bin_signal, cov=hist_cov)
  signal_matrix <- do.call(rbind, signals_list)
  log_mat <- log2(signal_matrix + pc)
  cors <- apply(log_mat, 2, function(x) cor(x, expression, use="complete.obs"))
  current_best_bin_idx <- which.max(abs(cors))
  list(pseudocount = pc, best_bin = current_best_bin_idx, correlations = cors)
}
