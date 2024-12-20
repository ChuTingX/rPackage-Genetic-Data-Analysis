#' Generate Example Data
#'
#' This function generates example data used to demonstrate the package functions.
#' It creates a set of genomic bins (`bins_list`), a coverage object (`example_cov`),
#' and an expression vector (`expression`), all from scratch each time it is called.
#'
#' @param seed A numeric seed for reproducibility. Set to a fixed value by default.
#' @return A list with elements: bins_list, example_cov, and expression.
#' @import GenomicRanges IRanges S4Vectors
#' @export
generate_example_data <- function(seed = 123) {
  set.seed(seed)
  num_genes <- 5
  num_bins <- 10
  bin_width <- 100

  bins_list <- vector("list", num_genes)
  for (g in seq_len(num_genes)) {
    start_bin <- (g-1)*(num_bins*bin_width) + 1
    end_bin <- g*(num_bins*bin_width)
    bin_starts <- seq(start_bin, end_bin, by=bin_width)
    bin_ends <- bin_starts + bin_width - 1
    bins_list[[g]] <- GRanges(
      seqnames = Rle("chr1"),
      ranges = IRanges(start=bin_starts, end=bin_ends),
      strand = Rle("*")
    )
  }

  total_length <- num_genes * num_bins * bin_width
  cov_chr1 <- Rle(runif(total_length, 0, 10))
  example_cov <- list(chr1=cov_chr1)

  # Use the package's mean_bin_signal function to compute signals
  # If mean_bin_signal is defined in the package
  temp_signals_list <- lapply(bins_list, mean_bin_signal, cov=example_cov)
  signal_matrix_temp <- do.call(rbind, temp_signals_list)
  expression <- 2 * signal_matrix_temp[,5] + rnorm(num_genes, 0, 0.5)

  list(
    bins_list = bins_list,
    example_cov = example_cov,
    expression = expression
  )
}
