% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_best_bin_fixed_pc.R
\name{find_best_bin_fixed_pc}
\alias{find_best_bin_fixed_pc}
\title{Find the Best Bin with a Fixed Pseudocount}
\usage{
find_best_bin_fixed_pc(hist_cov, bins_list, expression)
}
\arguments{
\item{hist_cov}{A named list of coverage Rle objects by chromosome.}

\item{bins_list}{A list of GRanges objects, one per gene, representing gene bins.}

\item{expression}{A numeric vector of gene expression values, one per gene.}
}
\value{
A list with:
\itemize{
\item \code{pseudocount}: The fixed pseudocount (0.1).
\item \code{best_bin}: The index of the bin with the strongest absolute correlation.
\item \code{correlations}: A numeric vector of correlations (one per bin).
}
}
\description{
This function identifies which bin has the strongest (absolute) correlation
with gene expression levels, given coverage data and a fixed pseudocount.
}
\details{
We first compute mean signals for each gene's bins using \code{\link{mean_bin_signal}}.
Then we add a pseudocount (0.1), log2-transform, and correlate each bin's coverage
profile with the expression vector.
}
