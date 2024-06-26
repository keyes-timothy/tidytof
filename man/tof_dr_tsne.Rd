% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensionality_reduction.R
\name{tof_dr_tsne}
\alias{tof_dr_tsne}
\title{Perform t-distributed stochastic neighborhood embedding on CyTOF data}
\usage{
tof_dr_tsne(
  tof_tibble,
  tsne_vars = where(tof_is_numeric),
  num_comp = 2,
  perplexity = 30,
  theta = 0.5,
  max_iterations = 1000,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{tof_tibble}{A `tof_tibble`.}

\item{tsne_vars}{Unquoted column names indicating which columns in `tof_tibble` to
use in computing the tSNE embedding. Defaults to all numeric columns
in `tof_tibble`. Supports tidyselect helpers.}

\item{num_comp}{The number of tSNE components to calculate for the embedding.
Defaults to 2.}

\item{perplexity}{A positive numeric value that represents represents the rough
balance between the input data’s local and global structure emphasized in
the embedding. Smaller values emphasize local structure; larger values emphasize
global structure. The recommended range is generally 5-50. Defaults to 30.}

\item{theta}{A numeric value representing the speed/accuracy tradeoff for the
embedding. Set to 0 for the exact tSNE; increase for a faster approximation.
Defaults to 0.5}

\item{max_iterations}{An integer number of iterations to use during embedding
calculation. Defaults to 1000.}

\item{verbose}{A boolean value indicating whether progress updates should be
printed during embedding calculation. Default is FALSE.}

\item{...}{Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.}
}
\value{
A tibble with the same number of rows as `tof_tibble`, each representing
a single cell. Each of the `num_comp` columns represents each cell's embedding
in the calculated tSNE space.
}
\description{
This function calculates a tSNE embedding using single-cell data from a `tof_tibble`.
}
\examples{
NULL

}
