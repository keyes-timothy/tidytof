% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tof_find_jsd}
\alias{tof_find_jsd}
\title{Find the Jensen-Shannon Divergence (JSD) between two numeric vectors}
\usage{
tof_find_jsd(vec_1, vec_2, num_bins = 100)
}
\arguments{
\item{vec_1}{A numeric vector.}

\item{vec_2}{A numeric vector.}

\item{num_bins}{An integer number of bins to use when binning
across the two vectors' combined range. Defaults to 100.}
}
\value{
A double (of length 1) representing the JSD between the two vectors.
}
\description{
Find the Jensen-Shannon Divergence (JSD) between two numeric vectors
}
