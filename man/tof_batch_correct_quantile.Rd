% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_correction.R
\name{tof_batch_correct_quantile}
\alias{tof_batch_correct_quantile}
\title{Batch-correct a tibble of high-dimensional cytometry data using quantile
normalization.}
\usage{
tof_batch_correct_quantile(
  tof_tibble,
  channel_cols,
  group_cols,
  augment = TRUE
)
}
\arguments{
\item{tof_tibble}{A `tof_tbl` or a `tibble`.}

\item{channel_cols}{Unquoted column names representing columns that contain
single-cell protein measurements. Supports tidyselect helpers.}

\item{group_cols}{Optional. Unquoted column names indicating which columns
should be used to group cells before batch correction. Batch correction is then
performed independently within each group. Supports tidyselect helpers.}

\item{augment}{A boolean value indicating if the output should replace the
`channel_cols` in `tof_tibble` with the new, batch corrected columns (TRUE, the default)
or if it should only return the batch-corrected columns (FALSE) with all other columns
omitted.}
}
\value{
If augment = TRUE, a tibble with the same number of rows and columns as
tof_tibble, with the columns specified by `channel_cols` batch-corrected. If
augment = FALSE, a tibble containing only the batch-corrected `channel_cols`.
}
\description{
This function performs quantile normalization on high-dimensional cytometry
data in tidy format using \code{\link[preprocessCore]{normalize.quantiles}}.
Optionally, groups can be specified and normalized separately.
}
\examples{
NULL

}
