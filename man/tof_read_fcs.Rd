% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{tof_read_fcs}
\alias{tof_read_fcs}
\title{Read high-dimensional cytometry data from an .fcs file into a tidy tibble.}
\usage{
tof_read_fcs(file_path = NULL, sep = "|")
}
\arguments{
\item{file_path}{A file path to a single .fcs file.}

\item{sep}{A string to use to separate the antigen name and its associated
metal in the column names of the output tibble. Defaults to "|".}
}
\value{
a `tof_tbl` in which each row represents a single cell and each
column represents a high-dimensional cytometry antigen channel.

A `tof_tbl` is an S3 class that extends the "tibble" class by storing
one additional attribute: "panel" (a tibble storing information about the
panel used during data acquisition).
}
\description{
This function reads high-dimensional cytometry data from a single .fcs file into a tidy data
structure called a `tof_tbl` ("tof_tibble"). tof_tibbles are identical to normal
tibbles except for an additional attribute ("panel") that stores information
about the high-dimensional cytometry panel used during data acquisition.
}
