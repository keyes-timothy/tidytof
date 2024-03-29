% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{tof_read_file}
\alias{tof_read_file}
\title{Read high-dimensional cytometry data from a single .fcs or .csv file into a tidy tibble.}
\usage{
tof_read_file(file_path = NULL, sep = "|", panel_info = dplyr::tibble())
}
\arguments{
\item{file_path}{A file path to a single .fcs or .csv file.}

\item{sep}{A string to use to separate the antigen name and its associated
metal in the column names of the output tibble. Defaults to "|". Only used
if the input file is an .fcs file.}

\item{panel_info}{Optional. A tibble or data.frame containing information about the
panel used during high-dimensional cytometry data acquisition. Two columns are required:
"metals" and "antigens". Only used if the input file is a .csv file.}
}
\value{
A `tof_tbl` in which each row represents a single cell and each
column represents a high-dimensional cytometry antigen channel.

A `tof_tbl` is an S3 class that extends the "tibble" class by storing
one additional attribute: "panel" (a tibble storing information about the
panel used during data acquisition). Because panel information isn't
obvious from data read as a .csv file, this information must be provided
manually by the user.
}
\description{
Read high-dimensional cytometry data from a single .fcs or .csv file into a tidy tibble.
}
