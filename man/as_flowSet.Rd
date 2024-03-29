% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tof_tbl.R
\name{as_flowSet}
\alias{as_flowSet}
\alias{as_flowSet.tof_tbl}
\title{Coerce an object into a \code{\link[flowCore]{flowSet}}}
\usage{
as_flowSet(x, ...)

\method{as_flowSet}{tof_tbl}(x, group_cols, ...)
}
\arguments{
\item{x}{A tof_tbl.}

\item{...}{Unused.}

\item{group_cols}{Unquoted names of the columns in `x` that should
be used to group cells into separate \code{\link[flowCore]{flowFrame}}s.
Supports tidyselect helpers. Defaults to
NULL (all cells are written into a single \code{\link[flowCore]{flowFrame}}).}
}
\value{
A \code{\link[flowCore]{flowSet}}

A \code{\link[flowCore]{flowSet}}. Note that all non-numeric
columns in `x` will be removed.
}
\description{
Coerce an object into a \code{\link[flowCore]{flowSet}}

Coerce a tof_tbl into a \code{\link[flowCore]{flowSet}}
}
\examples{
NULL

NULL
}
