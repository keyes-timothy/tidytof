% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tof_tbl.R
\name{tof_split_tidytof_reduced_dimensions}
\alias{tof_split_tidytof_reduced_dimensions}
\title{Split the dimensionality reduction data that tidytof combines during \code{\link[SingleCellExperiment]{SingleCellExperiment}} conversion}
\usage{
tof_split_tidytof_reduced_dimensions(sce)
}
\arguments{
\item{sce}{A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with an
entry named "tidytof_reduced_dimensions" in its \code{\link[SingleCellExperiment]{reducedDims}} slot.}
}
\value{
A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with separate entries
named "tidytof_pca", "tidytof_umap", and "tidytof_tsne" in its
\code{\link[SingleCellExperiment]{reducedDims}} slots (one for each of the
dimensionality reduction methods for which tidytof has native support).
}
\description{
Split the dimensionality reduction data that tidytof combines during \code{\link[SingleCellExperiment]{SingleCellExperiment}} conversion
}
\examples{
NULL
}
