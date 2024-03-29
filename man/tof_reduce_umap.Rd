% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensionality_reduction.R
\name{tof_reduce_umap}
\alias{tof_reduce_umap}
\title{Apply uniform manifold approximation and projection (UMAP) to single-cell data}
\usage{
tof_reduce_umap(
  tof_tibble,
  umap_cols = where(tof_is_numeric),
  num_comp = 2,
  neighbors = 5,
  min_dist = 0.01,
  learn_rate = 1,
  epochs = NULL,
  verbose = FALSE,
  n_threads = 1,
  return_recipe = FALSE,
  ...
)
}
\arguments{
\item{tof_tibble}{A `tof_tbl` or `tibble`.}

\item{umap_cols}{Unquoted column names indicating which columns in `tof_tibble` to
use in computing the UMAP embedding. Defaults to all numeric columns
in `tof_tibble`. Supports tidyselect helpers.}

\item{num_comp}{An integer for the number of UMAP components.}

\item{neighbors}{An integer for the number of nearest neighbors used to
construct the target simplicial set.}

\item{min_dist}{The effective minimum distance between embedded points.}

\item{learn_rate}{Positive number of the learning rate for the optimization
process.}

\item{epochs}{Number of iterations for the neighbor optimization.
See \code{\link[uwot]{umap}} for details.}

\item{verbose}{A boolean indicating if run details should be logged to the
console. Defaults to FALSE.}

\item{n_threads}{Number of threads to use during UMAP calculation. Defaults
to 1.}

\item{return_recipe}{A boolean value indicating if instead of the UMAP result, a
prepped \code{\link[recipes]{recipe}} object
containing the UMAP embedding should be
returned. Set this option to TRUE if you want to create the UMAP embedding using
one dataset but also want to project new observations onto the same embedding
space later.}

\item{...}{Optional. Other options to be passed as arguments to \code{\link[uwot]{umap}}.}
}
\value{
A tibble with the same number of rows as `tof_tibble`, each representing
a single cell. Each of the `num_comp` columns represents each cell's embedding
in the calculated UMAP space.
}
\description{
This function calculates a UMAP embedding from single-cell data in a `tof_tibble`.
}
\examples{
# simulate single-cell data
sim_data <-
    dplyr::tibble(
        cd45 = rnorm(n = 200),
        cd38 = rnorm(n = 200),
        cd34 = rnorm(n = 200),
        cd19 = rnorm(n = 200)
    )
new_data <-
    dplyr::tibble(
        cd45 = rnorm(n = 50),
        cd38 = rnorm(n = 50),
        cd34 = rnorm(n = 50),
        cd19 = rnorm(n = 50)
    )

# calculate umap
tof_reduce_umap(tof_tibble = sim_data)

# calculate umap with only 2 columns
tof_reduce_tsne(tof_tibble = sim_data, umap_cols = c(cd34, cd38))

# return recipe
umap_recipe <- tof_reduce_umap(tof_tibble = sim_data, return_recipe = TRUE)

# apply recipe to new data
recipes::bake(umap_recipe, new_data = new_data)

}
\seealso{
Other dimensionality reduction functions: 
\code{\link{tof_reduce_dimensions}()},
\code{\link{tof_reduce_pca}()},
\code{\link{tof_reduce_tsne}()}
}
\concept{dimensionality reduction functions}
