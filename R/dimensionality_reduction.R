# dimensionality_reduction.R
# This file contains functions relevant to performing dimensionality reduction
# on tof_tibble objects containing CyTOF data.

# tof_reduce_pca ----------------------------
#
#' Perform principal component analysis on CyTOF data
#'
#' This function calculates principal components using single-cell data from a `tof_tibble`.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param pca_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use for computing the principal components. Defaults to all numeric columns.
#' Supports tidyselect helpers.
#'
#' @param num_comp The number of PCA components to calculate. Defaults
#' to 5. See \code{\link[recipes]{step_pca}}.
#'
#' @param threshold A double between 0 and 1 representing the fraction of total
#' variance that should be covered by the components returned in the output. See
#' \code{\link[recipes]{step_pca}}.
#'
#' @param center A boolean value indicating if each column should be centered to
#' mean 0 before PCA analysis. Defaults to TRUE.
#'
#' @param scale A boolean value indicating if each column should be scaled to
#' standard deviation = 1 before PCA analysis. Defaults to TRUE.
#'
#' @return A tibble with the same number of rows as `tof_tibble`, each representing
#' a single cell. Each of the `num_comp` columns represents each cell's embedding
#' in the calculated principal component space.
#'
#' @export
#'
#'
#' @importFrom recipes recipe
#' @importFrom recipes all_numeric
#' @importFrom recipes prep
#' @importFrom recipes juice
#' @importFrom recipes step_pca
#'
#'
tof_reduce_pca <-
  function(
    tof_tibble,
    pca_cols = where(tof_is_numeric),
    num_comp = 5,
    threshold = NA,
    center = TRUE,
    scale = TRUE
  ) {

    recipes::recipe(~ ., data = select(tof_tibble, {{pca_cols}})) %>%
      # remove any variables that have 0 variance
      recipes::step_zv(recipes::all_numeric()) %>%
      recipes::step_pca(
        recipes::all_numeric(),
        num_comp = num_comp,
        threshold = threshold,
        options = list(center = center, scale. = scale)
      ) %>%
      recipes::prep() %>%
      recipes::juice() %>%
      dplyr::rename_with(.fn = ~ paste0(".", .x))
  }





# tof_reduce_tsne ----------------------------

#' Perform t-distributed stochastic neighborhood embedding on CyTOF data
#'
#' This function calculates a tSNE embedding using single-cell data from a `tof_tibble`.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param tsne_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the tSNE embedding. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_comp The number of tSNE components to calculate for the embedding.
#' Defaults to 2.
#'
#' @param perplexity A positive numeric value that represents represents the rough
#' balance between the input data’s local and global structure emphasized in
#' the embedding. Smaller values emphasize local structure; larger values emphasize
#' global structure. The recommended range is generally 5-50. Defaults to 30.
#'
#' @param theta A numeric value representing the speed/accuracy tradeoff for the
#' embedding. Set to 0 for the exact tSNE; increase for a faster approximation.
#' Defaults to 0.5
#'
#' @param max_iterations An integer number of iterations to use during embedding
#' calculation. Defaults to 1000.
#'
#' @param verbose A boolean value indicating whether progress updates should be
#' printed during embedding calculation. Default is FALSE.
#'
#' @param ... Additional arguments to pass to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A tibble with the same number of rows as `tof_tibble`, each representing
#' a single cell. Each of the `num_comp` columns represents each cell's embedding
#' in the calculated tSNE space.
#'
#' @export
#'
#'
#' @importFrom dplyr as_tibble
#' @importFrom purrr pluck
#'
tof_reduce_tsne <-
  function(
    tof_tibble,
    tsne_cols = where(tof_is_numeric),
    num_comp = 2,
    perplexity = 30,
    theta = 0.5,
    max_iterations = 1000,
    verbose = FALSE,
    ...
  ) {

    # check that Rtsne is installed
    has_rtsne <- requireNamespace(package = "Rtsne")
    if (!has_rtsne) {
      stop(
        "This function requires the {Rtsne} package. Install it with this code:\n
           install.packages(\"Rtsne\")"
      )
    }

    result <-
      Rtsne::Rtsne(
        X = as.matrix(select(tof_tibble, {{tsne_cols}})),
        dims = num_comp,
        perplexity = perplexity,
        theta = theta,
        check_duplicates = FALSE,
        max_iter = max_iterations,
        verbose = verbose,
        ...
      ) %>%
      purrr::pluck("Y") %>%
      dplyr::as_tibble(.name_repair = "minimal")

    colnames(result) <- paste0(".tsne_", 1:num_comp)

    return(result)

  }




# tof_reduce_umap --------------------------------------------------------------

#' Perform uniform manifold approzimation and projection on CyTOF data
#'
#' This function calculates a UMAP embedding from single-cell data in a `tof_tibble`.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param umap_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the UMAP embedding. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_comp An integer for the number of UMAP components.
#'
#' @param neighbors An integer for the number of nearest neighbors used to
#' construct the target simplicial set.
#'
#' @param min_dist The effective minimum distance between embedded points.
#'
#' @param learn_rate Positive number of the learning rate for the optimization
#' process.
#'
#' @param epochs Number of iterations for the neighbor optimization.
#' See \code{\link[uwot]{umap}} for details.
#'
#' @param verbose A boolean indicating if run details should be logged to the
#' console. Defaults to FALSE.
#'
#' @param n_threads Number of threads to use during UMAP calculation. Defaults
#' to 1.
#'
#' @param ... Optional. Other options to be passed as arguments to \code{\link[uwot]{umap}}.
#'
#' @return A tibble with the same number of rows as `tof_tibble`, each representing
#' a single cell. Each of the `num_comp` columns represents each cell's embedding
#' in the calculated UMAP space.
#'
#' @export
#'
#'
#' @importFrom embed step_umap
#'
#'
tof_reduce_umap <-
  function(
    tof_tibble,
    umap_cols = where(tof_is_numeric),
    num_comp = 2,
    neighbors = 5,
    min_dist = 0.01,
    learn_rate = 1,
    epochs = NULL,
    verbose = FALSE,
    n_threads = 1,
    ...
  ) {

    suppressWarnings(
      recipes::recipe(~ ., data = select(tof_tibble, {{umap_cols}})) %>%
        # remove any variables that have 0 variance
        recipes::step_zv(recipes::all_numeric()) %>%
        embed::step_umap(
          recipes::all_numeric(),
          num_comp = num_comp,
          neighbors = neighbors,
          min_dist = min_dist,
          learn_rate = learn_rate,
          epochs = epochs,
          options = list(verbose = verbose, n_threads = n_threads, ...)
        ) %>%
        recipes::prep() %>%
        recipes::juice() %>%
        dplyr::rename_with(.fn = ~ paste0(".", .x))
    )
  }



# tof_reduce_dimensions --------------------------------------------------------

#' Apply dimensionality reduction to a CyTOF dataset.
#'
#' This function is a wrapper around {tidytof}'s tof_reduce_* function family.
#' It performs dimensionality reduction on CyTOF data using a user-specified method
#' (of 3 choices) and each method's corresponding input parameters
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param method A method of dimensionality reduction. Currently, PCA, tSNE, and
#' UMAP embedding are supported.
#'
#' @param ... Arguments to be passed to the tof_reduce_* function corresponding to
#' the embedding method. See \code{\link{tof_reduce_pca}}, \code{\link{tof_reduce_tsne}}, and
#' \code{\link{tof_reduce_umap}}.
#'
#' @return A tibble with the same number of rows as `tof_tibble`, each representing
#' a single cell. Each of the `num_comp` columns represents each cell's embedding
#' in the calculated embedding space.
#'
#' @export
#'
#'
tof_reduce_dimensions <- function(tof_tibble, method = c("pca", "tsne", "umap"), add_cols = TRUE, ...) {

  # check validity of method
  method <- match.arg(method, choices = c("pca", "tsne", "umap"))

  if (method == "pca") {
    result <- tof_reduce_pca(tof_tibble, ...)
  } else if (method == "tsne") {
    result <- tof_reduce_tsne(tof_tibble, ...)
  } else if (method == "umap") {
    result <- tof_reduce_umap(tof_tibble, ...)
  } else {
      stop("Method no implemented")
  }

  if (add_cols == TRUE) {
    result <-
      bind_cols(tof_tibble, result)
  }

  return(result)

}





