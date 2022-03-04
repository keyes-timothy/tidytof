# downsampling.R
# This file contains functions relevant to downsampling cells within
# tof_tibble objects containing CyTOF data.

# tof_downsample_constant ------------------------------------------------------

#' Downsample CyTOF data by randomly selecting a constant number of cells per group.
#'
#' This function downsamples the number of cells in a `tof_tbl` by randomly selecting
#' `num_cells` cells from each unique combination of values in `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Unquoted names of the columns in `tof_tibble` that should
#' be used to define groups from which `num_cells` will be downsampled.
#' Supports tidyselect helpers. Defaults to `NULL` (no grouping).
#'
#' @param num_cells An integer number of cells that should be sampled from each
#' group defined by `group_cols`.
#'
#' @return A `tof_tbl` with the same number of columns as the input `tof_tibble`,
#' but fewer rows. Specifically, the number of rows will be `num_cells` multiplied
#' by the number of unique combinations of the values in `group_cols`. If any group
#' has fewer than `num_cells` number of cells, all cells from that group will be
#' kept.
#'
#' @family downsampling functions
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr ungroup
#'
#'
tof_downsample_constant <- function(tof_tibble, group_cols = NULL, num_cells) {

  result <-
    tof_tibble %>%
    dplyr::group_by(dplyr::across({{group_cols}})) %>%
    dplyr::filter(
      sample(x = 1:dplyr::n(), size = dplyr::n(), replace = FALSE) %in% 1:num_cells
    ) %>%
    dplyr::ungroup()

  if (inherits(tof_tibble, "tof_tbl")) {
    return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))
  } else {
    return(result)
  }
}


# tof_downsample_prop ----------------------------------------------------------

#' Downsample CyTOF data by randomly selecting a proportion of the cells in each group.
#'
#' This function downsamples the number of cells in a `tof_tbl` by randomly selecting
#' a `prop_cells` proportion of the total number of cells with each unique combination
#' of values in `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Unquoted names of the columns in `tof_tibble` that should
#' be used to define groups from which `prop_cells` will be downsampled.
#' Supports tidyselect helpers. Defaults to `NULL` (no grouping).
#'
#' @param prop_cells A proportion of cells (between 0 and 1) that should be sampled
#' from each group defined by `group_cols`.
#'
#' @return A `tof_tbl` with the same number of columns as the input `tof_tibble`,
#' but fewer rows. Specifically, the number of rows should be `prop_cells` times the
#' number of rows in the input `tof_tibble`.
#'
#' @family downsampling functions
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_sample
#' @importFrom dplyr ungroup
#'
tof_downsample_prop <- function(tof_tibble, group_cols = NULL, prop_cells) {
  result <-
    tof_tibble %>%
    dplyr::group_by(dplyr::across({{group_cols}})) %>%
    dplyr::slice_sample(prop = prop_cells) %>%
    dplyr::ungroup()

  if (inherits(tof_tibble, "tof_tbl")) {
    return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))
  } else {
    return(result)
  }
}


# tof_downsample_density -------------------------------------------------------

#' Downsample CyTOF data by randomly selecting a proportion of the cells in each group.
#'
#' This function downsamples the number of cells in a `tof_tbl` using the
#' density-dependent downsampling algorithm described in
#' \href{https://pubmed.ncbi.nlm.nih.gov/21964415/}{Qiu et al., (2011)}.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Unquoted names of the columns in `tof_tibble` that should
#' be used to define groups within which the downsampling will be performed.
#' Supports tidyselect helpers. Defaults to `NULL` (no grouping).
#'
#' @param density_cols Unquoted names of the columns in `tof_tibble` to use in the
#' density estimation for each cell. Defaults to all numeric columns in `tof_tibble`.
#'
#' @param outlier_percentile The local density percentile (i.e. a value between 0 and 1)
#' below which cells should be considered outliers (and discarded). Cells with a local
#' density below `outlier_percentile` will never be selected during the downsampling
#' proceure. Defaults to 0.01 (cells below the 1st local density percentile will be removed).
#'
#' @param target_percentile The local density percentile (i.e. a value between 0 and 1) to which the downsampling
#' procedure should adjust all cells. In short, the algorithm will continue to remove
#' cells from the input `tof_tibble` until the local densities of all remaining cells
#' is equal to `target_percentile`. Lower values will result in more cells being removed.
#' See \href{https://pubmed.ncbi.nlm.nih.gov/21964415/}{Qiu et al., (2011)}
#' for details. Defaults to 0.1 (the 10th percentile of local densities).
#'
#' @param num_neighbors An integer indicating the number of neighbors to use
#' in the k-nearest neighbor local density estimation. Defaults to 15.
#'
#' @param knn_distance_function A character vector indicating which
#' distance function should be used to compute the k-nearest neighbors of each
#' input cell in `tof_tibble`. Valid options are "euclidean" (the default) and "cosine".
#'
#' @param ... Optional additional arguments to pass to \code{\link[RANN]{nn2}}.
#'
#' @return A `tof_tbl` with the same number of columns as the input `tof_tibble`,
#' but fewer rows. The number of rows will depend on the chosen value of `target_percentile`,
#' with fewer cells selected with lower values of `target_percentile`.
#'
#' @family downsampling functions
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr any_of
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr transmute
#' @importFrom dplyr ungroup
#'
#' @importFrom rlang arg_match
#' @importFrom rlang enquo
#'
#' @importFrom tidyselect eval_select
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @importFrom purrr map
#'
#' @importFrom stats runif
#'
tof_downsample_density <-
  function(
    tof_tibble,
    group_cols = NULL,
    density_cols = where(tof_is_numeric),
    outlier_percentile = 0.01,
    target_percentile = 0.1,
    num_neighbors = 15,
    knn_distance_function = c("euclidean", "cosine"),
    ...
  ) {

    # check knn_distance_function
    knn_distance_function <- rlang::arg_match(knn_distance_function)

    #initialize result
    result <-
      tof_tibble %>%
      dplyr::mutate(cell_id = 1:nrow(tof_tibble))

    # extract group and density column names
    group_names <-
      #rlang::enquo(group_cols) %>%
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) %>%
      names()

    density_names <-
      #rlang::enquo(density_cols) %>%
      tidyselect::eval_select(
        expr = rlang::enquo(density_cols),
        data = tof_tibble
      ) %>%
      names()

    # nest data needed to compute densities for each group

    nested_data <-
      result %>%
      dplyr::select(
        .data$cell_id,
        dplyr::any_of(group_names),
        dplyr::any_of(density_names)
      ) %>%
      tidyr::nest(cell_ids = .data$cell_id, data = {{density_cols}})

    # find knn's for all samples
    knn_results <-
      nested_data %>%
      dplyr::group_by(dplyr::across({{group_cols}})) %>%
      dplyr::transmute(
        dplyr::across(dplyr::any_of(group_names)),
        .data$cell_ids,
        knn =
          purrr::map(
            .x = data,
            .f = tof_find_knn,
            k = num_neighbors,
            distance_function = knn_distance_function
          ),
        neighbor_ids = purrr::map(.x = .data$knn, ~.x$neighbor_ids),
        neighbor_distances = purrr::map(.x = .data$knn, ~.x$neighbor_distances)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$knn)

    # estimate knn densities for each cell
    densities <-
      purrr::map2(
        .x = knn_results$neighbor_ids,
        .y = knn_results$neighbor_distances,
        .f = ~
          tof_knn_density(
            neighbor_ids = .x,
            neighbor_distances = .y,
            method = "mean_distance",
            normalize = FALSE
          )
      )

    chosen_cells <-
      knn_results %>%
      dplyr::select(dplyr::any_of(group_names)) %>%
      dplyr::mutate(
        cell_ids = knn_results$cell_ids,
        densities = densities,
        target_density = purrr::map_dbl(.x = densities, .f = quantile, probs = target_percentile)
      ) %>%
      tidyr::unnest(cols = c(.data$cell_ids, .data$densities)) %>%
      dplyr::group_by(dplyr::across({{group_cols}})) %>%
      dplyr::arrange(.data$densities) %>%
      dplyr::mutate(
        rank = 1:dplyr::n(),
        percentile = rank / dplyr::n()
      ) %>%
      dplyr::mutate(
        sample_prob = (.data$densities / .data$target_density),
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sample_value = stats::runif(n = dplyr::n())) %>%
      dplyr::filter(
        .data$percentile > outlier_percentile,
        .data$sample_prob > .data$sample_value
      ) %>%
      dplyr::pull(.data$cell_id)

    # filter only selected cells out of the original tof_tibble
    result <-
      result %>%
      dplyr::filter(.data$cell_id %in% chosen_cells) %>%
      dplyr::select(-.data$cell_id)

    return(result)
  }

# tof_downsample ---------------------------------------------------------------

#' Downsample CyTOF data.
#'
#' This function downsamples the number of cells in a `tof_tbl` using the
#' one of three methods (randomly sampling a constant number of cells,
#' randomly sampling a proportion of cells, or performing density-dependent
#' downsampling per the algorithm in
#' \href{https://pubmed.ncbi.nlm.nih.gov/21964415/}{Qiu et al., (2011)}).
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Unquoted names of the columns in `tof_tibble` that should
#' be used to define groups within which the downsampling will be performed.
#' Supports tidyselect helpers. Defaults to `NULL` (no grouping).
#'
#'
#' @param ... Additional arguments to pass to the `tof_downsample_*` function
#' family member corresponding to the chosen method.
#'
#' @param method A string indicating which downsampling method to use: "constant"
#' (the default), "prop", or "density".
#'
#' @return A downsampled `tof_tbl` with the same number of columns as the input
#' `tof_tibble`, but fewer rows. The number of rows in the result will depend
#' on the chosen downsampling method.
#'
#' @family downsampling functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#'
tof_downsample <-
  function(
    tof_tibble,
    group_cols = NULL,
    ...,
    method = c("constant", "prop", "density")
  ) {

    # check method argument
    method <-
      rlang::arg_match(arg = method, values = c("constant", "prop", "density"))

    # perform the downsampling
    if (method == "constant") {
      result <-
        tof_tibble %>%
        tof_downsample_constant(group_cols = {{group_cols}}, ...)
    } else if (method == "prop") {
      result <-
        tof_tibble %>%
        tof_downsample_prop(group_cols = {{group_cols}}, ...)
    } else {
      result <-
        tof_tibble %>%
        tof_downsample_density(group_cols = {{group_cols}}, ...)
    }
    return(result)
  }


