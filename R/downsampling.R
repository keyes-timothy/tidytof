# downsampling.R
# This file contains functions relevant to downsampling cells within
# tof_tibble objects containing high-dimensional cytometry data.

# tof_downsample_constant ------------------------------------------------------

#' Downsample high-dimensional cytometry data by randomly selecting a constant number of cells per group.
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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE)
#'     )
#'
#' # sample 500 cells from the input data
#' tof_downsample_constant(
#'     tof_tibble = sim_data,
#'     num_cells = 500L
#' )
#'
#' # sample 20 cells per cluster from the input data
#' tof_downsample_constant(
#'     tof_tibble = sim_data,
#'     group_cols = cluster_id,
#'     num_cells = 20L
#' )
#'
#'
tof_downsample_constant <- function(tof_tibble, group_cols = NULL, num_cells) {

  result <-
    tof_tibble |>
    dplyr::group_by(dplyr::across({{group_cols}})) |>
    dplyr::filter(
      sample(x = 1:dplyr::n(), size = dplyr::n(), replace = FALSE) %in% 1:num_cells
    ) |>
    dplyr::ungroup()

  if (inherits(tof_tibble, "tof_tbl")) {
    return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))
  } else {
    return(result)
  }
}



# tof_downsample_prop ----------------------------------------------------------

#' Downsample high-dimensional cytometry data by randomly selecting a proportion of the cells in each group.
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
#'@examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE)
#'     )
#'
#' # sample 10% of all cells from the input data
#' tof_downsample_prop(
#'     tof_tibble = sim_data,
#'     prop_cells = 0.1
#' )
#'
#' # sample 10% of all cells from each cluster in the input data
#' tof_downsample_prop(
#'     tof_tibble = sim_data,
#'     group_cols = cluster_id,
#'     prop_cells = 0.1
#' )
#'
tof_downsample_prop <- function(tof_tibble, group_cols = NULL, prop_cells) {
  result <-
    tof_tibble |>
    dplyr::group_by(dplyr::across({{group_cols}})) |>
    dplyr::slice_sample(prop = prop_cells) |>
    dplyr::ungroup()

  if (inherits(tof_tibble, "tof_tbl")) {
    return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))
  } else {
    return(result)
  }
}


# tof_downsample_density -------------------------------------------------------

#' Downsample high-dimensional cytometry data by randomly selecting a proportion of the cells in each group.
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
#' @param target_num_cells An approximate constant number of cells (between 0 and 1)
#' that should be sampled from each group defined by `group_cols`. Slightly more
#' or fewer cells may be returned due to how the density calculation is performed.
#'
#' @param target_prop_cells An approximate proportion of cells (between 0 and 1)
#' that should be sampled from each group defined by `group_cols`. Slightly more
#' or fewer cells may be returned due to how the density calculation is performed.
#' Ignored if `target_num_cells` is specified.
#'
#' @param target_percentile The local density percentile (i.e. a value between 0 and 1) to which the downsampling
#' procedure should adjust all cells. In short, the algorithm will continue to remove
#' cells from the input `tof_tibble` until the local densities of all remaining cells
#' is equal to `target_percentile`. Lower values will result in more cells being removed.
#' See \href{https://pubmed.ncbi.nlm.nih.gov/21964415/}{Qiu et al., (2011)}
#' for details. Defaults to 0.1 (the 10th percentile of local densities). Ignored
#' if either `target_num_cells` or `target_prop_cells` are specified.
#'
#' @param outlier_percentile The local density percentile (i.e. a value between 0 and 1)
#' below which cells should be considered outliers (and discarded). Cells with a local
#' density below `outlier_percentile` will never be selected during the downsampling
#' procedure. Defaults to 0.01 (cells below the 1st local density percentile will be removed).
#'
#' @param distance_function A string indicating which distance function to use
#' for the cell-to-cell distance calculations. Options include "euclidean"
#' (the default) and "cosine" distances.
#'
#' @param density_estimation_method A string indicating which algorithm should be
#' used to calculate the local density estimate for each cell. Options include
#' k-nearest neighbor density estimation using the mean distance to a cell's
#' k-nearest neighbors ("mean_distance"; the default), k-nearest neighbor density
#' estimation using the summed distance to a cell's k nearest neighbors
#' ("sum_distance") and counting the number of neighboring cells within a
#' spherical radius around each cell as described in Qiu et al., 2011  ("spade").
#' While "spade" often produces the best results, it is slower than
#' knn-density estimation methods.
#'
#' @param ... Optional additional arguments to pass to
#' \code{\link{tof_knn_density}} or \code{\link{tof_spade_density}}.
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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#'
#' tof_downsample_density(
#'     tof_tibble = sim_data,
#'     density_cols = c(cd45, cd34, cd38),
#'     target_prop_cells = 0.5,
#'     density_estimation_method = "spade"
#' )
#'
#' tof_downsample_density(
#'     tof_tibble = sim_data,
#'     density_cols = c(cd45, cd34, cd38),
#'     target_num_cells = 200L,
#'     density_estimation_method = "spade"
#' )
#'
#' tof_downsample_density(
#'     tof_tibble = sim_data,
#'     density_cols = c(cd45, cd34, cd38),
#'     target_num_cells = 200L,
#'     density_estimation_method = "mean_distance"
#' )
#'
#'
tof_downsample_density <-
  function(
    tof_tibble,
    group_cols = NULL,
    density_cols = where(tof_is_numeric),
    target_num_cells,
    target_prop_cells,
    target_percentile = 0.03,
    outlier_percentile = 0.01,
    distance_function = c("euclidean", "cosine", "l2", "ip"),
    density_estimation_method = c("mean_distance", "sum_distance", "spade"),
    ...
  ) {

    # check distance_function
    distance_function <- rlang::arg_match(distance_function)

    #initialize result
    result <-
      tof_tibble |>
      dplyr::mutate(..cell_id = 1:nrow(tof_tibble))

    # extract group and density column names
    group_names <-
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) |>
      names()

    density_names <-
      tidyselect::eval_select(
        expr = rlang::enquo(density_cols),
        data = tof_tibble
      ) |>
      names()

    # nest data needed to compute densities for each group
    nested_data <-
      result |>
      dplyr::select(
        "..cell_id",
        dplyr::any_of(group_names),
        dplyr::any_of(density_names)
      ) |>
      tidyr::nest(cell_ids = "..cell_id", data = {{density_cols}})

    # find local density estimates for each cell in all groups
    nested_data <-
      nested_data |>
      dplyr::mutate(
        densities =
          purrr::map(
            .x = .data$data,
            .f = tof_estimate_density,
            method = density_estimation_method,
            augment = FALSE,
            ...
          )
      )

    # save the local density estimates as one vector per group
    densities <- purrr::map(.x = nested_data$densities, .f = ~.x[[1]])
    nested_data <-
      nested_data |>
      dplyr::select("cell_ids", dplyr::any_of(group_names)) |>
      dplyr::mutate(
        densities = densities
      )

    # if target_num_cells and target_prop_cells are not specified, directly use
    # target_percentile to do the downsampling
    if (missing(target_prop_cells) & missing(target_num_cells)) {

      # store a vector of ..cell_ids for which cells to keep from tof_tibble
      chosen_cells <-
        nested_data |>
        dplyr::mutate(
          target_density =
            purrr::map_dbl(
              .x = densities,
              # protect against returning 0 cells when there are many degenerate
              # densities by computing the target percentile after removing
              # the cells with 0 density (which will be removed during the
              # final filtering step anyway)
              .f = ~ quantile(subset(.x, .x > 0), probs = target_percentile),
            )
        ) |>
        tidyr::unnest(cols = c("cell_ids", "densities")) |>
        dplyr::group_by(dplyr::across({{group_cols}})) |>
        dplyr::arrange(.data$densities) |>
        dplyr::mutate(
          rank = 1:dplyr::n(),
          percentile = rank / dplyr::n()
        ) |>
        dplyr::mutate(
          sample_prob = (.data$target_density / .data$densities)
        ) |>
        dplyr::ungroup() |>
        dplyr::mutate(sample_value = stats::runif(n = dplyr::n())) |>
        dplyr::filter(
          .data$percentile > outlier_percentile,
          .data$sample_prob > .data$sample_value
        ) |>
        dplyr::pull(.data$..cell_id)

      # if either target_num_cells or target_prop_cells are specified, find
      # a threshold that approximates the requested number or proportion of
      # cells
    } else {
      # if target_num_cells is not constant, compute it for each group
      # using target_prop_cells
      if (missing(target_num_cells)) {
        target_num_cells_vector <-
          purrr::map_int(
            .x = nested_data$cell_ids,
            .f =
              ~ as.integer(ceiling(nrow(.x) * target_prop_cells))
          )
      # if target_num_cells is specified, shoot for that constant number
      # in each group
      } else {
        target_num_cells_vector <-
          rep(x = target_num_cells, times = nrow(nested_data))
      }

      chosen_cells <-
        nested_data |>
        dplyr::mutate(
          num_cells = target_num_cells_vector,
          # remove cells at a lower percentile than the outlier percentile
          # and also any cells whose local densities are exactly 0
          cells_to_remove =
            purrr::map(
              .x = .data$densities,
              .f = ~
                dplyr::tibble(
                  densities = .x,
                  ..cell_ids = 1:length(.x)
                ) |>
                dplyr::arrange(.data$densities) |>
                dplyr::mutate(
                  rank = 1:dplyr::n(),
                  percentile = rank / dplyr::n(),
                ) |>
                dplyr::filter(
                  .data$percentile <= outlier_percentile | .data$densities == 0
                ) |>
                dplyr::pull(.data$..cell_ids)
            ),
          densities =
            purrr::map2(
              .x = .data$densities,
              .y = .data$cells_to_remove,
              .f = ~ .x[-.y]
            ),
          cell_ids =
            purrr::map2(
              .x = .data$cell_ids,
              .y = .data$cells_to_remove,
              .f = ~ .x[-.y, ]
            )
        ) |>
        # use tof_spade_downsampling to find which cells should be retained
        # after downsampling
        dplyr::transmute(
          sampled_cells =
            purrr::pmap(
              .l = list(.data$densities, .data$cell_ids, .data$num_cells),
              .f = tof_spade_downsampling
            )
        ) |>
        dplyr::pull(.data$sampled_cells) |>
        c(recursive = TRUE)

    }

    # filter only selected cells out of the original tof_tibble
    result <-
      result |>
      dplyr::filter(.data$..cell_id %in% chosen_cells) |>
      dplyr::select(-"..cell_id")

    return(result)
  }

# tof_downsample ---------------------------------------------------------------

#' Downsample high-dimensional cytometry data.
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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE)
#'     )
#'
#' # sample 200 cells from the input data
#' tof_downsample(
#'     tof_tibble = sim_data,
#'     num_cells = 200L,
#'     method = "constant"
#' )
#'
#' # sample 10% of all cells from the input data
#' tof_downsample(
#'     tof_tibble = sim_data,
#'     prop_cells = 0.1,
#'     method = "prop"
#' )
#'
#' # sample ~10% of cells from the input data using density dependence
#' tof_downsample(
#'     tof_tibble = sim_data,
#'     target_prop_cells = 0.1,
#'     method = "density"
#' )
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
        tof_tibble |>
        tof_downsample_constant(group_cols = {{group_cols}}, ...)
    } else if (method == "prop") {
      result <-
        tof_tibble |>
        tof_downsample_prop(group_cols = {{group_cols}}, ...)
    } else {
      result <-
        tof_tibble |>
        tof_downsample_density(group_cols = {{group_cols}}, ...)
    }
    return(result)
  }


