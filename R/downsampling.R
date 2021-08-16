# downsampling.R
# This file contains functions relevant to downsampling cells within
# tof_tibble objects containing CyTOF data.

#' Title
#'
#' @param tof_tibble TO DO
#'
#' @param group_cols TO DO
#'
#' @param num_cells TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_downsample_constant <- function(tof_tibble, group_cols, num_cells) {

  result <-
    tof_tibble %>%
    dplyr::group_by({{group_cols}}) %>%
    dplyr::slice_sample(n = num_cells) %>%
    dplyr::ungroup()

  return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))
}


#' Title
#'
#' @param tof_tibble TO DO
#'
#' @param group_cols TO DO
#'
#' @param prop_cells TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_downsample_prop <- function(tof_tibble, group_cols, prop_cells) {
  result <-
    tof_tibble %>%
    dplyr::group_by({{group_cols}}) %>%
    dplyr::slice_sample(prop = prop_cells) %>%
    dplyr::ungroup()

  return(new_tof_tibble(x = result, panel = tof_get_panel(tof_tibble)))

}


tof_downsample_density <-
  function(
    tof_tibble,
    group_cols,
    density_cols = where(tof_is_numeric),
    outlier_percentile = 0.01,
    target_percentile = 0.1,
    num_neighbors = 15,
    knn_distance_function = c("euclidean", "cosine"),
    ...#optional additional arguments for RANN:nn2
  ) {
    # check knn_distance_function
    knn_distance_function <- rlang::arg_match(knn_distance_function)

    #initialize result
    result <-
      tof_tibble %>%
      dplyr::mutate(cell_id = 1:nrow(tof_tibble))

    # nest data needed to compute densities for each group
    nested_data <-
      result %>%
      dplyr::select(cell_id, {{group_cols}}, {{density_cols}}) %>%
      tidyr::nest(cell_ids = cell_id, data = {{density_cols}})

    # find knn's for all samples
    knn_results <-
      nested_data %>%
      dplyr::transmute(
        file_name,
        cell_ids,
        knn =
          purrr::map(
            .x = data,
            .f = tof_find_knn,
            k = num_neighbors,
            distance_function = knn_distance_function
          ),
        neighbor_ids = purrr::map(.x = knn, ~.x$neighbor_ids),
        neighbor_distances = purrr::map(.x = knn, ~.x$neighbor_distances)
      ) %>%
      dplyr::select(-knn)

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
      dplyr::tibble(
        file_name = knn_results$file_name,
        cell_ids = knn_results$cell_ids,
        densities = densities,
        target_density = purrr::map_dbl(.x = densities, .f = quantile, probs = target_percentile)
      ) %>%
      tidyr::unnest(cols = c(cell_ids, densities)) %>%
      dplyr::group_by({{group_cols}}) %>%
      dplyr::arrange(densities) %>%
      dplyr::mutate(
        rank = 1:dplyr::n(),
        percentile = rank / dplyr::n()
      ) %>%
      mutate(
        sample_prob = (densities / target_density),
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sample_value = runif(n = nrow(.))) %>%
      dplyr::filter(percentile > outlier_percentile, sample_prob > sample_value) %>%
      dplyr::pull(cell_id)

    # filter only selected cells out of the original tof_tibble
    result <-
      result %>%
      dplyr::filter(cell_id %in% chosen_cells) %>%
      dplyr::select(-cell_id)

    return(result)
  }

