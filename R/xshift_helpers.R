# xshift_helpers.R
# This file contains helper functions for performing X-Shift clustering
# per Samusik et al 2016 (Nature Methods).
#
# All functions in this file are subroutines of the function tof_cluster_xshift
# in clustering.R


# step 1
xshift_compute_local_densities <-
  function(
    neighbor_ids,
    neighbor_distances,
    distance_function = c("cosine", "euclidean")
  ) {

    # check distance_function argument
    distance_function <-
      match.arg(distance_function, choices = c("cosine", "euclidean"))

    if (distance_function == "cosine") {

      # convert distances from l2-normalized-euclidean to angular
      neighbor_cosine_distances <-
        # first convert to cosine distance
        (neighbor_distances^2) / 2

      neighbor_cosine_similarities <-
        1 - neighbor_cosine_distances

      # very little difference between using this and using cosine distance...
      # correlation in the final densities is very high (>0.99)
      angular_distances <-
        acos(neighbor_cosine_similarities)

      neighbor_distances <- angular_distances
    }

  # compute densities
  densities <-
    tof_knn_density(
      neighbor_ids = neighbor_ids,
      neighbor_distances = neighbor_distances,
      method = "sum_distance",
      normalize = FALSE
    )

  return(densities)
}

# step 2
xshift_find_candidate_centroids <- function(neighbor_ids, densities) {
  # check each datapoint for a neighbor with a higher density value
  # and connect the two points if such neighbor exists

  # find which neighbors of a given cell have a higher density estimate
  # than the cell itself and store in the `higher_density_neighbors` column
  candidate_centroids <-
    dplyr::tibble(
      cell_id = 1:length(densities),
      higher_density_neighbors =
        purrr::map(
          .x = cell_id,
          .f = ~return(neighbor_ids[.x, (densities[neighbor_ids[.x,]] > densities[.x])])
        ),
      closest_neighbor =
        purrr::map_int(
          .x = higher_density_neighbors,
          .f = ~ purrr::pluck(.x, 1, .default = NA_integer_)
        ),
      length = purrr::map_int(.x = higher_density_neighbors, length),
      is_candidate = dplyr::if_else(length == 0, TRUE, FALSE)
    ) %>%
    dplyr::select(-length)

  # returns a tibble with 3 columns
  # cell_id: integer id of each cell (row in the original dataset)
  # higher_density_neighbors: a list of integer vectors representing the neighbors
  #                           of each cell with a higher density than the cell itself
  # is_candidate: boolean indicating if a cell is a candidate centroid (no z-nearest
  #               neighbors have a higher density than the cell itself)
  return(candidate_centroids)
}

# step 3
xshift_find_midpoint_distances <-
  function(
    midpoint, # numeric vector representing the midpoint's location in high-dimensional space
    centroids # tibble in which each row represents a candidate centroid (with ids given by ..cell_id)
  ) {
    my_centroids <-
      centroids %>%
      dplyr::select(-..cell_id) %>%
      t()

    # find euclidean distance from the midpoint to all centroids (because
    # `centroids` has already been l2-normalized if cosine distance is being
    # used, this generalizes to all cases)
    distances <-
      ((my_centroids - midpoint)^2) %>%
      base::colSums() %>%
      base::sqrt() %>%
      stats::setNames(nm = centroids$..cell_id)

    # return the result as a tibble
    result <-
      dplyr::tibble(
        centroid = names(distances),
        distance = distances
      )

    return(result)
  }

xshift_find_close_centroids <-
  function(distances_to_all_centroids, distance_to_midpoint) {
    result <-
      distances_to_all_centroids %>%
      dplyr::filter(
        base::round(distance, 5) < base::round(distance_to_midpoint, 5)
      ) %>%
      dplyr::pull(centroid)

    return(result)
  }

xshift_find_gabriel_neighbors <-
  function(candidates, distance_function) {
    candidate_cells <-
      candidates %>%
      dplyr::filter(is_candidate) %>%
      dplyr::pull(cell_id)

    centroids <- tof_tibble[candidate_cells, ]

    if(distance_function == "cosine") {
      centroids <-
        t(apply(X = centroids, MARGIN = 1, FUN = l2_normalize)) %>%
        dplyr::as_tibble()
    }

    centroids <-
      centroids %>%
      dplyr::mutate(..cell_id = candidate_cells)

    centroid_distances <-
      centroids %>%
      dplyr::select(-..cell_id) %>%
      stats::dist() %>%
      as.matrix()

    base::row.names(centroid_distances) <- candidate_cells
    base::colnames(centroid_distances) <- candidate_cells

    gabriel_tibble <-
      utils::combn(x = candidate_cells, m = 2) %>%
      base::t() %>%
      dplyr::as_tibble(.name_repair = "minimal")

    colnames(gabriel_tibble) <- c("centroid_1", "centroid_2")

    gabriel_tibble <-
      gabriel_tibble %>%
      dplyr::mutate(
        midpoint =
          purrr::map2(
            .x = centroid_1,
            .y = centroid_2,
            .f = ~
              ((centroids[centroids$..cell_id == .x, ] + centroids[centroids$..cell_id == .y, ]) / 2) %>%
              dplyr::select(-..cell_id) %>%
              base::as.numeric()
          ),
        distance_to_midpoint =
          purrr::map2_dbl(
            .x = base::as.character(centroid_1),
            .y = base::as.character(centroid_2),
            .f = ~ centroid_distances[.x, .y] / 2
          ),
        distances_to_all_centroids =
          purrr::map(
            .x = midpoint,
            .f = xshift_find_midpoint_distances,
            centroids = centroids
          ),
        close_clusters =
          purrr::map2(
            .x = distances_to_all_centroids,
            .y = distance_to_midpoint,
            .f = xshift_find_close_centroids
          ),
        are_gabriel = purrr::map_lgl(.x = close_clusters, .f = ~length(.x) == 0)
      )

    gabriel_pairs <-
      gabriel_tibble %>%
      dplyr::filter(
        are_gabriel
      ) %>%
      dplyr::select(centroid_1, centroid_2)

    result <- list(centroids = centroids, gabriel_pairs = gabriel_pairs)

    return(result)
  }

xshift_make_path <-
  function(
    centroid_1, # integer indicating which centroid is first
    centroid_2, # integer indicating which centroid is second
    centroids # tibble of all protein measurements for each centroid
  ) {
    centroid_1_vector <-
      centroids[centroids$..cell_id == centroid_1, colnames(centroids) != "..cell_id"] %>%
      as.numeric() %>%
      stats::setNames(colnames(centroids)[colnames(centroids) != "..cell_id"])

    centroid_2_vector <-
      centroids[centroids$..cell_id == centroid_2, colnames(centroids) != "..cell_id"] %>%
      as.numeric() %>%
      stats::setNames(colnames(centroids)[colnames(centroids) != "..cell_id"])

    result <-
      purrr::map_dfr(
        .x = seq(0, 1, 0.1),
        .f = ~ (centroid_1_vector * .x) + (centroid_2_vector * (1 - .x))
      )

    return(result)
  }


xshift_test_density_minima <- function(gabriel_pairs, centroids, k, tof_tibble) {
  result <-
    gabriel_pairs %>%
    dplyr::mutate(
      path_cells =
        purrr::map2(
          .x = centroid_1,
          .y = centroid_2,
          .f = ~ xshift_make_path(centroid_1 = .x, centroid_2 = .y, centroids = centroids)
        ),
      nn_results =
        purrr::map(
          .x = path_cells,
          .f = ~
            tof_find_knn(
              .data = tof_tibble,
              k = k,
              distance_function = "euclidean",
              query = .x,
              eps = 0.5
            )
        ),
      nn_neighbors = purrr::map(.x = nn_results, .f = ~.x$neighbor_ids),
      nn_distances = purrr::map(.x = nn_results, .f = ~.x$neighbor_distances),
      densities =
        purrr::map2(
          .x = nn_neighbors,
          .y = nn_distances,
          .f = tof_knn_density,
          method = "sum_distance",
          normalize = FALSE
        ),
      min_density_index = purrr::map_int(.x = densities, .f = ~which(.x == min(.x)))
    ) %>%
    dplyr::transmute(
      centroid_1,
      centroid_2,
      should_merge = min_density_index %in% c(1, 11),
      merge_into =
        dplyr::case_when(
          min_density_index == 1 ~ centroid_2,
          min_density_index == 11 ~ centroid_1,
          TRUE                   ~ NA_integer_
        )
    )

  return(result)
}

xshift_merge_centroids <- function(candidates, centroids, distance_function) {
  # Determine which candidate centroids are Gabriel neighbors - aka do not have
  # another centroid closer to their midpoint
  gabriel_pairs <-
    xshift_find_gabriel_neighbors(candidates, distance_function)


}

# step 4
xshift_initialize_assignments <- function(.data) {
  # TO DO
  return(NULL)
}

# step 5
xshift_mahalanobis_merge <- function(.data) {
  # TO DO
  return(NULL)
}

xshift_finalize_assignments <- function(.data) {
  # TO DO
  return(NULL)
}



