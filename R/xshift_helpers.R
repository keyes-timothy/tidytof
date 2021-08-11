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
      method = "sum_distance"
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
    tibble::tibble(
      cell_id = 1:length(densities),
      higher_density_neighbors =
        purrr::map(
          .x = cell_id,
          .f = ~return(neighbor_ids[.x, (densities[neighbor_ids[.x,]] > densities[.x])])
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
xshift_merge_centroids <- function(candidate_centroids) {
  # Determine which candidate centroids are Gabriel neighbors - aka have another
  # centroid closer to their midpoint
  return(NULL)
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



