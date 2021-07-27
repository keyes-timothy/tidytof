# xshift_helpers.R
# This file contains helper functions for performing X-Shift clustering
# per Samusik et al 2016 (Nature Methods).
#
# All functions in this file are subroutines of the function tof_cluster_xshift
# in clustering.R


# step 1
xshift_find_knn <- function(.data, k, distance_function = c("euclidean", "cosine")) {
  # TO DO

  # RANN::nn2 for euclidean distance
  # l2 normalize all rows and then RANN::nn2 for cosine distance

  return(NULL)
}

xshift_compute_local_densities <- function(.data) {
  # TO DO
  return(NULL)
}

# step 2
xshift_find_candidate_centroids <- function(.data) {
  # TO DO
}

# step 3
xshift_merge_centroids <- function(.data) {
  # TO DO
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



