# upsample.R
# This file contains functions relevant to performing upsampling
# on tof_tbl objects containing CyTOF data.

# tof_upsample_distance --------------------------------------------------------

#' Upsample cells into the closest cluster in a reference dataset
#'
#' This function performs distance-based upsampling on CyTOF data
#' by sorting single cells (passed into the function as `tof_tibble`) into
#' their most phenotypically similar cell subpopulation in a reference dataset
#' (passed into the function as `reference_tibble`). It does so by calculating
#' the distance (either mahalanobis, cosine, or pearson) between each cell in
#' `tof_tibble` and the centroid of each cluster in `reference_tibble`, then
#' sorting cells into the cluster corresponding to their closest centroid.
#'
#' @param tof_tibble A `tibble` or `tof_tbl` containing cells to be upsampled
#' into their nearest reference subpopulation.
#'
#' @param reference_tibble A `tibble` or `tof_tibble` containing cells that have
#' already been clustered or manually gated into subpopulations.
#'
#' @param reference_cluster_col An unquoted column name indicating which column in
#' `reference_tibble` contains the subpopulation label (or cluster id) for
#' each cell in `reference_tibble`.
#'
#' @param upsample_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the distances used for upsampling. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param parallel_cols Optional. Unquoted column names indicating which columns in `tof_tibble` to
#' use for breaking up the data in order to parallelize the upsampling using
#' `foreach` on a `doParallel` backend.
#' Supports tidyselect helpers.
#'
#' @param distance_function A string indicating which distance function should
#' be used to perform the upsampling. Options are "mahalanobis" (the default),
#' "cosine", and "pearson".
#'
#' @param num_cores An integer indicating the number of CPU cores used to parallelize
#' the classification. Defaults to 1 (a single core).
#'
#' @param return_distances A boolean value indicating whether or not the returned
#' result should include only one column, the cluster ids corresponding to each row
#' of `tof_tibble` (return_distances = FALSE, the default), or if the returned
#' result should include additional columns representing the distance between each
#' row of `tof_tibble` and each of the reference subpopulation centroids
#' (return_distances = TRUE).
#'
#' @return If `return_distances = FALSE`, a tibble with one column named
#' `.upsample_cluster`, a character vector of length `nrow(tof_tibble)`
#' indicating the id of the reference cluster to which each cell
#' (i.e. each row) in `tof_tibble` was assigned.
#'
#' If `return_distances = TRUE`, a tibble with `nrow(tof_tibble)` rows and num_clusters + 1
#' columns, where num_clusters is the number of clusters in `reference_tibble`.
#' Each row represents a cell from `tof_tibble`, and num_clusters
#' of the columns represent the distance between the cell and each of the reference
#' subpopulations' cluster centroids. The final column represents the cluster id of
#' the reference subpopulation with the minimum distance to the cell represented
#' by that row.
#'
#' @export
#'
#'
tof_upsample_distance <-
  function(
    tof_tibble,
    reference_tibble,
    reference_cluster_col,
    upsample_cols = where(tof_is_numeric),
    parallel_cols,
    distance_function = c("mahalanobis", "cosine", "pearson"),
    num_cores = 1L,
    return_distances = FALSE
  ) {

    # if computed on 1 core
    if(missing(parallel_cols)) {
      result <-
        tof_cluster_ddpr(
          tof_tibble = tof_tibble,
          healthy_tibble = reference_tibble,
          healthy_label_col = {{reference_cluster_col}},
          cluster_cols = {{upsample_cols}},
          distance_function = distance_function,
          num_cores = num_cores,
          return_distances = return_distances,
          verbose = FALSE
        )

    # if computed in parallel
    } else {
      result <-
        tof_cluster_ddpr(
          tof_tibble = tof_tibble,
          healthy_tibble = reference_tibble,
          healthy_label_col = {{reference_cluster_col}},
          cluster_cols = {{upsample_cols}},
          parallel_cols = {{parallel_cols}},
          distance_function = distance_function,
          num_cores = num_cores,
          return_distances = return_distances,
          verbose = FALSE
        )
    }
    result_colnames <- colnames(result)
    cluster_colname <- result_colnames[grepl("_cluster$", result_colnames)]
    upsample_clusters <- result[[cluster_colname]]

    new_result <- result[, !(result_colnames %in% cluster_colname)]
    new_result$`.upsample_cluster` <- upsample_clusters
    result <- new_result
    return(result)
  }

# tof_upsample_neighbor --------------------------------------------------------

#' Upsample cells into the cluster of their nearest neighbor a reference dataset
#'
#' This function performs upsampling on CyTOF data
#' by sorting single cells (passed into the function as `tof_tibble`) into
#' their most phenotypically similar cell subpopulation in a reference dataset
#' (passed into the function as `reference_tibble`). It does so by finding
#' each cell in `tof_tibble`'s nearest neighbor in `reference_tibble` and assigning
#' it to the cluster to which its nearest neighbor belongs. The nearest neighbor
#' calculation can be performed with either euclidean or cosine distance.
#'
#' @param tof_tibble A `tibble` or `tof_tbl` containing cells to be upsampled
#' into their nearest reference subpopulation.
#'
#' @param reference_tibble A `tibble` or `tof_tibble` containing cells that have
#' already been clustered or manually gated into subpopulations.
#'
#' @param reference_cluster_col An unquoted column name indicating which column in
#' `reference_tibble` contains the subpopulation label (or cluster id) for
#' each cell in `reference_tibble`.
#'
#' @param upsample_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the distances used for upsampling. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param distance_function A string indicating which distance function should
#' be used to perform the upsampling. Options are "euclidean" (the default) and
#' "cosine".
#'
#' @return A tibble with one column named
#' `.upsample_cluster`, a character vector of length `nrow(tof_tibble)`
#' indicating the id of the reference cluster to which each cell
#' (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr tibble
#'
tof_upsample_neighbor <-
  function(
    tof_tibble,
    reference_tibble,
    reference_cluster_col,
    upsample_cols = where(tof_is_numeric),
    distance_function = c("euclidean", "cosine")
    ) {
    query_matrix <-
      tof_tibble %>%
      dplyr::select({{upsample_cols}}) %>%
      as.matrix()

    nn_result <-
      reference_tibble %>%
      dplyr::select({{upsample_cols}}) %>%
      tof_find_knn(
        k = 1,
        distance_function = distance_function,
        query = query_matrix
      )

    nn_ids <- nn_result$neighbor_ids

    upsampled_clusters <-
      dplyr::pull(reference_tibble, {{reference_cluster_col}})[nn_ids]

    result <-
      dplyr::tibble(.upsample_cluster = upsampled_clusters)

    return(result)
  }

# tof_upsample -----------------------------------------------------------------

#' Upsample cells into the closest cluster in a reference dataset
#'
#' This function performs distance-based upsampling on CyTOF data
#' by sorting single cells (passed into the function as `tof_tibble`) into
#' their most phenotypically similar cell subpopulation in a reference dataset
#' (passed into the function as `reference_tibble`). It does so by calculating
#' the distance (either mahalanobis, cosine, or pearson) between each cell in
#' `tof_tibble` and the centroid of each cluster in `reference_tibble`, then
#' sorting cells into the cluster corresponding to their closest centroid.
#'
#' @param tof_tibble A `tibble` or `tof_tbl` containing cells to be upsampled
#' into their nearest reference subpopulation.
#'
#' @param reference_tibble A `tibble` or `tof_tibble` containing cells that have
#' already been clustered or manually gated into subpopulations.
#'
#' @param reference_cluster_col An unquoted column name indicating which column in
#' `reference_tibble` contains the subpopulation label (or cluster id) for
#' each cell in `reference_tibble`.
#'
#' @param upsample_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the distances used for upsampling. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param ... Additional arguments to pass to the `tof_upsample_*`
#' function family member corresponding to the chosen method.
#'
#' @param add_col A boolean value indicating if the output should column-bind the
#' cluster ids of each cell as a new column in `tof_tibble` (TRUE, the default) or if
#' a single-column tibble including only the cluster ids should be returned (FALSE).
#'
#' @param method A string indicating which clustering methods should be used. Valid
#' values include "distance" (default) and "neighbor".
#'
#' @return A `tof_tbl` or `tibble` If add_col = FALSE, it will have a single column encoding
#' the upsampled cluster ids for each cell in `tof_tibble`.
#' If add_col = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the cluster ids.
#'
#' @export
#'
#' @importFrom dplyr bind_cols
#'
tof_upsample <-
  function(
    tof_tibble,
    reference_tibble,
    reference_cluster_col,
    upsample_cols,
    ...,
    add_col = TRUE,
    method = c("distance", "neighbor")
  ) {

    method <- rlang::arg_match(arg = method, values = c("distance", "neighbor"))

    if (method == "distance") {
      result <-
        tof_upsample_distance(
          tof_tibble = tof_tibble,
          reference_tibble = reference_tibble,
          reference_cluster_col = {{reference_cluster_col}},
          upsample_cols = {{upsample_cols}},
          ...
        )
    } else if (method == "neighbor") {
      result <-
        tof_upsample_neighbor(
          tof_tibble = tof_tibble,
          reference_tibble = reference_tibble,
          reference_cluster_col = {{reference_cluster_col}},
          upsample_cols = {{upsample_cols}},
          ...
        )
    } else {
      stop("Not a valid method.")
    }

    if (add_col) {
      result <-
        dplyr::bind_cols(tof_tibble, result)
    }
    return(result)
  }
