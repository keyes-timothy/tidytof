# clustering.R
# This file contains functions relevant to performing single-cell clustering
# on tof_tibble objects containing CyTOF data.

# tof_cluster_flowsom -----------------

#' Perform flowSOM clustering on CyTOF data
#'
#' This function performs FlowSOM clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[FlowSOM]{BuildSOM}} and \code{\link[FlowSOM]{MetaClustering}}.
#'
#' For additional details about the FlowSOM algorithm,
#' see \href{https://pubmed.ncbi.nlm.nih.gov/25573116/}{this paper}.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param cluster_vars Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the flowSOM clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param som_xdim The width of the grid used by the self-organizing map. The
#' total number of clusters returned by FlowSOM will be som_xdim * som_ydim,
#' so adjust this value to affect the final number of clusters. Defaults to 10.
#'
#' @param som_ydim The height of the grid used by the self-organizing map. The
#' total number of clusters returned by FlowSOM will be som_xdim * som_ydim,
#' so adjust this value to affect the final number of clusters. Defaults to 10.
#'
#' @param som_distance_function The distance function used during self-organizing
#' map calculations. Options are "euclidian" (the default), "manhattan", "chebyshev",
#' and "cosine".
#'
#' @param perform_metaclustering A boolean value indicating if metaclustering
#' should be performed on the initial clustering result returned by FlowSOM.
#' Defaults to TRUE.
#'
#' @param num_metaclusters An integer indicating the maximum number of metaclusters
#' that should be returned after metaclustering. Defaults to 20.
#'
#' @param ... Optional additional parameters that can be passed to the \code{\link[FlowSOM]{BuildSOM}}
#' function.
#'
#' @return An integer vector of length `nrow(tof_tibble)` indicating the id of
#' the flowSOM cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @examples
#'
#' tof_tibble <- tof_read_data(tidytof_example_data("phenograph")[[1]])
#'
#' flowsom_clusters <-
#'    tof_cluster_flowsom(
#'      tof_tibble,
#'      cluster_vars = contains("CD", ignore.case = FALSE)
#'    )
#'
#'
#'
tof_cluster_flowsom <-
  function(
    tof_tibble = NULL,
    cluster_vars = where(tof_is_numeric),
    som_xdim = 10,
    som_ydim = 10,
    som_distance_function = c("euclidean", "manhattan", "chebyshev", "cosine"),
    perform_metaclustering = TRUE,
    num_metaclusters = 20,
    ...
  ) {
    som_distance_function <-
      match.arg(
        arg = som_distance_function,
        choices = c("euclidean", "manhattan", "chebyshev", "cosine")
      )

    # extract string indicating which markers should be used for clustering
    clustering_markers <-
      tof_tibble %>%
      dplyr::select({{cluster_vars}}) %>%
      colnames()

    # build the flowsom object
    fsom <-
      list(
        data =
          tof_tibble %>%
          select(all_of(clustering_markers)) %>%
          data.matrix(),
        compensate = FALSE,
        spillover = NULL,
        transform = FALSE,
        scale = NULL,
        prettyColnames = clustering_markers
      )

    # build self-organizing map and extract cluster labels
    distf <-
      switch(
        som_distance_function,
        manhattan = 1,
        euclidean = 2,
        chebyshev = 3,
        cosine = 4
      )

    som <-
      FlowSOM::BuildSOM(
        fsom = fsom,
        colsToUse = clustering_markers,
        silent = TRUE,
        xdim = som_xdim,
        ydim = som_ydim,
        distf = distf,
        ...
      )

    # if no metaclustering, return flowSOM cluster labels
    if (!perform_metaclustering) {
      flowsom_clusters <- som$map$mapping[,1]
      return(tibble::tibble(flowsom_cluster = flowsom_clusters))

    # otherwise, perform metaclustering
    } else {

      mst <- FlowSOM::BuildMST(som, silent = TRUE, tSNE = FALSE)

      flowsom_metacluster_object <-
        FlowSOM::MetaClustering(
          data = mst$map$codes,
          method = "metaClustering_consensus",
          max = num_metaclusters
        )

      flowsom_metaclusters <-
        flowsom_metacluster_object[mst$map$mapping[,1]] %>%
        as.integer()

      return(tibble::tibble(flowsom_metacluster = flowsom_metaclusters))
    }
  }


# tof_cluster_phenograph --------------

#' Perform Phenograph clustering on CyTOF data.
#'
#'This function performs Phenograph clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[Rphenograph]{Rphenograph}}.
#'
#' For additional details about the Phenograph algorithm,
#' see \href{https://pubmed.ncbi.nlm.nih.gov/25573116/}{this paper}.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param cluster_vars Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the flowSOM clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_neighbors An integer indicating the number of neighbors to use when
#' constructing Phenograph's k-nearest-neighbor graph. Smaller values emphasize
#' local graph structure; larger values emphasize global graph structure (and
#' will add time to the computation). Defaults to 30.
#'
#' @param ... Optional additional parameters that can be passed to \code{\link[Rphenograph]{Rphenograph}}.
#'
#' @return An integer vector of length `nrow(tof_tibble)` indicating the id of
#' the Phenograph cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @examples
#' tof_tibble <- tof_read_data(tidytof_example_data("phenograph")[[1]])
#'
#' phenograph_clusters <-
#'    tof_cluster_phenograph(
#'      tof_tibble,
#'      cluster_vars = contains("CD", ignore.case = FALSE)
#'    )
#'
#'
tof_cluster_phenograph <-
  function(
    tof_tibble,
    cluster_vars = where(tof_is_numeric),
    num_neighbors = 30,
    ...
  ) {
    invisible(
      utils::capture.output(
        suppressMessages(
          phenograph_result <-
            Rphenograph::Rphenograph(
              data = dplyr::select(tof_tibble, {{cluster_vars}}),
              k = num_neighbors
            )
        )
      )
    )
    phenograph_clusters <- igraph::membership(phenograph_result[[2]])

    # deal with cells not assigned to any cluster
    if(length(phenograph_clusters) != nrow(tof_tibble)) {
      # assign all cells an id
      cell_ids <- as.character(as.numeric(1:nrow(tof_tibble)))

      # find the cells not assigned to a cluster and assign them to NA
      unassigned_cells <- which(!(cell_ids %in% names(phenograph_clusters)))
      unassigned_cell_clusters <- rep(NA, length(unassigned_cells))
      names(unassigned_cell_clusters) <- as.character(unassigned_cells)

      # concatenate
      phenograph_clusters <- c(phenograph_clusters, unassigned_cell_clusters)
      phenograph_clusters <-
        phenograph_clusters[order(as.numeric(names(phenograph_clusters)))]
    }

    #return final result
    phenograph_clusters <- as.integer(phenograph_clusters)
    return(tibble(phenograph_cluster = phenograph_clusters))
  }


# tof_cluster_kmeans ------------------

#' Perform k-means clustering on CyTOF data.
#'
#' This function performs k-means clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[stats]{kmeans}}.
#'
#' @param tof_tibbleA `tof_tibble`.
#'
#' @param cluster_vars Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the flowSOM clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_clusters An integer indicating the maximum number of clusters
#' that should be returned Defaults to 20.
#'
#' @param ... Optional additional arguments that can be passed to
#' \code{\link[stats]{kmeans}}.
#'
#' @return An integer vector of length `nrow(tof_tibble)` indicating the id of
#' the k-means cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @examples
#' tof_tibble <- tof_read_data(tidytof_example_data("phenograph")[[1]])
#'
#' kmeans_clusters <-
#'    tof_cluster_kmeans(
#'      tof_tibble,
#'      cluster_vars = contains("CD", ignore.case = FALSE)
#'    )
#'
tof_cluster_kmeans <-
  function(
    tof_tibble,
    cluster_vars = where(tof_is_numeric),
    num_clusters = 20,
    ...
  ) {

    kmeans_clusters <-
      stats::kmeans(
        x = select(tof_tibble, {{cluster_vars}}),
        centers = num_clusters,
        ...
      ) %>%
      purrr::pluck("cluster")

    return(tibble(kmeans_cluster = kmeans_clusters))
  }



# tof_cluster_ddpr --------------------

#' Perform k-means clustering on CyTOF data.
#'
#' This function performs k-means clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[stats]{kmeans}}.
#'
#' @param healthy_tibble A `tibble` or `tof_tibble` containing cells from only
#' healthy control samples (i.e. not disease samples).
#'
#' @param cancer_tibble
#'
#' @param healthy_cell_labels A character or integer vector of length `nrow(healthy_tibble)`.
#' Each entry in this vector should represent the cell subpopulation label (or cluster id) for
#' the corresponding row in `healthy_tibble`.
#'
#' @param cluster_vars Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the DDPR clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_cores
#'
#' @param distance_function
#'
#' @param return_distances
#'
#' @param verbose
#'
#' @param ...
#'
#' @return
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_cluster_ddpr <-
  function(
    healthy_tibble,
    cancer_tibble,
    healthy_cell_labels,
    cluster_vars = where(tof_is_numeric),
    distance_function = c("mahalanobis", "cosine", "pearson"),
    num_cores = 1,
    parallel_vars = NULL,
    return_distances = FALSE,
    verbose = FALSE,
    ...
  ) {
    distance_function <-
      match.arg(distance_function, c("mahalanobis", "cosine", "pearson"))

    # build classifier
    classifier_fit <-
      tof_build_classifier(
        healthy_tibble = healthy_tibble,
        healthy_cell_labels = healthy_cell_labels,
        classifier_markers = {{cluster_vars}},
        verbose = verbose
      )

    # apply classifier
    result <-
      tof_apply_classifier(
        cancer_tibble = cancer_tibble,
        classifier_fit = classifier_fit,
        distance_function = distance_function,
        num_cores = num_cores,
        parallel_vars = parallel_vars
      )

    if (!return_distances) {
      result <-
        result %>%
        dplyr::select(tidyselect::all_of(paste0(distance_function, "_cluster")))
    }

    return(result)

  }




# tof_cluster_spade -----------------
# TO DO


# tof_cluster -------------------------

tof_cluster <- function(tof_tibble, method, ...) {
  # TO DO: Fill this in
  return(NULL)
}



