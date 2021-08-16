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
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
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
#' map calculations. Options are "euclidean" (the default), "manhattan", "chebyshev",
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
#' @return A tibble with one column named `flowsom_cluster` or `flowsom_metacluster`
#' depending on the value of `perform_metaclustering`. The column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the flowSOM cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @importFrom FlowSOM BuildSOM
#' @importFrom FlowSOM BuildMST
#'
#'
#' @export
#'
#'
tof_cluster_flowsom <-
  function(
    tof_tibble = NULL,
    cluster_cols = where(tof_is_numeric),
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
      dplyr::select({{cluster_cols}}) %>%
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
      # convert character distance function name to a number that BuildSOM understands
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
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
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
#' \dontrun{
#' tof_tibble <-
#'    tof_read_data(tidytof_example_data("aml")) %>%
#'    dplyr::slice_sample(n = 10000)
#'
#' phenograph_clusters <-
#'    tof_cluster_phenograph(
#'      tof_tibble,
#'      cluster_cols = contains("CD", ignore.case = FALSE)
#'    )
#'}
#'
tof_cluster_phenograph <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_neighbors = 30,
    ...
  ) {
    # check that Rphenograph is installed
    has_phenograph <- requireNamespace(package = "Rphenograph")
    if (!has_phenograph) {
      stop(
        "This function requires the {Rphenograph} package. Install it with this code:\n
           if(!require(devtools)){
              install.packages(\"devtools\")
           }
           devtools::install_github(\"JinmiaoChenLab/Rphenograph\")\n"
      )
    }

    invisible(
      utils::capture.output(
        suppressMessages(
          phenograph_result <-
            Rphenograph::Rphenograph(
              data = dplyr::select(tof_tibble, {{cluster_cols}}),
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
#' @param tof_tibble A `tof_tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the kmeans clusters. Defaults to all numeric columns
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
#'
tof_cluster_kmeans <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_clusters = 20,
    ...
  ) {

    kmeans_clusters <-
      stats::kmeans(
        x = select(tof_tibble, {{cluster_cols}}),
        centers = num_clusters,
        ...
      ) %>%
      purrr::pluck("cluster")

    return(tibble(kmeans_cluster = kmeans_clusters))
  }



# tof_cluster_ddpr --------------------

#' Perform developmental clustering on CyTOF data.
#'
#' This function performs distance-based clustering on CyTOF data
#' by sorting cancer cells (passed into the function as `cancer_tibble`) with
#' their most phenotypically similar healthy cell subpopulation (passed into the
#' function using `healthy_tibble` and `healthy_cell_labels`). For details about
#' the algorithm used to perform the clustering, see \href{https://pubmed.ncbi.nlm.nih.gov/29505032/}{this paper}.
#'
#' @param healthy_tibble A `tibble` or `tof_tibble` containing cells from only
#' healthy control samples (i.e. not disease samples).
#'
#' @param cancer_tibble A `tibble` or `tof_tibble` containing cells to be classified
#' into their nearest healthy subpopulation (generally cancer cells).
#'
#' @param healthy_cell_labels A character or integer vector of length `nrow(healthy_tibble)`.
#' Each entry in this vector should represent the cell subpopulation label (or cluster id) for
#' the corresponding row in `healthy_tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the DDPR clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param distance_function A string indicating which distance function should
#' be used to perform the classification. Options are "mahalanobis" (the default),
#' "cosine", and "pearson".
#'
#' @param num_cores An integer indicating the number of CPU cores used to parallelize
#' the classification. Defaults to 1 (a single core).
#'
#' @param parallel_cols Optional. Unquoted column names indicating which columns in `cancer_tibble` to
#' use for breaking up the data in order to parallelize the classification using
#' `foreach` on a `doParallel` backend.
#' Supports tidyselect helpers.
#'
#' @param return_distances A boolean value indicating whether or not the returned
#' result should include only one column, the cluster ids corresponding to each row
#' of `cancer_tibble` (return_distances = FALSE, the default), or if the returned
#' result should include additional columns representing the distance between each
#' row of `cancer_tibble` and each of the healthy subpopulation centroids
#' (return_distances = TRUE).
#'
#' @param verbose  A boolean value indicating whether progress updates should be
#' printed during developmental classification. Default is FALSE.
#'
#' @return  If `return_distances = FALSE`, a tibble with one column named
#' `{distance_function}_cluster`, a character vector of length `nrow(cancer_tibble)`
#' indicating the id of the developmental cluster to which each cell
#' (i.e. each row) in `cancer_tibble` was assigned.
#'
#' If `return_distances = TRUE`, a tibble with `nrow(cancer_tibble)` rows and `nrow(classifier_fit) + 1`
#' columns. Each row represents a cell from `cancer_tibble`, and `nrow(classifier_fit)`
#' of the columns represent the distance between the cell and each of the healthy
#' subpopulations' cluster centroids. The final column represents the cluster id of
#' the healthy subpopulation with the minimum distance to the cell represented
#' by that row.
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
    cluster_cols = where(tof_is_numeric),
    distance_function = c("mahalanobis", "cosine", "pearson"),
    num_cores = 1L,
    parallel_cols,
    return_distances = FALSE,
    verbose = FALSE
  ) {

    # check distance function
    distance_function <-
      match.arg(distance_function, c("mahalanobis", "cosine", "pearson"))

    # build classifier
    classifier_fit <-
      tof_build_classifier(
        healthy_tibble = healthy_tibble,
        healthy_cell_labels = healthy_cell_labels,
        classifier_markers = {{cluster_cols}},
        verbose = verbose
      )

    # apply classifier
    if(missing(parallel_cols)) {
      result <-
        tof_apply_classifier(
          cancer_tibble = cancer_tibble,
          classifier_fit = classifier_fit,
          distance_function = distance_function,
          num_cores = num_cores
        )
    } else {
      result <-
        tof_apply_classifier(
          cancer_tibble = cancer_tibble,
          classifier_fit = classifier_fit,
          distance_function = distance_function,
          num_cores = num_cores,
          parallel_vars = {{parallel_cols}}
        )
    }

    # return desired result
    if (!return_distances) {
      result <-
        result %>%
        dplyr::select(tidyselect::all_of(paste0(distance_function, "_cluster")))
    }

    return(result)
  }



# tof_cluster_xshift -----------------
tof_cluster_xshift <-
  function(
    tof_tibble,
    k = max(20, nrow(tof_tibble)),
    distance_function = c("cosine", "euclidean"),
    p_value = 0.01
  ) {
    # check distance_function argument
    distance_function <-
      match.arg(distance_function, choices = c("cosine", "euclidean"))

    # calculate "z" value for nearest-neighbor search during step 2
    z <-
      -log(p_value / nrow(tof_tibble), base = 2) %>%
      floor()

    # step 1a - find nearest neighbors of each cell in tof_tibble
    nn_result <-
      tof_tibble %>%
      tof_find_knn(k = max(z, k), distance_function = distance_function)

    # step 1b - compute local densities of each cell using knn density estimation
    densities <-
      xshift_compute_local_densities(
        neighbor_ids = nn_result$neighbor_ids[, 1:k],
        neighbor_distances = nn_result$neighbor_distances[, 1:k]
      )

    # step 2 -



  }


# tof_cluster -------------------------

tof_cluster <- function(tof_tibble, method, ...) {
  # TO DO: Fill this in
  return(NULL)
}



