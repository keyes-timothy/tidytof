# clustering.R
# This file contains functions relevant to performing single-cell clustering
# on tof_tibble objects containing CyTOF data.

# tof_cluster_flowsom ----------------------------------------------------------

#' Perform FlowSOM clustering on CyTOF data
#'
#' This function performs FlowSOM clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[FlowSOM]{BuildSOM}} and \code{\link[FlowSOM]{MetaClustering}}.
#'
#' For additional details about the FlowSOM algorithm,
#' see \href{https://pubmed.ncbi.nlm.nih.gov/25573116/}{this paper}.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
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
#' @param seed An integer used to set the random seed for the FlowSOM clustering.
#' Setting this argument explicitly can be useful for reproducibility purposes.
#' Defaults to a random integer
#'
#' @param ... Optional additional parameters that can be passed to the \code{\link[FlowSOM]{BuildSOM}}
#' function.
#'
#' @return A tibble with one column named `.flowsom_cluster` or `.flowsom_metacluster`
#' depending on the value of `perform_metaclustering`. The column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the flowSOM cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
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
    seed,
    ...
  ) {

    # check that flowSOM is installed
    has_flowsom <- requireNamespace(package = "FlowSOM")
    if (!has_flowsom) {
      stop(
        "This function requires the {FlowSOM} package. Install it with this code:\n
          if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\n
             install.packages(\"BiocManager\")\n
          }\n
          BiocManager::install(\"FlowSOM\")\n"
      )
    }

    # set random seed
    if(!missing(seed)) {
      set.seed(seed)
    }

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
        prettyColnames = setNames(clustering_markers, clustering_markers)
      )

    class(fsom) <- "FlowSOM"

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
        silent = TRUE,
        xdim = som_xdim,
        ydim = som_ydim,
        distf = distf,
        ...
      )

    # if no metaclustering, return flowSOM cluster labels
    if (!perform_metaclustering) {
      flowsom_clusters <- som$map$mapping[,1]
      return(dplyr::tibble(.flowsom_cluster = flowsom_clusters))

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
        as.integer() %>%
        as.character()

      return(dplyr::tibble(.flowsom_metacluster = flowsom_metaclusters))
    }
  }


# tof_cluster_phenograph -------------------------------------------------------

#' Perform PhenoGraph clustering on CyTOF data.
#'
#' This function performs PhenoGraph clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[Rphenograph]{Rphenograph}}.
#'
#' For additional details about the Phenograph algorithm,
#' see \href{https://pubmed.ncbi.nlm.nih.gov/26095251/}{this paper}.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the PhenoGraph clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_neighbors An integer indicating the number of neighbors to use when
#' constructing PhenoGraph's k-nearest-neighbor graph. Smaller values emphasize
#' local graph structure; larger values emphasize global graph structure (and
#' will add time to the computation). Defaults to 30.
#'
#' @param distance_function A string indicating which distance function to use
#' for the nearest-neighbor calculation. Options include "euclidean"
#' (the default) and "cosine" distances.
#'
#' @param seed An integer used to set the random seed for the FlowSOM clustering.
#' Setting this argument explicitly can be useful for reproducibility purposes.
#' Defaults to a random integer
#'
#' @param legacy A boolean value indicating if the Chen Lab implementation of
#' phenograph (Rphenograph) should be used (TRUE) or if tidytof's native
#' implementation of phenograph should be used (FALSE; the default).
#'
#' @param ... Optional additional parameters that can be passed to \code{\link[Rphenograph]{Rphenograph}}.
#'
#' @return A tibble with one column named `.phenograph_cluster`. This column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the PhenoGraph cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#'
tof_cluster_phenograph <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_neighbors = 30,
    distance_function = c("euclidean", "cosine"),
    seed,
    legacy = FALSE,
    ...
  ) {
    if (legacy) {
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


      # set random seed
      if(!missing(seed)) {
        set.seed(seed)
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

      phenograph_clusters <- phenograph_result[[2]]$membership

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
      return(dplyr::tibble(.phenograph_cluster = as.character(phenograph_clusters)))

    } else {
      result <-
        phenograph_cluster(
          tof_tibble,
          cluster_cols = {{cluster_cols}},
          num_neighbors = num_neighbors,
          distance_function = distance_function,
          seed = seed,
          ...
        )

      return(result)
    }
  }


# tof_cluster_kmeans -----------------------------------------------------------

#' Perform k-means clustering on CyTOF data.
#'
#' This function performs k-means clustering on CyTOF data using a user-specified
#' selection of input variables/CyTOF measurements. It is mostly a convenient
#' wrapper around \code{\link[stats]{kmeans}}.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the k-means clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_clusters An integer indicating the maximum number of clusters
#' that should be returned Defaults to 20.
#'
#' @param seed An integer used to set the random seed for the FlowSOM clustering.
#' Setting this argument explicitly can be useful for reproducibility purposes.
#' Defaults to a random integer
#'
#' @param ... Optional additional arguments that can be passed to
#' \code{\link[stats]{kmeans}}.
#'
#' @return A tibble with one column named `.kmeans_cluster`. This column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the k-means cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#'
#' @export
#'
#' @importFrom stats kmeans
#' @importFrom purrr pluck
#'
#'
tof_cluster_kmeans <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_clusters = 20,
    seed,
    ...
  ) {

    # set random seed
    if(!missing(seed)) {
      set.seed(seed)
    }

    kmeans_clusters <-
      stats::kmeans(
        x = select(tof_tibble, {{cluster_cols}}),
        centers = num_clusters,
        ...
      ) %>%
      purrr::pluck("cluster")

    return(dplyr::tibble(.kmeans_cluster = as.character(kmeans_clusters)))
  }



# tof_cluster_ddpr -------------------------------------------------------------

#' Perform developmental clustering on CyTOF data.
#'
#' This function performs distance-based clustering on CyTOF data
#' by sorting cancer cells (passed into the function as `tof_tibble`) with
#' their most phenotypically similar healthy cell subpopulation (passed into the
#' function using `healthy_tibble`). For details about
#' the algorithm used to perform the clustering, see \href{https://pubmed.ncbi.nlm.nih.gov/29505032/}{this paper}.
#'
#'
#' @param tof_tibble A `tibble` or `tof_tbl` containing cells to be classified
#' into their nearest healthy subpopulation (generally cancer cells).
#'
#' @param healthy_tibble A `tibble` or `tof_tibble` containing cells from only
#' healthy control samples (i.e. not disease samples).
#'
#' @param healthy_label_col An unquoted column name indicating which column in
#' `healthy_tibble` contains the subpopulation label (or cluster id) for
#' each cell in `healthy_tibble`.
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
#' @param parallel_cols Optional. Unquoted column names indicating which columns in `tof_tibble` to
#' use for breaking up the data in order to parallelize the classification using
#' `foreach` on a `doParallel` backend.
#' Supports tidyselect helpers.
#'
#' @param return_distances A boolean value indicating whether or not the returned
#' result should include only one column, the cluster ids corresponding to each row
#' of `tof_tibble` (return_distances = FALSE, the default), or if the returned
#' result should include additional columns representing the distance between each
#' row of `tof_tibble` and each of the healthy subpopulation centroids
#' (return_distances = TRUE).
#'
#' @param verbose  A boolean value indicating whether progress updates should be
#' printed during developmental classification. Default is FALSE.
#'
#' @return  If `return_distances = FALSE`, a tibble with one column named
#' `{distance_function}_cluster`, a character vector of length `nrow(tof_tibble)`
#' indicating the id of the developmental cluster to which each cell
#' (i.e. each row) in `tof_tibble` was assigned.
#'
#' If `return_distances = TRUE`, a tibble with `nrow(tof_tibble)` rows and `nrow(classifier_fit) + 1`
#' columns. Each row represents a cell from `tof_tibble`, and `nrow(classifier_fit)`
#' of the columns represent the distance between the cell and each of the healthy
#' subpopulations' cluster centroids. The final column represents the cluster id of
#' the healthy subpopulation with the minimum distance to the cell represented
#' by that row.
#'
#' If `return_distances = FALSE`, a tibble with one column named `.{distance_function}_cluster`.
#' This column will contain an integer vector of length `nrow(tof_tibble)` indicating the id of
#' the developmental cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @importFrom tidyselect all_of
#'
#'
tof_cluster_ddpr <-
  function(
    tof_tibble,
    healthy_tibble,
    healthy_label_col,
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

    # check that healthy_tibble exists
    if (missing(healthy_tibble)) {
      stop("DDPR clustering requires the specification of a healthy_tibble.")
    }

    # build classifier
    classifier_fit <-
      tof_build_classifier(
        healthy_tibble = dplyr::select(healthy_tibble, -{{healthy_label_col}}),
        healthy_cell_labels = pull(healthy_tibble, {{healthy_label_col}}),
        classifier_markers = {{cluster_cols}},
        verbose = verbose
      )

    # apply classifier
    if(missing(parallel_cols)) {
      result <-
        tof_apply_classifier(
          cancer_tibble = tof_tibble,
          classifier_fit = classifier_fit,
          distance_function = distance_function,
          num_cores = num_cores
        )

    } else {
      result <-
        tof_apply_classifier(
          cancer_tibble = tof_tibble,
          classifier_fit = classifier_fit,
          distance_function = distance_function,
          num_cores = num_cores,
          parallel_vars = {{parallel_cols}}
        )

    }

    # return result
    result <-
      result %>%
      dplyr::rename_with(.fn = function(x) paste0(".", x))

    if (!return_distances) {
      result <-
        result %>%
        dplyr::select(tidyselect::all_of(paste0(".", distance_function, "_cluster")))
    }

    return(result)
  }



# tof_cluster ------------------------------------------------------------------

#' Cluster CyTOF data.
#'
#' This function is a wrapper around {tidytof}'s tof_cluster_* function family.
#' It performs clustering on CyTOF data using a user-specified method (of 5 choices)
#' and each method's corresponding input parameters
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param method A string indicating which clustering methods should be used. Valid
#' values include "flowsom", "phenograph", "kmeans", "ddpr", and "xshift".
#'
#' @param ... Additional arguments to pass onto the `tof_cluster_*`
#' function family member corresponding to the chosen method.
#'
#' @param add_col A boolean value indicating if the output should column-bind the
#' cluster ids of each cell as a new column in `tof_tibble` (TRUE, the default) or if
#' a single-column tibble including only the cluster ids should be returned (FALSE).
#'
#' @return A `tof_tbl` or `tibble` If add_col = FALSE, it will have a single column encoding
#' the cluster ids for each cell in `tof_tibble`. If add_col = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the cluster ids.
#'
#' @export
#'
tof_cluster <- function(tof_tibble, method, ..., add_col = TRUE) {

  if (method == "flowsom") {
    clusters <-
      tof_tibble %>%
      tof_cluster_flowsom(...)

  } else if (method == "phenograph") {
    clusters <-
      tof_tibble %>%
      tof_cluster_phenograph(...)

  } else if (method == "kmeans") {
    clusters <-
      tof_tibble %>%
      tof_cluster_kmeans(...)

  } else if (method == "ddpr") {
    clusters <-
      tof_tibble %>%
      tof_cluster_ddpr(...)

  } else if (method == "xshift") {
    # clusters <-
    #   tof_tibble %>%
    #   tof_cluster_xshift(...)
    stop("X-shift is not yet implemented.")

  } else {
    stop("Not a valid clustering method.")
  }

  if (add_col == TRUE) {
    result <-
      dplyr::bind_cols(tof_tibble, clusters)
  } else {
    result <- clusters
  }

  return(result)
}



