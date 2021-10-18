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
#' @param seed An integer used to set the random seed for the FlowSOM clustering.
#' Setting this argument explicitly can be useful for reproducibility purposes.
#' Defaults to a random integer
#'
#' @param ... Optional additional parameters that can be passed to \code{\link[Rphenograph]{Rphenograph}}.
#'
#' @return A tibble with one column named `.phenograph_cluster`. This column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the PhenoGraph cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @export
#'
#' @importFrom igraph membership
#'
tof_cluster_phenograph <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_neighbors = 30,
    seed,
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
    return(dplyr::tibble(.phenograph_cluster = as.character(phenograph_clusters)))
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
        ) %>%
        dplyr::rename_with(.fn = ~ paste0(".", .x))

    } else {
      result <-
        tof_apply_classifier(
          cancer_tibble = tof_tibble,
          classifier_fit = classifier_fit,
          distance_function = distance_function,
          num_cores = num_cores,
          parallel_vars = {{parallel_cols}}
        ) %>%
        dplyr::rename_with(.fn = ~ paste0(".", .x))

    }

    # return desired result
    if (!return_distances) {
      result <-
        result %>%
        dplyr::select(tidyselect::all_of(paste0(".", distance_function, "_cluster")))
    }

    return(result)
  }



# tof_cluster_xshift -----------------------------------------------------------
tof_cluster_xshift <-
  function(
    tof_tibble,
    # cluster_cols = some_default, - have to add this!!!
    k = max(20, nrow(tof_tibble)),
    distance_function = c("cosine", "euclidean"),
    p_value = 0.01
  ) {
    # check distance_function argument
    distance_function <-
      match.arg(distance_function, choices = c("cosine", "euclidean"))

    # calculate "z" value for nearest-neighbor search during step 2 according
    # to statistical criteria in the x-shift methods section
    z <- floor(-log(p_value / nrow(tof_tibble), base = 2))


    # step 1a - find nearest neighbors of each cell in tof_tibble
    # note that the number of neighbors we find is the larger of k and z,
    # as this allows us not to repeat the calculation needlessly in step 2
    nn_result <-
      tof_tibble %>%
      tof_find_knn(k = max(z, k), distance_function = distance_function)

    # step 1b - compute local densities of each cell using knn density estimation
    # using only k neighbors (if z > k, only use the first k neighbors)
    densities <-
      xshift_compute_local_densities(
        neighbor_ids = nn_result$neighbor_ids[, 1:k],
        neighbor_distances = nn_result$neighbor_distances[, 1:k]
      )

    # step 2 - identify candidate centroids from all cells in the dataset
    candidates <-
      xshift_find_candidate_centroids(
        neighbor_ids = nn_result$neighbor_ids,
        densities = densities
      )

    ## step 3 - merge clusters

    # step 3a - find which candidate centroids are gabriel neighbors
    gabriel_results <-
      xshift_find_gabriel_neighbors(
        candidates = candidates,
        distance_function = distance_function
      )

    # candidate centroids' protein measurements extracted from tof_tibble
    centroids <- gabriel_results$centroids

    # tibble containing information about all gabriel pairs from all 2-wise
    # combinations of candidate nodes
    gabriel_pairs <- gabriel_results$gabriel_pairs


    # step 3b - for all gabriel pairs,
    #           test if there is a minimum of density on the segment connecting
    #           the centroids. If there is not, the lower-density centroid should
    #           be merged with the larger-density centroid.
    gabriel_pairs <-
      xshift_test_density_minima(
        gabriel_pairs = gabriel_pairs,
        centroids = centroids,
        k = k,
        tof_tibble = tof_tibble
      )

    # step 4 - Merge centroids based on density connections
    merge_information <-
      gabriel_pairs %>%
      dplyr::filter(should_merge) %>%
      dplyr::select(-should_merge) %>%
      tidyr::pivot_longer(
        cols = c(centroid_1, centroid_2),
        names_to = "centroid",
        values_to = "merge_from"
      ) %>%
      dplyr::select(-centroid) %>%
      # potentially remove this line
      dplyr::filter(merge_into != merge_from)

    #########################
    merge_information_from <-
      merge_information %>%
      dplyr::group_by(merge_from) %>%
      dplyr::summarize(connected_to = list(unique(merge_into)))

    merge_information_into <-
      merge_information %>%
      dplyr::group_by(merge_into) %>%
      dplyr::summarize(connected_to = list(unique(merge_from)))
    #########################

    candidates <-
      candidates %>%
      dplyr::left_join(
        merge_information_from,
        by = c("cell_id" = "merge_from")
      ) %>%
      dplyr::mutate(
        connected_to =
          purrr::map2(
            .x = closest_neighbor,
            .y = connected_to,
            .f = ~ as.numeric(na.omit(unique(c(.x, .y))))
          )
      ) %>%
      dplyr::select(-higher_density_neighbors)

    eliminated_candidates <-
      merge_information_from %>%
      dplyr::pull(merge_from)

    candidates <-
      candidates %>%
      dplyr::mutate(
        still_candidate =
          dplyr::if_else(cell_id %in% eliminated_candidates, FALSE, is_candidate)
      )

    cluster_ids <-
      candidates %>%
      dplyr::filter(still_candidate) %>%
      dplyr::transmute(cell_id, cluster_id = 1:dplyr::n())

    # density clusters for candidate centroids
    density_clusters <-
      candidates %>%
      dplyr::left_join(cluster_ids, by = "cell_id")

    edge_list <-
      candidates %>%
      select(cell_id, connected_to) %>%
      mutate(
        connected_to = purrr::map(connected_to, ~ tibble(to = .x))
      ) %>%
      tidyr::unnest(cols = connected_to) %>%
      rename(from = cell_id)

    node_list <-
      candidates %>%
      select(cell_id)

    my_graph <-
      tidygraph::tbl_graph(
        nodes = node_list,
        edges = edge_list,
        directed = FALSE
      )

    xshift_clusters <-
      my_graph %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(
        group = tidygraph::group_components(type = "strong")
      ) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(xshift_cluster = as.character(group)) %>%
      dplyr::arrange(cell_id) %>%
      dplyr::select(xshift_cluster)



    # step 5 - Merge centroids based on mahalanobis distance (until all clusters
    # have mahalanobis distance of 2.0 or more between them)

    # TO DO

    # return result
    return(xshift_clusters)



    # result <-
    #   list(
    #     z = z,
    #     nn_result = nn_result,
    #     densities = densities,
    #     candidates = candidates,
    #     merge_information = merge_information,
    #     merge_information_from = merge_information_from,
    #     merge_information_into = merge_information_into,
    #     gabriel_pairs = gabriel_pairs,
    #     cluster_ids = cluster_ids,
    #     density_clusters = density_clusters
    #   )





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
#' values include "flowsom", "phenograph", "kmeans", "ddpr" (although this will
#' throw an error and redirect the user), and "xshift".
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
    clusters <-
      tof_tibble %>%
      tof_cluster_xshift(...)

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



