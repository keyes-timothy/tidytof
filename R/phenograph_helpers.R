# phenograph_helpers.R
# This file contains functions relevant to performing phenograph clustering
# using a tidygraph implementation



#'
#' @importFrom rlang arg_match
#'
#' @importFrom dplyr select
#' @importFrom dplyr as_tibble
#' @importFrom dplyr filter
#' @importFrom dplyr tibble
#' @importFrom dplyr transmute
#'
#' @importFrom tidygraph tbl_graph
#' @importFrom tidygraph mutate
#' @importFrom tidygraph group_louvain
#' @importFrom tidygraph graph_modularity
#' @importFrom tidygraph as_tibble
#'
phenograph_cluster <-
  function(
    tof_tibble,
    cluster_cols = where(tof_is_numeric),
    num_neighbors = 30,
    distance_function = c("euclidean", "cosine"),
    seed,
    ...# optional additional arguments to tof_find_knn
  ) {
    # check distance_function
    distance_function <-
      rlang::arg_match(distance_function, values = c("euclidean", "cosine"))

    ## find knn for each cell in tof_tibble
    nn_result <-
      tof_tibble %>%
      dplyr::select({{cluster_cols}}) %>%
      tof_find_knn(k = num_neighbors, distance_function = distance_function, ...)

    ## extract knn_ids (the ids of neighbors for each row/cell in tof_tibble)
    knn_ids <- nn_result$neighbor_ids
    colnames(knn_ids) <- 1:ncol(knn_ids)

    # construct the second graph - a graph in which edges between cells are weighted
    # based on the number of neighbors they share in the first graph
    # this is called the Jaccard Similarity Coefficient between cells

    # I have to write an RCpp function for this
    jaccards <- find_jaccard_coefficients(knn_ids)
    colnames(jaccards) <- c("from", "to", "jaccard")
    jaccard_edges <-
      jaccards %>%
      dplyr::as_tibble() %>%
      dplyr::filter(.data$jaccard > 0)

    jaccard_graph <-
      tidygraph::tbl_graph(
        nodes = dplyr::tibble(name = 1:nrow(tof_tibble)),
        edges = jaccard_edges,
        directed = FALSE
      )

    # perform louvain clustering on the jaccard graph
    jaccard_graph <-
      jaccard_graph %>%
      tidygraph::mutate(
        .phenograph_cluster = tidygraph::group_louvain(weights = .data$jaccard)
      )

    modularity <-
      jaccard_graph %>%
      tidygraph::mutate(
        modularity =
          tidygraph::graph_modularity(group = .data$.phenograph_cluster)
      ) %>%
      dplyr::as_tibble() %>%
      dplyr::pull(.data$modularity) %>%
      unique()

    clusters <-
      jaccard_graph %>%
      tidygraph::as_tibble() %>%
      dplyr::transmute(
        .phenograph_cluster = as.character(.data$.phenograph_cluster)
      )

    attr(x = clusters, which = "modularity") <- modularity

    # return result
    return(clusters)
  }
