

# single-cell visualizations ----------------------------

tof_plot_sc_histograms <-
  function(
    tof_tibble,
    x_col,
    y_col,
    facet_cols,
    theme = ggplot2::theme_bw(),
    ...
  ) {
    stop("This function is not yet implemented!")
  }

#' Title
#'
#' description
#'
#' @param tof_tibble
#'
#' @param dr_cols
#'
#' @param color_col
#'
#' @param facet_cols
#'
#' @param dr_method
#'
#' @param point_alpha
#'
#' @param theme
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
tof_plot_sc_dr <-
  function(
    tof_tibble,
    dr_cols,
    color_col,
    facet_cols,
    dr_method = c("pca", "tsne", "umap"),
    point_alpha = 1,
    theme = ggplot2::theme_bw(),
    ...# optional additional arguments to the specified tof_dr_* family function
  ) {

    # if no dr_cols are specified, use the dr_method to compute them
    if (missing(dr_cols)) {
      # if there's no dr_method specified, just use PCA (for speed)
      if (dr_method == c("pca", "tsne", "umap")) {
        message("No dr_cols were specified, and no dr_method was specified.
                Performing PCA as the default dimensionality reduction method.")
      }
      # check dr_method columns
      dr_method <- rlang::arg_match(dr_method)
      dr_tibble <-
        tof_tibble %>%
        tof_dr(method = dr_method, ...)

    # if there are dr_cols specified, just use those
    } else {
      # check dr_cols - there should only be two
      dr_tibble <-
        tof_tibble %>%
        dplyr::select({{dr_cols}})

      num_dr_cols <-
        dr_tibble %>%
        ncol()

      if (num_dr_cols != 2) {
        stop("2 dimensionality reduction columns must be selected.")
      }
    }

    if (missing(color_col)) {
      shape = 16
    } else {
      shape = 21
    }

    # create plot
    result <-
      tof_tibble %>%
      ggplot2::ggplot(aes(x = dr_tibble[[1]], y = dr_tibble[[2]], fill = {{color_col}})) +
      ggplot2::geom_point(shape = shape, alpha = point_alpha) +
      labs(
        x = colnames(dr_tibble)[[1]],
        y = colnames(dr_tibble)[[2]]
      )

    if (!missing(facet_cols)) {
      result <-
        result +
        ggplot2::facet_wrap(facets = vars({{facet_cols}}))
    }

    return(result)

  }

#' Title
#'
#' description
#'
#' @param tof_tibble
#'
#' @param knn_cols
#'
#' @param color_col
#'
#' @param facet_cols
#'
#' @param num_neighbors
#'
#' @param graph_type
#'
#' @param graph_layout
#'
#' @param distance_function
#'
#' @param knn_error
#'
#' @param edge_alpha
#'
#' @param node_size
#'
#' @param theme
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
#'
tof_plot_sc_layout <-
  function(
    tof_tibble,
    knn_cols = where(tof_is_numeric),
    color_col,
    facet_cols,
    num_neighbors = 5,
    graph_type = c("weighted", "unweighted"),
    graph_layout = "fr",
    distance_function = c("euclidean", "cosine"),
    knn_error = 0,
    edge_alpha = 0.25,
    node_size = 2,
    theme = ggplot2::theme_void(),
    ...
  ) {


    # check distance function
    distance_function <- rlang::arg_match(distance_function)

    # check graph type
    graph_type = rlang::arg_match(graph_type)

    knn_data <-
      tof_tibble %>%
      # select only knn_cols
      dplyr::select({{knn_cols}}) %>%
      tof_find_knn(
        k = num_neighbors,
        distance_function = distance_function,
        eps = knn_error
      )

    # extract knn_ids and put them into long format
    knn_ids <-
      knn_data %>%
      purrr::pluck("neighbor_ids")
    colnames(knn_ids) <- 1:ncol(knn_ids)

    knn_ids <-
      knn_ids %>%
      tibble::as_tibble() %>%
      dplyr::mutate(from = seq(from = 1, to = nrow(knn_ids), by = 1)) %>%
      tidyr::pivot_longer(
        cols = -from,
        names_to = "neighbor_index",
        values_to = "to"
      )

    if (graph_type == "weighted") {
      # extract knn distances and put them into long format
      knn_dists <-
        knn_data %>%
        purrr::pluck("neighbor_distances")
      colnames(knn_dists) <- 1:ncol(knn_dists)

      knn_dists <-
        knn_dists %>%
        tibble::as_tibble() %>%
        dplyr::mutate(from = seq(from = 1, to = nrow(knn_dists), by = 1)) %>%
        tidyr::pivot_longer(
          cols = -from,
          names_to = "neighbor_index",
          values_to = "distance"
        )

      # join knn distances with knn ids for final edge tibble
      edge_tibble <-
        knn_ids %>%
        dplyr::left_join(knn_dists, by = (c("from", "neighbor_index")))

      if (distance_function == "euclidean") {
        edge_tibble <-
          edge_tibble %>%
          dplyr::mutate(weight = 1 / (1 + distance))
      } else {
        edge_tibble <-
          edge_tibble %>%
          dplyr::mutate(weight = 1 - distance)
      }

    } else {
      edge_tibble <-
        knn_ids
    }

    # make the knn_graph
    knn_graph <-
      tidygraph::tbl_graph(
        nodes = tof_tibble,
        edges = edge_tibble,
        directed = FALSE
      )

    # make the initial ggraph call with or without weights
    if (graph_type == "weighted") {
      knn_plot <-
        ggraph::ggraph(
          graph = knn_graph,
          layout = graph_layout,
          weights = weight,
          ...
        )
    } else {
      knn_plot <-
        ggraph::ggraph(
          graph = knn_graph,
          layout = graph_layout,
          ...
        )
    }

    knn_plot <-
      knn_plot +
      ggraph::geom_edge_link(alpha = edge_alpha) +
      ggraph::geom_node_point(
        aes(fill = {{color_col}}),
        shape = 21,
        size = node_size
      )

    if (!missing(facet_cols)) {
      knn_plot <-
        knn_plot +
        ggraph::facet_nodes(facets = vars({{facet_cols}}))
    }

    return(knn_plot + theme)
  }



# community-level visualizations ----------------------------

tof_plot_community_layout <-
  function(
    ...
  ) {
    stop("This function is not yet implemented!")
  }

tof_plot_community_volcano <-
  function(
    ...
  ) {
    stop("This function is not yet implemented!")

  }


# patient-level visualizations --------------------------

tof_plot_patient_features <-
  function(
    ...
  ) {
    stop("This function is not yet implemented!")

  }

tof_plot_patient_model <-
  function(
    model_fit,
    new_data,
    ...
  ) {
    stop("This function is not yet implemented!")

    # find model type from model_fit
    NULL

    # make plot depending on the input model_fit
    if (model_type == "regression") {
      # make scatterplot of real y values vs. predictions
      NULL
    } else if (model_type == "classification") {
      # make an ROC curve
      NULL
    } else {
      # make some kind of survival curve
      NULL
    }
  }


