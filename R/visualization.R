

# single-cell visualizations ----------------------------

tof_plot_cells_histogram <-
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
#' @param tof_tibble TO DO
#'
#' @param color_col TO DO
#'
#' @param facet_cols TO DO
#'
#' @param point_alpha TO DO
#'
#' @param theme TO DO
#'
#' @param ... TO DO
#'
#' @param embedding_cols TO DO
#'
#' @param embedding_method TO DO
#'
#' @return TO DO
#'
#' @family visualization functions
#'
#' @export
#'
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 vars
#'
#'
tof_plot_cells_embedding <-
  function(
    tof_tibble,
    embedding_cols,
    color_col,
    facet_cols,
    embedding_method = c("pca", "tsne", "umap"),
    point_alpha = 1,
    theme = ggplot2::theme_bw(),
    ...# optional additional arguments to the specified tof_dr_* family function
  ) {

    # if no dr_cols are specified, use the embedding_method to compute them
    if (missing(embedding_cols)) {
      # if there's no embedding_method specified, just use PCA (for speed)
      if (identical(embedding_method, c("pca", "tsne", "umap"))) {
        message("No embedding_cols were specified, and no embedding_method was specified.
                Performing PCA as the default dimensionality reduction method.")
      }
      # check embedding_method columns
      embedding_method <- rlang::arg_match(embedding_method)
      embedding_method <-
        tof_tibble %>%
        tof_reduce_dimensions(method = embedding_method, ...)

    # if there are embedding_cols specified, just use those
    } else {
      # check embedding_cols - there should only be two
      embed_tibble <-
        tof_tibble %>%
        dplyr::select({{embedding_cols}})

      num_embed_cols <-
        embed_tibble %>%
        ncol()

      if (num_embed_cols != 2) {
        stop("2 dimensionality reduction columns must be selected.")
      }
    }

    if (missing(color_col)) {
      shape = 16
    } else {
      shape = 21
    }

    # create plot tibble for memory efficiency
    if (!missing(facet_cols)) {
      plot_tibble <-
        tof_tibble %>%
        dplyr::select({{color_col}}, {{facet_cols}})
    } else {
      plot_tibble <-
        tof_tibble %>%
        dplyr::select({{color_col}})
    }

    # create plot
    result <-
      plot_tibble %>%
      ggplot2::ggplot(
        ggplot2::aes(x = embed_tibble[[1]], y = embed_tibble[[2]], fill = {{color_col}})
      ) +
      ggplot2::geom_point(shape = shape, alpha = point_alpha) +
      ggplot2::labs(
        x = colnames(embed_tibble)[[1]],
        y = colnames(embed_tibble)[[2]]
      )

    if (!missing(facet_cols)) {
      result <-
        result +
        ggplot2::facet_wrap(facets = ggplot2::vars({{facet_cols}}))
    }
    return(result + theme)
  }

#' Title
#'
#' description
#'
#' @param tof_tibble TO DO
#'
#' @param knn_cols TO DO
#'
#' @param color_col TO DO
#'
#' @param facet_cols TO DO
#'
#' @param num_neighbors TO DO
#'
#' @param graph_type TO DO
#'
#' @param graph_layout TO DO
#'
#' @param distance_function TO DO
#'
#' @param knn_error TO DO
#'
#' @param edge_alpha TO DO
#'
#' @param node_size TO DO
#'
#' @param theme TO DO
#'
#' @param ... TO DO
#'
#' @return TO DO
#'
#' @family visualization functions
#'
#' @export
#'
#' @importFrom dplyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 vars
#'
#' @importFrom ggraph facet_nodes
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#'
#' @importFrom purrr pluck
#'
#' @importFrom rlang arg_match
#'
#' @importFrom tidygraph tbl_graph
#' @importFrom tidygraph select
#'
#' @importFrom tidyr pivot_longer
#'
#'
tof_plot_cells_layout <-
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

    # throw error if color_col is missing
    if (missing(color_col)) {
      stop("color_col must be specified.")
    }

    knn_graph <-
      tof_tibble %>%
      tof_make_knn_graph(
        knn_cols = {{knn_cols}},
        num_neighbors = num_neighbors,
        distance_function = distance_function,
        knn_error = knn_error,
        graph_type = graph_type
      )

    # retain only the needed columns for memory purposes
    if (missing(facet_cols)) {
    plot_graph <-
      knn_graph %>%
      tidygraph::select({{color_col}})
    } else {
      plot_graph <-
        knn_graph %>%
        tidygraph::select({{color_col}}, {{facet_cols}})
    }

    # make the initial ggraph call with or without weights
    if (graph_type == "weighted") {
      knn_plot <-
        ggraph::ggraph(
          graph = plot_graph,
          layout = graph_layout,
          weights = .data$weight,
          ...
        )
    } else if (graph_type == "unweighted") {
      knn_plot <-
        ggraph::ggraph(
          graph = plot_graph,
          layout = graph_layout,
          ...
        )
    } else {
      stop("Not a valid graph_type")
    }

    knn_plot <-
      knn_plot +
      ggraph::geom_edge_link(alpha = edge_alpha) +
      ggraph::geom_node_point(
        ggplot2::aes(fill = {{color_col}}),
        shape = 21,
        size = node_size
      )

    if (!missing(facet_cols)) {
      knn_plot <-
        knn_plot +
        ggraph::facet_nodes(facets = ggplot2::vars({{facet_cols}}))
    }

    return(knn_plot + theme)
  }



# community-level visualizations ----------------------------

#' Visualize clusters using a minimum spanning tree (MST).
#'
#' @param tof_tibble TO DO
#' @param cluster_col TO DO
#' @param group_cols TO DO
#' @param knn_cols TO DO
#' @param color_col TO DO
#' @param num_neighbors TO DO
#' @param graph_type TO DO
#' @param graph_layout TO DO
#' @param central_tendency_function TO DO
#' @param distance_function TO DO
#' @param knn_error TO DO
#' @param edge_alpha TO DO
#' @param node_size TO DO
#' @param theme TO DO
#' @param ... TO DO
#'
#' @return TO DO
#'
#' @export
#'
#' @importFrom dplyr count
#' @importFrom dplyr left_join
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_void
#'
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#'
#' @importFrom rlang arg_match
#'
#' @importFrom tidygraph convert
#' @importFrom tidygraph to_minimum_spanning_tree
#'
tof_plot_cluster_mst <-
  function(
    tof_tibble,
    cluster_col,
    group_cols,
    knn_cols = where(tof_is_numeric),
    color_col,
    #facet_cols,
    num_neighbors = 5L,
    graph_type = c("weighted", "unweighted"),
    graph_layout = "fr",
    central_tendency_function = stats::median,
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

    # throw error if color_col is missing
    if (missing(color_col)) {
      stop("color_col must be specified.")
    }

    if (missing(group_cols)) {
      group_cols <- NULL
    }

    # if (missing(facet_cols)) {
    #   has_facets <- FALSE
    #   facet_cols <- NULL
    # } else {
    #   has_facets <- TRUE
    # }

    # summarize the clusters
    cluster_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = {{knn_cols}},
        group_cols = c({{group_cols}}, {{color_col}}), #, {{facet_cols}}),
        central_tendency_function = central_tendency_function
      )

    cluster_sizes <-
      tof_tibble %>%
      dplyr::count(
        {{group_cols}}, {{cluster_col}}, {{color_col}},
        name = ".cluster_size"
      )

    cluster_tibble <-
      suppressMessages(
        cluster_tibble %>%
          dplyr::left_join(cluster_sizes)
      )

    # make the knn graph
    knn_graph <-
      cluster_tibble %>%
      tof_make_knn_graph(
        knn_cols = {{knn_cols}},
        num_neighbors = num_neighbors,
        distance_function = distance_function,
        knn_error = knn_error,
        graph_type = graph_type
      )

    # make the mst -------------------------------------------------------------

    # make the edges depending on whether the graph is weighted or unweighted
    if (graph_type == "weighted") {
      mst <-
        knn_graph %>%
        tidygraph::convert(
          .f = tidygraph::to_minimum_spanning_tree,
          weights = .data$weight
        )

      mst_plot <-
        mst %>%
        ggraph::ggraph(layout = graph_layout, weights = .data$weight) +
        ggraph::geom_edge_link(alpha = edge_alpha)

    } else if (graph_type == "unweighted") {
      mst <-
        knn_graph %>%
        tidygraph::convert(.f = tidygraph::to_minimum_spanning_tree)

      mst_plot <-
        mst %>%
        ggraph::ggraph(layout = graph_layout) +
        ggraph::geom_edge_link(alpha = edge_alpha)
    }

    # make the nodes depending on whether a constant size or a size
    # proportional to the cluster sizes was requested
    if (node_size == "cluster_size") {
      mst_plot <-
        mst_plot +
        ggraph::geom_node_point(
          ggplot2::aes(fill = {{color_col}}, size = .data$.cluster_size),
          shape = 21
        )
    } else {
      mst_plot <-
        mst_plot +
        ggraph::geom_node_point(
          ggplot2::aes(fill = {{color_col}}),
          shape = 21,
          size = node_size
        )
    }

    # if (has_facets) {
    #   mst_plot <-
    #     mst_plot +
    #     ggraph::facet_nodes(facets = ggplot2::vars({{facet_cols}}))
    # }

    # return result
    result <- mst_plot + theme
    return(result)
  }

tof_plot_cluster_volcano <-
  function(
    ...
  ) {
    stop("This function is not yet implemented!")

  }


# sample-level visualizations --------------------------------------------------

tof_plot_sample_features <-
  function(
    ...
  ) {
    stop("This function is not yet implemented!")

  }

tof_plot_sample_model <-
  function(
    model_fit,
    new_data,
    ...
  ) {
    stop("This function is not yet implemented!")

    # find model type from model_fit
    model_type <- NULL

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


