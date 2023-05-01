

# single-cell visualizations ---------------------------------------------------

#' Plot marker expression density plots
#'
#' This function plots marker expression density plots for a user-specified
#' column in a tof_tbl. Optionally, cells can be grouped to plot multiple
#' vertically-arranged density plots
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param marker_col An unquoted column name representing which column in `tof_tibble`
#' (i.e. which CyTOF protein measurement) should be included in the feature extraction
#' calculation.
#'
#' @param group_col Unquoted column names representing which column in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups to be plotted
#' as separate histograms. Defaults to plotting without subgroups.
#'
#' @param num_points The number of points along the full range of `marker_col` at
#' which the density should be calculated
#'
#' @param theme The ggplot2 theme for the plot. Defaults to
#' \code{\link[ggplot2]{theme_bw}}
#'
#' @param use_ggridges A boolean value indicting if
#' \code{\link[ggridges]{geom_ridgeline}} should be used to plot overlain
#' histograms. Defaults to FALSE. If TRUE, the ggridges package must be installed.
#'
#' @param scale Use to set the `scale` argument in \code{\link[ggridges]{geom_ridgeline}},
#' which controls how far apart (vertically) density plots are arranged along the
#' y-axis. Defaults to 1.
#'
#' @param ... Additional optional arguments to send to \code{\link[ggridges]{geom_ridgeline}}.
#'
#' @return A ggplot object
#'
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 vars
#'
#' @importFrom purrr map
#'
#' @importFrom rlang check_installed
#'
#' @importFrom stats density
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @importFrom tidyselect everything
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(c("a", "b"), size = 1000, replace = TRUE)
#'     )
#'
#' density_plot <-
#'     tof_plot_cells_density(
#'         tof_tibble = sim_data,
#'         marker_col = cd45,
#'         group_col = cluster_id
#'     )
#'
#'
tof_plot_cells_density <-
  function(
    tof_tibble,
    marker_col,
    group_col,
    num_points = 512,
    theme = ggplot2::theme_bw(),
    use_ggridges = FALSE,
    scale = 1,
    ...
  ) {
    # collect marker column name as a string
    marker_colname <-
      tof_tibble %>%
      dplyr::select({{marker_col}}) %>%
      colnames()

    # calculate the density for each group independently
    # ggplot2 can do this, but will store every cell in the resulting plot
    # so this is more memory (and therefore speed) efficient
    marker_tibble <-
      tof_tibble %>%
      dplyr::select({{marker_col}}, {{group_col}}) %>%
      tidyr::nest(data = {{marker_col}}) %>%
      dplyr::mutate(
        densities =
          purrr::map(
            .x = .data$data,
            .f = ~
              stats::density(dplyr::pull(.x, {{marker_col}}), n = num_points)
          ),
        expression = purrr::map(.x = .data$densities, .f = ~ .x$x),
        density = purrr::map(.x = .data$densities, .f = ~ .x$y)
      ) %>%
      dplyr::select({{group_col}}, "expression", "density") %>%
      tidyr::unnest(cols = tidyselect::everything())

    # if ggridges requested
    if (use_ggridges) {
      # check to see if ggridges is installed
      rlang::check_installed(pkg = "ggridges")

      if (!requireNamespace(package = "ggridges")) {
        stop("if use_ggridges == TRUE, the ggridges package must be installed")
      }

      # if no group_col is provided, just plot without ggridges
      if (missing(group_col)) {
        return(
          tof_plot_cells_density(
            tof_tibble,
            marker_col = {{marker_col}},
            num_points = num_points,
            theme = theme,
            use_ggridges = FALSE,
            scale = scale,
            ...
          )
        )
      }

      result <-
        ggplot2::ggplot(
          ggplot2::aes(
            x = expression,
            y = {{group_col}},
            fill = {{group_col}},
            height = density
          ),
          data = marker_tibble
        ) +
        ggridges::geom_ridgeline(scale = scale, ...)

    } else {
      # no ggridges
      result <-
        ggplot2::ggplot(
          ggplot2::aes(x = expression, y = density),
          data = marker_tibble
        ) +
        ggplot2::geom_line()

      if (!missing(group_col)) {
        result <-
          result +
          ggplot2::facet_grid(
            rows = ggplot2::vars({{group_col}}),
            scales = "free_y"
          )
      }
    }

    # add x axis label
    result <-
      result +
      ggplot2::labs(x = paste0(marker_colname, " expression"))

    return(result + theme)
  }



#' Plot scatterplots of single-cell data.
#'
#' This function makes scatterplots of single-cell data using user-specified
#' x- and y-axes. Additionally, each point in the scatterplot can be colored
#' using a user-specified variable.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param x_col An unquoted column name specifying which column in
#' `tof_tibble` should be used as the x-axis.
#'
#' @param y_col An unquoted column name specifying which column in
#' `tof_tibble` should be used as the y-axis.
#'
#' @param color_col An unquoted column name specifying which column in
#' `tof_tibble` should be used to color each point in the scatterplot.
#'
#' @param facet_cols An unquoted column name specifying which column in
#' `tof_tibble` should be used to break the scatterplot into facets using
#' \code{\link[ggplot2]{facet_wrap}}.
#'
#' @param theme A ggplot2 theme to apply to the scatterplot. Defaults to
#' \code{\link[ggplot2]{theme_bw}}.
#'
#' @param ... Optional additional arguments to pass to \code{\link[ggplot2]{geom_point}}
#' if \code{method = "ggplot2"} or \code{\link[scattermore]{geom_scattermore}} if
#' \code{method = "scattermore"}.
#'
#' @param method A string indicating which plotting engine should be used. Valid
#' values include "ggplot2" (the default) and "scattermore" (recommended if more than
#' 100K cells are being plotted). Note that \code{method = "scattermore"} requires the
#' scattermore package to be installed.
#'
#' @return A ggplot object.
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
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 vars
#'
#' @importFrom rlang arg_match
#' @importFrom rlang check_installed
#' @importFrom rlang is_installed
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = c(rnorm(n = 500), rnorm(n = 500, mean = 2)),
#'         cd34 = c(rnorm(n = 500), rnorm(n = 500, mean = 4)),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = c(rep("a", 500), rep("b", 500))
#'     )
#'
#'
#'
tof_plot_cells_scatter <-
  function(
    tof_tibble,
    x_col,
    y_col,
    color_col,
    facet_cols,
    theme = ggplot2::theme_bw(),
    ...,# other arguments to ggplot2::geom_point() or scattermore::geom_scattermore()
    method = c("ggplot2", "scattermore")
  ) {
    # check arguments
    method <- rlang::arg_match(method)

    if (missing(x_col) | missing(y_col)) {
      stop("Both x_col and y_col are required.")
    }

    # create plot tibble for memory efficiency
    if (!missing(facet_cols)) {
      plot_tibble <-
        tof_tibble |>
        dplyr::select({{x_col}}, {{y_col}}, {{color_col}}, {{facet_cols}})
    } else {
      plot_tibble <-
        tof_tibble |>
        dplyr::select({{x_col}}, {{y_col}}, {{color_col}})
    }

    # set shape of points for scatterplot
    if (missing(color_col)) {
      shape = 16
    } else {
      shape = 21
    }

    # create point geom
    if (method == "ggplot2") {
      cell_geom <- ggplot2::geom_point(shape = shape, ...)
    } else if (method == "scattermore") {

      # check for scattermore package
      rlang::check_installed(pkg = "scattermore")

      if (!rlang::is_installed(pkg = "scattermore")) {
        stop("`method = scattermore` requires the scattermore package to be installed.")
      }

      cell_geom <-
        scattermore::geom_scattermore(aes(color = {{color_col}}), ...)
    } else {
      stop("Method must be ggplot2 or scattermore.")
    }

    # create plot
    result <-
      plot_tibble %>%
      ggplot2::ggplot(
        ggplot2::aes(x = {{x_col}}, y = {{y_col}}, fill = {{color_col}})
      ) +
      cell_geom +
      ggplot2::labs(
        x = colnames(dplyr::select(plot_tibble, {{x_col}}))[[1]],
        y = colnames(dplyr::select(plot_tibble, {{y_col}}))[[1]]
      )

    if (!missing(facet_cols)) {
      result <-
        result +
        ggplot2::facet_wrap(facets = ggplot2::vars({{facet_cols}}))
    }

    # return result
    return(result + theme)

  }



#' Plot scatterplots of single-cell data using low-dimensional feature embeddings
#'
#' This function makes scatterplots using single-cell data embedded in a
#' low-dimensional space (such as that generated by
#' \code{\link{tof_reduce_dimensions}}, with each point colored using a
#' user-specified variable.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param embedding_cols  Unquoted column names indicating which columns in
#' `tof_tibble` should be used as the x and y axes of the scatterplot. Supports
#' tidyselect helpers. Must select exactly 2 columns. If not provided, a
#' feature embedding can be computed from scratch using the method provided
#' using the `embedding_method` argument and the
#' \code{\link{tof_reduce_dimensions}} arguments passed to `embedding_args`.
#'
#' @param color_col An unquoted column name specifying which column in
#' `tof_tibble` should be used to color each point in the scatterplot.
#'
#' @param facet_cols An unquoted column name specifying which column in
#' `tof_tibble` should be used to break the scatterplot into facets using
#' \code{\link[ggplot2]{facet_wrap}}.
#'
#' @param compute_embedding_cols Unquoted column names indicating which columns
#' in 'tof_tibble' to use for computing the embeddings with the method specified
#' by `embedding_method`. Defaults to all numeric columns in 'tof_tibble'.
#' Supports tidyselect helpers.
#'
#' @param embedding_method A string indicating which method should be used for
#' the feature embedding (if `embedding_cols` are not provided). Options
#' (which are passed to \code{\link{tof_reduce_dimensions}}) are "pca" (the default),
#' "tsne", and "umap".
#'
#' @param embedding_args Optional additional arguments to pass to
#' \code{\link{tof_reduce_dimensions}}. For example, for `method = "tsne"`, these
#' might include `num_comp`, `perplexity`, and `theta`.
#'
#' @param theme A ggplot2 theme to apply to the scatterplot. Defaults to
#' \code{\link[ggplot2]{theme_bw}}.
#'
#' @param ... Optional additional arguments to pass to
#' \code{\link{tof_plot_cells_scatter}}.
#'
#' @param method A string indicating which plotting engine should be used. Valid
#' values include "ggplot2" (the default) and "scattermore" (recommended if more than
#' 100K cells are being plotted). Note that \code{method = "scattermore"} requires the
#' scattermore package to be installed.
#'
#' @return A ggplot object.
#'
#' @family visualization functions
#'
#' @export
#'
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 theme_bw
#'
#' @importFrom rlang arg_match
#' @importFrom rlang sym
#'
#' @examples
#'
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = c(rnorm(n = 500), rnorm(n = 500, mean = 2)),
#'         cd34 = c(rnorm(n = 500), rnorm(n = 500, mean = 4)),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = c(rep("a", 500), rep("b", 500))
#'     )
#'
#' # embed with pca
#' pca_plot <-
#'     tof_plot_cells_embedding(
#'         tof_tibble = sim_data,
#'         color_col = cd38,
#'         embedding_method = "pca",
#'         compute_embedding_cols = starts_with("cd")
#'     )
#'
#' # embed with tsne
#' tsne_plot <-
#'     tof_plot_cells_embedding(
#'         tof_tibble = sim_data,
#'         color_col = cluster_id,
#'         embedding_method = "tsne",
#'         compute_embedding_cols = starts_with("cd")
#'     )
#'
tof_plot_cells_embedding <-
  function(
    tof_tibble,
    embedding_cols,
    color_col,
    facet_cols,
    compute_embedding_cols = where(tof_is_numeric),
    embedding_method = c("pca", "tsne", "umap"),
    embedding_args = list(), # list of arguments for embedding function
    theme = ggplot2::theme_bw(),
    ...,
    method = c("ggplot2", "scattermore")
  ) {

    # if no embedding_cols are specified, use the embedding_method to compute them
    if (missing(embedding_cols)) {
      # if there's no embedding_method specified, just use PCA (for speed)
      if (identical(embedding_method, c("pca", "tsne", "umap"))) {
        message("No embedding_cols were specified, and no embedding_method was specified.
                Performing PCA as the default dimensionality reduction method.")
      }
      # check embedding_method columns
      embedding_method <- rlang::arg_match(embedding_method)

      # compute de novo embedding if needed
      embed_tibble <-
        do.call(
          what = tof_reduce_dimensions,
          args =
            c(
              list(
                tof_tibble = dplyr::select(tof_tibble, {{compute_embedding_cols}}),
                augment = FALSE,
                method = embedding_method),
              embedding_args
            )
        )

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
        stop("2 embedding columns must be selected.")
      }
    }

    # remove any shared columns between tof_tibble and embed_tibble
    cols_to_remove <- intersect(colnames(embed_tibble), colnames(tof_tibble))
    tof_tibble <-
      tof_tibble |>
      dplyr::select(-any_of(cols_to_remove))

    x_col <- rlang::sym(colnames(embed_tibble)[[1]])
    y_col <- rlang::sym(colnames(embed_tibble)[[2]])

    embed_tibble <- dplyr::bind_cols(embed_tibble, tof_tibble)

    # make plot
    if (!missing(facet_cols)) {
      result <-
        tof_plot_cells_scatter(
          tof_tibble = embed_tibble,
          x_col = {{x_col}},
          y_col = {{y_col}},
          color_col = {{color_col}},
          facet_cols = {{facet_cols}},
          theme = theme,
          ...,
          method = method
        )
    } else {
      result <-
        tof_plot_cells_scatter(
          tof_tibble = embed_tibble,
          x_col = {{x_col}},
          y_col = {{y_col}},
          color_col = {{color_col}},
          theme = theme,
          ...,
          method = method
        )
    }

    return(result)
  }




#' Plot force-directed layouts of single-cell data
#'
#' This function makes force-directed layouts using single-cell data embedded in
#' a 2-dimensional space representing a k-nearest-neighbor graph constructed
#' using cell-to-cell similarities. Each node in the force-directed layout
#' represents a single cell colored using a user-specified variable.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param knn_cols Unquoted column names indicating which columns in `tof_tibble`
#' should be used to compute the cell-to-cell distances used to construct
#' the k-nearest-neighbor graph. Supports tidyselect helpers. Defaults to all
#' numeric columns.
#'
#' @param color_col Unquoted column name indicating which column in `tof_tibble`
#' should be used to color the nodes in the force-directed layout.
#'
#' @param facet_cols Unquoted column names indicating which columns in `tof_tibble`
#' should be used to separate nodes into different force-directed layouts.
#'
#' @param num_neighbors An integer specifying how many neighbors should be used
#' to construct the k-nearest neighbor graph.
#'
#' @param graph_type A string specifying if the k-nearest neighbor graph should
#' be "weighted" (the default) or "unweighted".
#'
#' @param graph_layout A string specifying which algorithm should be used to
#' compute the force-directed layout. Passed to \code{\link[ggraph]{ggraph}}.
#' Defaults to "fr", the Fruchterman-Reingold algorithm. Other examples include
#' "nicely", "gem", "kk", and many others. See
#' \code{\link[ggraph]{layout_tbl_graph_igraph}} for other examples.
#'
#' @param distance_function A string indicating which distance function to use
#' in computing the cell-to-cell distances. Valid options include "euclidean"
#' (the default) and "cosine".
#'
#' @param knn_error A value > 0 used in the k-nearest-neighbor approximation.
#' `knn_error` is an error bound such that the ratio between the a cell's reported
#' ith nearest neighbor and its true ith nearest neighbor will be at most
#' 1 + `knn_error`. Defaults to 0 (exact nearest neighbor calculation).
#' Will generally be between 0 and 1. Larger values will result in more
#' approximate KNN calculations (i.e. a higher likelihood of small errors) but
#' will also decrease the computational time of the algorithm significantly.
#'
#' @param edge_alpha A numeric value between 0 and 1 specifying the transparency
#' of the edges drawn in the force-directed layout. Defaults to 0.25.
#'
#' @param node_size A numeric value specifying the size of the nodes in the
#' force-directed layout. Defaults to 2.
#'
#' @param theme A ggplot2 theme to apply to the force-directed layout.
#' Defaults to \code{\link[ggplot2]{theme_void}}
#'
#' @param ... Optional additional arguments to pass to
#' \code{\link[ggraph]{ggraph}}
#'
#' @return A ggraph/ggplot object.
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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = c(rnorm(n = 500), rnorm(n = 500, mean = 2)),
#'         cd34 = c(rnorm(n = 500), rnorm(n = 500, mean = 4)),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = c(rep("a", 500), rep("b", 500))
#'     )
#'
#' # make a layout colored by a marker
#' layout_cd38 <-
#'     tof_plot_cells_layout(
#'         tof_tibble = sim_data,
#'         color_col = cd38
#'     )
#'
#' # make a layout colored by cluster id
#' layout_cluster <-
#'     tof_plot_cells_layout(
#'         tof_tibble = sim_data,
#'         color_col = cluster_id,
#'     )
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



# cluster-level visualizations -------------------------------------------------

#' Visualize clusters in CyTOF data using a minimum spanning tree (MST).
#'
#' This function plots a minimum-spanning tree using clustered single-cell data
#' in order to summarize cluster-level characteristics. Each node in the MST
#' represents a single cluster colored using a user-specified variable (either
#' continuous or discrete).
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param knn_cols Unquoted column names indicating which columns in `tof_tibble`
#' should be used to compute the cluster-to-cluster distances used to construct
#' the k-nearest-neighbor graph. Supports tidyselect helpers. Defaults to all
#' numeric columns.
#'
#' @param color_col Unquoted column name indicating which column in `tof_tibble`
#' should be used to color the nodes in the MST.
#'
#' @param num_neighbors An integer specifying how many neighbors should be used
#' to construct the k-nearest neighbor graph.
#'
#' @param graph_type A string specifying if the k-nearest neighbor graph should
#' be "weighted" (the default) or "unweighted".
#'
#' @param graph_layout This argument specifies a layout for the MST in one of two ways.
#' Option 1: Provide a string specifying which algorithm should be used to
#' compute the force-directed layout. Passed to \code{\link[ggraph]{ggraph}}.
#' Defaults to "nicely", which tries to automatically select a visually-appealing
#' layout. Other examples include "fr", "gem", "kk", and many others. See
#' \code{\link[ggraph]{layout_tbl_graph_igraph}} for other examples.
#' Option 2: Provide a ggraph object previously generated with this
#' function. The layout used to plot this ggraph object will then be used as a
#' template for the new plot. Using this option, number of clusters (and their
#' labels) must be identical to the template. This option is useful if you want
#' to make multiple plots of the same tof_tibble colored by different protein
#' markers, for example.
#'
#' @param central_tendency_function A function to use for computing the
#' measure of central tendency that will be aggregated from each cluster in
#' cluster_col. Defaults to the median.
#'
#' @param distance_function  A string indicating which distance function to use
#' in computing the cluster-to-clusters distances in constructing the MST.
#' Valid options include "euclidean" (the default) and "cosine".
#'
#' @param knn_error A value > 0 used in the k-nearest-neighbor approximation.
#' `knn_error` is an error bound such that the ratio between the a cell's reported
#' ith nearest neighbor and its true ith nearest neighbor will be at most
#' 1 + `knn_error`. Defaults to 0 (exact nearest neighbor calculation).
#' Will generally be between 0 and 1. Larger values will result in more
#' approximate KNN calculations (i.e. a higher likelihood of small errors) but
#' will also decrease the computational time of the algorithm significantly.
#'
#' @param edge_alpha A numeric value between 0 and 1 specifying the transparency
#' of the edges drawn in the force-directed layout. Defaults to 0.25.
#'
#' @param node_size Either a numeric value specifying the size of the nodes in the
#' MST or the string "cluster_size", in which case the size of the node representing
#' each cluster will be scaled according to the number of cells in that cluster
#' (the default).
#'
#' @param theme A ggplot2 theme to apply to the force-directed layout.
#' Defaults to \code{\link[ggplot2]{theme_void}}
#'
#' @param ... Not currently used.
#'
#' @return A ggraph/ggplot object.
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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE)
#'     )
#'
#' # make a layout colored by a marker
#' layout_cd38 <-
#'     tof_plot_cluster_mst(
#'         tof_tibble = sim_data,
#'         cluster_col = cluster_id,
#'         color_col = cd38
#'     )
#'
#'  # use the same layout as the plot above to color the same
#'  # tree using a different marker
#'  layout_cd45 <-
#'     tof_plot_cluster_mst(
#'         tof_tibble = sim_data,
#'         cluster_col = cluster_id,
#'         color_col = cd45,
#'         graph_layout = layout_cd38
#'     )
#'
tof_plot_cluster_mst <-
  function(
    tof_tibble,
    cluster_col,
    knn_cols = where(tof_is_numeric),
    color_col, # each value in cluster_col must map onto only 1 value of group_cols
    num_neighbors = 5L,
    graph_type = c("unweighted", "weighted"),
    graph_layout = "nicely",
    central_tendency_function = stats::median,
    distance_function = c("euclidean", "cosine"),
    knn_error = 0,
    edge_alpha = 0.4,
    node_size = "cluster_size",
    theme = ggplot2::theme_void(),
    ...
  ) {
    # check arguments ----------------------------------------------------------
    # check distance_function argument
    distance_function <- rlang::arg_match(distance_function)

    # check graph_type argument
    graph_type = rlang::arg_match(graph_type)

    # throw error if color_col is missing
    if (missing(color_col)) {
      stop("color_col must be specified.")
    }

    # summarize the clusters ---------------------------------------------------
    color_vector <- dplyr::pull(tof_tibble, {{color_col}})

    # if color_col is a numeric vector
    if (tof_is_numeric(color_vector)) {
      # use a continuous fill scale
      scale_fill <- ggplot2::scale_fill_viridis_c()

      # compute cluster-wise summary statistics
      cluster_tibble <-
        tof_tibble %>%
        dplyr::select(
          {{ cluster_col }},
          {{ color_col }},
          {{ knn_cols }}
        ) %>%
        # compute one summary statistic for each cluster across all knn_cols
        tof_summarize_clusters(
          cluster_col = {{cluster_col}},
          metacluster_cols = c({{knn_cols}}, {{color_col}}),
          central_tendency_function = central_tendency_function
        )

      # compute the size of each cluster
      cluster_sizes <-
        tof_tibble %>%
        dplyr::count(
          {{cluster_col}},
          name = ".cluster_size"
        )

      # if color_col is a character or factor vector
    } else if (is.character(color_vector) | is.factor(color_vector)) {
      # check that each cluster maps to exactly one color
      cluster_groups <-
        tof_tibble %>%
        dplyr::distinct({{cluster_col}}, {{color_col}}) %>%
        dplyr::count({{cluster_col}})

      if (any(cluster_groups$n > 1)) {
        stop(
          "If color_col is a character vector or factor, each cluster must map to exactly one color (i.e. cluster IDs must be nested within color IDs)"
        )
      } else {
        # use a discrete fill scale
        scale_fill <- ggplot2::scale_fill_discrete()

        # compute summary statistics
        cluster_tibble <-
          tof_tibble %>%
          dplyr::select(
            {{cluster_col}},
            {{color_col}},
            {{knn_cols}}
          ) %>%
          # compute one summary statistic for each cluster across all knn_cols
          # but also hold onto each cluster's color_col for plotting
          tof_summarize_clusters(
            cluster_col = {{cluster_col}},
            metacluster_cols = {{knn_cols}},
            group_cols = {{color_col}},
            central_tendency_function = central_tendency_function
          )


        # compute cluster sizes
        cluster_sizes <-
          tof_tibble %>%
          dplyr::count(
            {{cluster_col}},
            {{color_col}},
            name = ".cluster_size"
          )
      }
    }

    # save the names of the clusters to use for calculating each cluster's KNNs
    knn_names <-
      cluster_tibble %>%
      dplyr::select({{knn_cols}}) %>%
      colnames()

    # add the sizes of each cluster to the summary statistics for each cluster
    cluster_tibble <-
      suppressMessages(
        cluster_tibble %>%
          dplyr::left_join(cluster_sizes)
      )

    # make the knn graph -------------------------------------------------------

    # if graph_layout is a previously-plotted mst
    # extract coordinates for each cluster in the mst
    if (inherits(graph_layout, "ggraph")) {
      # save the names of cluster_col and color_col as strings
      cluster_colname <-
        colnames(dplyr::select(cluster_tibble, {{cluster_col}}))

      if (!(cluster_colname %in% colnames(graph_layout$data))) {
        stop("The original layout must have been computed using the same cluster_col as the new plot")
      }

      color_colname <-
        colnames(dplyr::select(cluster_tibble, {{color_col}}))

      # find columns that are shared between the original layout and cluster_tibble
      common_columns <-
        intersect(colnames(cluster_tibble), colnames(graph_layout$data)) %>%
        purrr::discard(.p = ~ .x %in% c(cluster_colname))

      # join any new columns in the cluster_tibble that weren't in the original

      layout_attributes <-
        attributes(graph_layout$data)

      new_layout <-
        graph_layout$data %>%
        dplyr::select(-dplyr::any_of(common_columns)) %>%
        # join
        dplyr::left_join(cluster_tibble, by = cluster_colname)

      # use the new layout to create a new knn_graph from the old one, plus
      # any new information
      knn_graph <-
        layout_attributes[["graph"]] %>%
        tidygraph::activate("nodes")

      if (color_colname %in% colnames(tidygraph::as_tibble(knn_graph)) &
          color_colname %in% colnames(new_layout)) {
        # avoid duplicating the color_col - introduces a bug in the discrete case
        knn_graph <-
          knn_graph %>%
          tidygraph::select(-{{color_col}})

      }

      # make sure that clusters are encoded as character vectors in both
      # representations
      if (!is.character(dplyr::pull(new_layout, {{cluster_col}}))) {
        new_layout[[cluster_colname]] <- as.character(new_layout[[cluster_colname]])
      }

      graph_cluster_vector <-
        knn_graph %>%
        tidygraph::pull({{cluster_col}})

      if (!is.character(graph_cluster_vector)) {
        knn_graph <-
          knn_graph %>%
          tidygraph::mutate(
            "{{cluster_col}}" := as.character({{cluster_col}})
          )
      }

      knn_graph <-
        knn_graph %>%
        tidygraph::select(
          {{cluster_col}},
        ) %>%
        tidygraph::left_join(
          new_layout %>%
            dplyr::select(
              -"x",
              -"y",
              -".ggraph.index",
              -".ggraph.orig_index",
              -"circular"
            ),
          by = cluster_colname
        )

      graph_layout <-
        knn_graph %>%
        tidygraph::activate("nodes") %>%
        tidygraph::as_tibble() %>%
        dplyr::left_join(
          new_layout %>%
            dplyr::select({{cluster_col}}, "x", "y"),
          by = cluster_colname
        ) %>%
        dplyr::select("x", "y")

    } else {
      #calculate the KNN graph from scratch
      knn_graph <-
        cluster_tibble %>%
        tof_make_knn_graph(
          knn_cols = dplyr::any_of(knn_names),
          num_neighbors = num_neighbors,
          distance_function = distance_function,
          knn_error = knn_error,
          graph_type = graph_type
        )
    }

    # make the mst plot --------------------------------------------------------

    # create the edges depending on whether the graph is weighted or unweighted
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

    # return result
    result <- mst_plot + scale_fill + theme
    return(result)
  }


#' Create a volcano plot from differential expression analysis results
#'
#' This function makes a volcano plot using the results of a differential
#' expression analysis (DEA) produced by one of the `tof_dea_*` verbs. Each
#' point in the volcano plot represents a single cluster-marker pair, colored by
#' significance level and the direction of the marker expression difference.
#'
#' @param dea_result A tibble containing the differential expression analysis (DEA)
#' results produced by one of the members of the `tof_dea_*` function family.
#'
#' @param num_top_pairs An integer representing the number of most significant
#' cluster-marker pairs that should be labeled in the volcano plot.
#'
#' @param alpha A numeric value between 0 and 1 representing the significance
#' level below which a p-value should be considered
#' statistically significant. Defaults to 0.05.
#'
#' @param point_size A numeric value specifying the size of the points in the
#' volcano plot.
#'
#' @param label_size A numeric value specifying the size of the text labeling
#' cluster-marker pairs.
#'
#' @param nudge_x A numeric value specifying how far cluster-marker pair labels
#' should be adjusted to the left (if `nudge_x` is negative) or to the right
#' (if `nudge_x` is positive) to avoid overlap with the plotted points.
#' Passed to  \code{\link[ggplot2]{geom_text}}, and ignored if
#' `use_ggrepel` = TRUE. Defaults to 0.
#'
#' @param nudge_y A numeric value specifying how far cluster-marker pair labels
#' should be adjusted downwards (if `nudge_y` is negative) or upwards
#' (if `nudge_y` is positive) to avoid overlap with the plotted points.
#' Passed to  \code{\link[ggplot2]{geom_text}}, and ignored if
#' `use_ggrepel` = TRUE. Defaults to 0.25.
#'
#' @param increase_color A hex code specifying which fill color should
#' be used for points corresponding to cluster-marker pairs where significant
#' increases were detected.
#'
#' @param decrease_color A hex code specifying which fill color should
#' be used for points corresponding to cluster-marker pairs where significant
#' decreases were detected.
#'
#' @param insignificant_color A hex code specifying which fill color should
#' be used for points corresponding to cluster-marker pairs where no significant
#' differences were detected.
#'
#' @param use_ggrepel A boolean value indicting if
#' \code{\link[ggrepel]{geom_text_repel}} should be used to plot labels for
#' cluster-marker pairs. Defaults to FALSE.
#' If TRUE, the ggrepel package must be installed.
#'
#' @param theme A ggplot2 theme to apply to the volcano plot.
#' Defaults to \code{\link[ggplot2]{theme_bw}}
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr transmute
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_vline
#'
#' @importFrom rlang check_installed
#'
#' @importFrom tidyr drop_na
#' @importFrom tidyr unnest
#'
#' @examples
#'
#' # create a mock differential expression analysis result
#' sim_dea_result <-
#'     dplyr::tibble(
#'         cluster_id = rep(letters, 2),
#'         marker = rep(c("cd45", "cd34"), times = length(letters)),
#'         p_adj = runif(n = 2 * length(letters), min = 0, max = 0.5),
#'         mean_fc = runif(n = 2 * length(letters), min = 0.01, max = 10),
#'         significant = dplyr::if_else(p_adj < 0.05, "*", "")
#'     )
#'
#'  attr(sim_dea_result, which = "dea_method") <- "t_unpaired"
#'
#'  # create the volcano plot
#'  volcano <- tof_plot_cluster_volcano(dea_result = sim_dea_result)
#'
#'
tof_plot_cluster_volcano <-
  function(
    dea_result,
    num_top_pairs = 10L,
    alpha = 0.05,
    point_size = 2,
    label_size = 3,
    nudge_x = 0,
    nudge_y = 0.25,
    increase_color = "#207394",
    decrease_color = "#cd5241",
    insignificant_color = "#cdcdcd",
    use_ggrepel = FALSE,
    theme = ggplot2::theme_bw()
  ) {
    # extract dea method from dea_result object
    dea_method <- attr(dea_result, which = "dea_method")

    # if there are multiple results, plot the omnibus
    if("dea_results" %in% colnames(dea_result)) {
      plot_tibble <-
        dea_result %>%
        dplyr::filter(.data$tested_effect == "omnibus") %>%
        tidyr::unnest(cols = "dea_results")
    } else {
      plot_tibble <- dea_result
    }

    num_top_pairs <- min(num_top_pairs, nrow(tidyr::drop_na(dea_result)))

    cluster_index <-
      switch(
        dea_method,
        "lmm" = 2L,
        "t_unpaired" = 1L,
        "t_paired" = 1L,
        "diffcyt_lmm" = 2L,
        "diffcyt_limma" = 2L
      )

    colnames(plot_tibble)[[cluster_index]] <- "cluster"

    if (dea_method %in% c("lmm", "t_unpaired", "t_paired")) {
      plot_tibble <-
        plot_tibble %>%
        dplyr::transmute(
          .data$cluster,
          .data$marker,
          log2_fc = log(.data$mean_fc, base = 2),
          log_p = -log(.data$p_adj),
          significance = .data$significant,
          direction =
            dplyr::case_when(
              .data$significance != "*" ~ "No change",
              .data$mean_fc > 1         ~ "Increase",
              .data$mean_fc < 1         ~ "Decrease"
            )
        )
    } else if (dea_method %in% "diffcyt_limma") {
      plot_tibble <-
        plot_tibble %>%
        dplyr::transmute(
          .data$cluster,
          .data$marker,
          log2_fc = .data$logFC,
          log_p = -log(.data$p_adj),
          significance = .data$significant,
          direction =
            dplyr::case_when(
              .data$significance != "*" ~ "No change",
              .data$logFC > 0           ~ "Increase",
              .data$logFC < 0           ~ "Decrease"
            )
        )
    } else if (dea_method %in% "diffcyt_lmm") {
      stop("diffcyt doesn't report enough information about model fitting to make a volcano plot when diffcyt_method == \"lmm\". Try using tof_dea_lmm()")
    }

    plot_tibble <-
      plot_tibble %>%
      dplyr::arrange(-.data$log_p) %>%
      dplyr::mutate(label = paste(.data$marker, .data$cluster, sep = "@"))

    volcano_plot <-
      plot_tibble %>%
      ggplot2::ggplot(aes(x = .data$log2_fc, y = .data$log_p, fill = .data$direction)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_hline(yintercept = -log(alpha), linetype = "dashed", color = "red") +
      ggplot2::geom_point(shape = 21, size = point_size)

    if (use_ggrepel) {
      # if ggrepel requested
      # check to see if ggridges is installed
      rlang::check_installed(pkg = "ggridges")

      if (!requireNamespace(package = "ggridges")) {
        stop("if use_ggridges == TRUE, the ggridges package must be installed")
      }
      volcano_plot <-
        volcano_plot +
        ggrepel::geom_text_repel(
          ggplot2::aes(label = .data$label),
          data = dplyr::slice_head(plot_tibble, n = num_top_pairs),
          size = label_size
        )
    } else {
      volcano_plot <-
        volcano_plot +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$label),
          data = dplyr::slice_head(plot_tibble, n = num_top_pairs),
          nudge_x = nudge_x,
          nudge_y = nudge_y,
          size = label_size
        )
    }

    volcano_plot <-
      volcano_plot +
      ggplot2::scale_fill_manual(
        values =
          c(
            "Decrease" = decrease_color,
            "Increase" = increase_color,
            "No change" = insignificant_color
          )
      ) +
      ggplot2::labs(
        x = "log2(Fold-change)",
        y = "-log10(p-value)",
        fill = NULL,
        caption =
          paste0(
            "Labels indicate the ",
            as.character(num_top_pairs),
            " most significant cluster-marker pairs"
          )
      )
    return(volcano_plot + theme)
  }


#' Make a heatmap summarizing cluster marker expression patterns in CyTOF data
#'
#' This function makes a heatmap of cluster-to-cluster marker expression patterns
#' in single-cell data. Markers are plotted along the horizontal (x-) axis of
#' the heatmap and cluster IDs are plotted along the vertical (y-) axis of the
#' heatmap.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param marker_cols Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be plotted along the x-axis of the heatmap.
#' Supports tidyselect helpers.
#'
#' @param central_tendency_function A function to use for computing the
#' measure of central tendency that will be aggregated from each cluster in
#' cluster_col. Defaults to the median.
#'
#' @param scale_markerwise A boolean value indicating if the heatmap should
#' rescale the columns of the heatmap such that the maximum value for each
#' marker is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param scale_clusterwise A boolean value indicating if the heatmap should
#' rescale the rows of the heatmap such that the maximum value for each
#' cluster is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param line_width A numeric value indicating how thick the lines separating
#' the tiles of the heatmap should be. Defaults to 0.25.
#'
#' @param theme A ggplot2 theme to apply to the heatmap.
#' Defaults to \code{\link[ggplot2]{theme_minimal}}
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom ggplot2 theme_minimal
#'
#' @importFrom stats median
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE)
#'     )
#'
#' heatmap <-
#'     tof_plot_cluster_heatmap(
#'         tof_tibble = sim_data,
#'         cluster_col = cluster_id
#'     )
#'
tof_plot_cluster_heatmap <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    scale_markerwise = FALSE,
    scale_clusterwise = FALSE,
    line_width = 0.25,
    theme = ggplot2::theme_minimal()
  ) {

    result <-
      tof_tibble %>%
      tof_plot_heatmap(
        y_col = {{cluster_col}},
        marker_cols = {{marker_cols}},
        central_tendency_function = central_tendency_function,
        scale_markerwise = scale_markerwise,
        scale_ywise = scale_clusterwise,
        line_width = line_width,
        theme = theme
      )

    return(result)

  }


# sample-level visualizations --------------------------------------------------

#' Make a heatmap summarizing sample marker expression patterns in CyTOF data
#'
#' This function makes a heatmap of sample-to-sample marker expression patterns
#' in single-cell data. Markers are plotted along the horizontal (x-) axis of
#' the heatmap and sample IDs are plotted along the vertical (y-) axis of the
#' heatmap.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' stores the ids for the sample to which each cell belongs.
#'
#' @param marker_cols Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be plotted along the x-axis of the heatmap.
#' Supports tidyselect helpers.
#'
#' @param central_tendency_function A function to use for computing the
#' measure of central tendency that will be aggregated from each sample in
#' cluster_col. Defaults to the median.
#'
#' @param scale_markerwise A boolean value indicating if the heatmap should
#' rescale the columns of the heatmap such that the maximum value for each
#' marker is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param scale_samplewise A boolean value indicating if the heatmap should
#' rescale the rows of the heatmap such that the maximum value for each
#' sample is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param line_width A numeric value indicating how thick the lines separating
#' the tiles of the heatmap should be. Defaults to 0.25.
#'
#' @param theme A ggplot2 theme to apply to the heatmap.
#' Defaults to \code{\link[ggplot2]{theme_minimal}}
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom ggplot2 theme_minimal
#'
#' @importFrom stats median
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         sample_id = sample(paste0("sample", 1:5), size = 1000, replace = TRUE)
#'     )
#'
#' heatmap <-
#'     tof_plot_sample_heatmap(
#'         tof_tibble = sim_data,
#'         sample_col = sample_id
#'     )
#'
tof_plot_sample_heatmap <-
  function(
    tof_tibble,
    sample_col,
    marker_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    scale_markerwise = FALSE,
    scale_samplewise = FALSE,
    line_width = 0.25,
    theme = ggplot2::theme_minimal()
  ) {
    result <-
      tof_tibble %>%
      tof_plot_heatmap(
        y_col = {{sample_col}},
        marker_cols = {{marker_cols}},
        central_tendency_function = central_tendency_function,
        scale_markerwise = scale_markerwise,
        scale_ywise = scale_samplewise,
        line_width = line_width,
        theme = theme
      )
    return(result)
  }

#' Make a heatmap summarizing sample marker expression patterns in CyTOF data
#'
#' This function makes a heatmap of sample-to-sample marker expression patterns
#' in single-cell data. Markers are plotted along the horizontal (x-) axis of
#' the heatmap and sample IDs are plotted along the vertical (y-) axis of the
#' heatmap.
#'
#' @param feature_tibble A tbl_df or data.frame of aggregated sample-level features,
#' such as that generated by \code{\link{tof_extract_features}}.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' stores the IDs for each sample. If no sample IDs are present, a numeric ID
#' will be assigned to each row of `feature_tibble` based on its row index.
#'
#' @param feature_cols Unquoted column names indicating which column in `feature_tibble`
#' should be interpreted as features to be plotted along the x-axis of the heatmap.
#' Supports tidyselect helpers.
#'
#' @param scale_featurewise A boolean value indicating if the heatmap should
#' rescale the columns of the heatmap such that the maximum value for each
#' marker is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param scale_samplewise A boolean value indicating if the heatmap should
#' rescale the rows of the heatmap such that the maximum value for each
#' sample is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param line_width A numeric value indicating how thick the lines separating
#' the tiles of the heatmap should be. Defaults to 0.25.
#'
#' @param theme A ggplot2 theme to apply to the heatmap.
#' Defaults to \code{\link[ggplot2]{theme_minimal}}
#'
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_minimal
#'
#' @importFrom rlang quo
#'
#' @examples
#'
#' # simulate single-cell data
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         sample_id = sample(paste0("sample", 1:5), size = 1000, replace = TRUE)
#'     )
#'
#' # extract cluster proportions in each simulated patient
#' feature_data <-
#'     tof_extract_proportion(
#'        tof_tibble = sim_data,
#'        cluster_col = cluster_id,
#'        group_cols = sample_id
#'     )
#'
#' # plot the heatmap
#' heatmap <- tof_plot_sample_features(feature_tibble = feature_data)
#'
tof_plot_sample_features <-
  function(
    feature_tibble,
    sample_col,
    feature_cols = where(tof_is_numeric),
    scale_featurewise = FALSE,
    scale_samplewise = FALSE,
    line_width = 0.25,
    theme = ggplot2::theme_minimal()
  ) {
    if (missing(sample_col)) {
      feature_tibble$sample_id <- paste0("sample_", 1:nrow(feature_tibble))
      sample_id <- NULL
      sample_col <- rlang::quo(sample_id)
    }

    result <-
      feature_tibble %>%
      tof_plot_heatmap(
        y_col = {{sample_col}},
        marker_cols = {{feature_cols}},
        scale_markerwise = scale_featurewise,
        scale_ywise = scale_samplewise,
        line_width = line_width,
        theme = theme
      ) +
      ggplot2::labs(x = "feature")

    return(result)
  }

#' Plot the results of a glmnet model fit on sample-level data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which a plot should be made.
#' If new_data isn't provided, the plot will be made using the training data used to
#' fit the model. Alternatively, the string "tuning_data" can be provided, and the
#' plot will be generated using the predictions generated during model tuning.
#'
#' @param theme A ggplot2 theme to apply to the plot
#' Defaults to \code{\link[ggplot2]{theme_bw}}
#'
#' @return A ggplot object. If the `tof_model` is a linear model, a scatterplot
#' of the predicted outcome vs. the true outcome will be returned. If the `tof_model`
#' is a two-class model, an ROC curve will be returned. If the `tof_model` is a
#' multiclass model, a one-versus-all ROC curve will be returned for each class.
#' If `tof_model` is a survival model, a Kaplan-Meier curve will be returned.
#'
#' @export
#'
#' @importFrom ggplot2 theme_bw
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             )
#'    )
#'
#' new_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:20),
#'         cd45 = runif(n = 20),
#'         pstat5 = runif(n = 20),
#'         cd34 = runif(n = 20),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(20),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             )
#'    )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'    )
#'
#' # make the plot
#' plot_1 <- tof_plot_model(tof_model = regression_model, new_data = new_tibble)
#'
#' # train a logistic regression classifier
#' logistic_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = class,
#'         model_type = "two-class"
#'     )
#'
#' # make the plot
#'
#' plot_2 <- tof_plot_model(tof_model = logistic_model, new_data = new_tibble)
#'
tof_plot_model <-
  function(
    tof_model,
    new_data,
    theme = ggplot2::theme_bw()
  ) {

    # check that the tof_model is a tof_model
    if (!inherits(tof_model, "tof_model")) {
      stop("the input `tof_model` must be a tof_model object")
    }

    # find model type from model_fit
    model_type <- tof_get_model_type(tof_model)

    # if new_data is not provided, use training data
    if (missing(new_data)) {
      new_data <- tof_get_model_training_data(tof_model)
    }

    # make plot depending on the input model_fit
    if (model_type == "linear") {
      # make scatterplot of real y values vs. predictions
      result <-
        tof_plot_model_linear(
          tof_model = tof_model,
          new_data = new_data,
          theme = theme
        )

    } else if (model_type == "two-class") {
      # make an ROC curve
      result <-
        tof_plot_model_logistic(
          tof_model = tof_model,
          new_data = new_data,
          theme = theme
        )

    } else if (model_type == "multiclass") {
      # make an ROC curve for each class
      result <-
        tof_plot_model_multinomial(
          tof_model = tof_model,
          new_data = new_data,
          theme = theme
        )
    } else {
      # make a survival curve using the optimal split point
      result <-
        tof_plot_model_survival(
          tof_model = tof_model,
          new_data = new_data,
          theme = theme
        )
    }
    return(result)
  }

#' Plot the results of a linear glmnet model fit on sample-level data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which a plot should be made.
#' If new_data isn't provided, the plot will be made using the training data used to
#' fit the model. Alternatively, the string "tuning_data" can be provided, and the
#' plot will be generated using the predictions generated during model tuning.
#'
#' @param theme A ggplot2 theme to apply to the plot
#' Defaults to \code{\link[ggplot2]{theme_bw}}
#'
#' @return A ggplot object. Specifically, a scatterplot
#' of the predicted outcome vs. the true outcome will be returned.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr tibble
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#' @importFrom stats cor
#'
#' @importFrom tidyr unnest
#'
tof_plot_model_linear <-
  function(tof_model, new_data, theme = ggplot2::theme_bw()) {
    if (is.character(new_data)) {
      if (new_data == "tuning") {
        plot_df <-
          tof_model$tuning_metrics

        plot_df <-
          plot_df %>%
          tidyr::unnest(cols = ".predictions") %>%
          dplyr::mutate(predictions = .data$response)
      }
    } else {
      predictions <-
        tof_predict(
          tof_model = tof_model,
          new_data = new_data,
          prediction_type = "response"
        )

      plot_df <-
        dplyr::tibble(
          truth = new_data[[tof_model$outcome_colnames]],
          predictions = predictions$.pred
        )
    }

    correlation <-
      stats::cor(plot_df$truth, plot_df$predictions) %>%
      round(3)

    result <-
      plot_df %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$truth, y = .data$predictions)) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, ) +
      ggplot2::geom_point() +
      theme +
      ggplot2::labs(caption = paste0("Correlation = ", correlation))

    return(result)
  }

#' Plot the results of a two-class glmnet model fit on sample-level data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which a plot should be made.
#' If new_data isn't provided, the plot will be made using the training data used to
#' fit the model. Alternatively, the string "tuning_data" can be provided, and the
#' plot will be generated using the predictions generated during model tuning.
#'
#' @param theme A ggplot2 theme to apply to the plot.
#' Defaults to \code{\link[ggplot2]{theme_bw}}
#'
#' @return A ggplot object. Specifically, an ROC curve..
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_equal
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
tof_plot_model_logistic <-
  function(tof_model, new_data, theme = ggplot2::theme_bw()) {

    assessment <-
      tof_assess_model(tof_model = tof_model, new_data = new_data)

    roc_auc <-
      assessment$model_metrics %>%
      dplyr::filter(.data$metric == "roc_auc") %>%
      dplyr::pull(.data$value) %>%
      round(3)

    result <-
      assessment$roc_curve %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
      ggplot2::geom_abline(slope = 1, linetype = "dotted", alpha = 0.8) +
      ggplot2::geom_path() +
      ggplot2::coord_equal() +
      theme +
      ggplot2::labs(caption = paste0("AUC = ", roc_auc))

    return(result)
  }


#' Plot the results of a multiclass glmnet model fit on sample-level data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which a plot should be made.
#' If new_data isn't provided, the plot will be made using the training data used to
#' fit the model. Alternatively, the string "tuning_data" can be provided, and the
#' plot will be generated using the predictions generated during model tuning.
#'
#' @param theme A ggplot2 theme to apply to the plot.
#' Defaults to \code{\link[ggplot2]{theme_bw}}.
#'
#' @return A ggplot object. Specifically, a one-versus-all ROC curve
#' (one for each class).
#'
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_equal
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#'
tof_plot_model_multinomial <-
  function(tof_model, new_data, theme = ggplot2::theme_bw()) {

    assessment <-
      tof_assess_model(tof_model = tof_model, new_data = new_data)

    roc_auc <-
      assessment$model_metrics %>%
      dplyr::filter(.data$metric == "roc_auc") %>%
      dplyr::pull(.data$value) %>%
      round(3)

    result <-
      assessment$roc_curve %>%
      ggplot2::ggplot(ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
      ggplot2::geom_abline(slope = 1, linetype = "dotted", alpha = 0.8) +
      ggplot2::geom_path() +
      ggplot2::coord_equal() +
      ggplot2::facet_wrap(facets = ggplot2::vars(.data$.level)) +
      theme +
      ggplot2::labs(caption = paste0("Hand-Till AUC = ", roc_auc))

    return(result)
  }

#' Plot the results of a survival glmnet model fit on sample-level data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which a plot should be made.
#' If new_data isn't provided, the plot will be made using the training data used to
#' fit the model. Alternatively, the string "tuning_data" can be provided, and the
#' plot will be generated using the predictions generated during model tuning.
#'
#' @param theme A ggplot2 theme to apply to the plot.
#' Defaults to \code{\link[ggplot2]{theme_bw}}
#'
#' @param censor_size A numeric value indicating how large to plot the tick marks
#' representing censored values in the Kaplan-Meier curve.
#'
#' @return A ggplot object. Specifically, a Kaplan-Meier curve.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_step
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#' @importFrom purrr map
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
tof_plot_model_survival <-
  function(tof_model, new_data, censor_size = 2.5, theme = ggplot2::theme_bw()) {
    assessment <-
      tof_assess_model(tof_model = tof_model, new_data = new_data)

    p_value <-
      assessment$model_metrics %>%
      dplyr::filter(.data$metric == "log_rank_p_value") %>%
      dplyr::pull(.data$value) %>%
      round(3)

    km_curves <-
      assessment$survival_curves %>%
      dplyr::group_by(.data$risk_group) %>%
      tidyr::nest() %>%
      dplyr::mutate(km_curves = purrr::map(.x = data, .f = tof_compute_km_curve)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"data") %>%
      tidyr::unnest(cols = "km_curves")

    censor_dat <-
      km_curves %>%
      dplyr::filter(.data$is_censored)

    result <-
      km_curves %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data$time_to_event,
          y = .data$survival_probability,
          color = .data$risk_group
        )
      ) +
      ggplot2::geom_step() +
      ggplot2::geom_point(shape = "|", size = censor_size, data = censor_dat) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      theme +
      ggplot2::labs(
        x = "Time to event",
        y = "Survival probability",
        color = "Risk group",
        caption = paste0("log-rank test p-value: ", p_value)
      )

    return(result)

  }


#' Generate a color palette using tidytof.
#'
#' This function generates a color palette based on the color palette of the
#' author's favorite pokemon.
#'
#' @param num_colors An integer specifying the number of colors you'd like to generate.
#'
#' @return A character vector of hex codes specifying the colors in the palette.
#'
#' @export
#'
#' @examples
#' tof_generate_palette(num_colors = 5L)
#'
tof_generate_palette <- function(num_colors) {
  charizard <-
    c(
      "#D86020", "#28A8B8", "#F89040", "#D0D0D0", "#903000",
      "#184068", "#E85040", "#F8D068", "#F8E098",
      "#207890", "#F8A058", "#C03020", "#F8C060", "#F8F8F8"
    )
  result <- charizard[1:num_colors]
  return(result)
}

