# metaclustering.R
# This file contains functions relevant to metaclustering clusters within
# tof_tibble objects containing CyTOF data.

# tof_metacluster_hierarchical -------------------------------------------------

#' Metacluster clustered CyTOF data using hierarchical agglomerative clustering
#'
#' This function performs hierarchical metaclustering on a `tof_tbl` containing
#' CyTOF data using a user-specified selection of input variables/CyTOF
#' measurements and
#' the number of desired metaclusters. See \code{\link[stats]{hclust}}.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param num_metaclusters An integer indicating the number of clusters
#' that should be returned. Defaults to 10.
#'
#' @param distance_function A string indicating which distance function should
#' be used to compute the distances between clusters during the hierarchical
#' metaclustering. Options are "euclidean" (the default),
#' "manhattan", "minkowski", "maximum", "canberra", and "binary". See
#' \code{\link[stats]{dist}} for additional details.
#'
#' @param agglomeration_method A string indicating which agglomeration algorithm
#' should be used during hierarchical cluster combination. Options are
#' "complete" (the default), "single", "average", "median", "centroid", "ward.D",
#' "ward.D2", and "mcquitty". See \code{\link[stats]{hclust}} for details.
#'
#' @return A tibble with a single column (`.hierarchical_metacluster`) and
#' the same number of rows as the input `tof_tibble`. Each entry in the column
#' indicates the metacluster label assigned to the same row in `tof_tibble`.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr rename
#' @importFrom dplyr select
#'
#' @importFrom rlang arg_match
#'
#' @importFrom stats dist
#' @importFrom stats cutree
#' @importFrom stats hclust
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
#' tof_metacluster_hierarchical(tof_tibble = sim_data, cluster_col = cluster_id)
#'
tof_metacluster_hierarchical <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    num_metaclusters = 10L,
    distance_function = c("euclidean", "manhattan", "minkowski", "maximum", "canberra", "binary"),
    agglomeration_method = c("complete", "single", "average", "median", "centroid", "ward.D", "ward.D2", "mcquitty")
  ) {
    # check arguments
    distance_function <- rlang::arg_match(distance_function)
    agglomeration_method <- rlang::arg_match(agglomeration_method)

    if (missing(cluster_col)) {
      stop("cluster_col must be specified")
    }

    # extract metacluster colnames
    metacluster_colnames <-
      tof_tibble %>%
      dplyr::select({{metacluster_cols}}) %>%
      colnames()

    # find centroids of all input clusters in tof_tibble
    meta_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = dplyr::any_of(metacluster_colnames),
        central_tendency_function = central_tendency_function
      )

    # distance object
    dist_object <-
      meta_tibble %>%
      dplyr::select(-{{cluster_col}}) %>%
      as.matrix() %>%
      stats::dist(method = distance_function)

    # hierarchical clustering
    hclust_object <-
      stats::hclust(d = dist_object, method = agglomeration_method)
    hclusts <-
      stats::cutree(tree = hclust_object, k = num_metaclusters) %>%
      as.character()

    # return result
    result <-
      tof_tibble %>%
      tof_join_metacluster(
        cluster_col = {{cluster_col}},
        meta_tibble = meta_tibble,
        metacluster_vector = hclusts
      ) %>%
      dplyr::rename(.hierarchical_metacluster = .data$.metacluster)
    return(result)
  }

# tof_metacluster_kmeans--------------------------------------------------------


#' Metacluster clustered CyTOF data using k-means clustering
#'
#' This function performs k-means metaclustering on a `tof_tbl` containing CyTOF data
#' using a user-specified selection of input variables/CyTOF measurements and
#' the number of desired metaclusters. See \code{\link[stats]{hclust}}.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param num_metaclusters An integer indicating the number of clusters
#' that should be returned. Defaults to 10.
#'
#' @param ... Optional additional method specifications to pass to
#' \code{\link{tof_cluster_kmeans}}.
#'
#' @return A tibble with a single column (`.kmeans_metacluster`) and
#' the same number of rows as the input `tof_tibble`. Each entry in the column
#' indicates the metacluster label assigned to the same row in `tof_tibble`.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr rename
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
#' tof_metacluster_kmeans(tof_tibble = sim_data, cluster_col = cluster_id)
#'
tof_metacluster_kmeans <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    num_metaclusters = 10L,
    ...
  ) {

    if (missing(cluster_col)) {
      stop("cluster_col must be specified")
    }

    # extract metacluster colnames
    metacluster_colnames <-
      tof_tibble %>%
      dplyr::select({{metacluster_cols}}) %>%
      colnames()

    # find centroids of all input clusters in tof_tibble
    meta_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = dplyr::any_of(metacluster_colnames),
        central_tendency_function = central_tendency_function
      )

    # k-means clustering
    kmeans_metaclusters <-
      meta_tibble %>%
      tof_cluster_kmeans(
        cluster_cols = dplyr::any_of(metacluster_colnames),
        num_clusters = num_metaclusters,
        ...
      ) %>%
      dplyr::pull(.data$.kmeans_cluster)

    # return result
    result <-
      tof_tibble %>%
      tof_join_metacluster(
        meta_tibble = meta_tibble,
        cluster_col = {{cluster_col}},
        metacluster_vector = kmeans_metaclusters
      ) %>%
      dplyr::rename(.kmeans_metacluster = .data$.metacluster)

    return(result)
  }

# tof_metacluster_phenograph ---------------------------------------------------

#' Metacluster clustered CyTOF data using PhenoGraph clustering
#'
#' This function performs PhenoGraph metaclustering on a `tof_tbl` containing CyTOF data
#' using a user-specified selection of input variables/CyTOF measurements. The number
#' of metaclusters is automatically detected by the PhenoGraph algorithm.
#' See \code{\link{tof_cluster_phenograph}}.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param num_neighbors An integer indicating the number of neighbors to use when
#' constructing PhenoGraph's k-nearest-neighbor graph. Smaller values emphasize
#' local graph structure; larger values emphasize global graph structure (and
#' will add time to the computation). Defaults to 5.
#'
#' @param ... Optional additional method specifications to pass to
#' \code{\link{tof_cluster_phenograph}}.
#'
#' @return A tibble with a single column (`.phenograph_metacluster`) and
#' the same number of rows as the input `tof_tibble`. Each entry in the column
#' indicates the metacluster label assigned to the same row in `tof_tibble`.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr rename
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
#' tof_metacluster_phenograph(tof_tibble = sim_data, cluster_col = cluster_id)
#'
tof_metacluster_phenograph <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    num_neighbors = 5L,
    ...
  ) {

    if (missing(cluster_col)) {
      stop("cluster_col must be specified")
    }

    # extract metacluster colnames
    metacluster_colnames <-
      tof_tibble %>%
      dplyr::select({{metacluster_cols}}) %>%
      colnames()

    # find centroids of all input clusters in tof_tibble
    meta_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = dplyr::any_of(metacluster_colnames),
        central_tendency_function = central_tendency_function
      )

    # phenograph metaclustering
    pheno_metaclusters <-
      meta_tibble %>%
      tof_cluster_phenograph(
        cluster_cols = dplyr::any_of(metacluster_colnames),
        num_neighbors = num_neighbors,
        ...
      ) %>%
      dplyr::pull(.data$.phenograph_cluster)

    # return result
    result <-
      tof_tibble %>%
      tof_join_metacluster(
        meta_tibble = meta_tibble,
        cluster_col = {{cluster_col}},
        metacluster_vector = pheno_metaclusters
      ) %>%
      dplyr::rename(.phenograph_metacluster = .data$.metacluster)

    return(result)

  }

# tof_metacluster_consensus ----------------------------------------------------

#' Metacluster clustered CyTOF data using consensus clustering
#'
#' This function performs consensus metaclustering on a `tof_tbl` containing CyTOF data
#' using a user-specified selection of input variables/CyTOF measurements and
#' the number of desired metaclusters.
#' See \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}} for additional
#' details.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param num_metaclusters An integer indicating the number of clusters
#' that should be returned. Defaults to 10.
#'
#' @param proportion_clusters A numeric value between 0 and 1 indicating the
#' proportion of clusters to subsample (from the total number of clusters in
#' `cluster_col`) during each iteration of the consensus clustering. Defaults
#' to 0.9
#'
#' @param proportion_features A numeric value between 0 and 1 indicating the
#' proportion of features (i.e. the proportion of columns specified by
#' `metacluster_cols`) to subsample during each iteration of the consensus
#' clustering. Defaults to 1 (all features are included).
#'
#' @param num_reps An integer indicating how many subsampled replicates to run
#' during consensus clustering. Defaults to 20.
#'
#' @param clustering_algorithm A string indicating which clustering algorithm
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}} should use to metacluster
#' the subsampled clusters during each resampling. Options are "hierarchical"
#' (the default), "pam" (partitioning around medoids), and "kmeans".
#'
#' @param distance_function A string indicating which distance function should
#' be used to compute the distances between clusters during consensus clustering.
#' Options are "euclidean" (the default),
#' "manhattan", "minkowski", "pearson", "spearman", "maximum", "binary", and
#' "canberra". See \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#'
#' @param ... Optional additional arguments to pass to
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}.
#'
#' @return A tibble with a single column (`.consensus_metacluster`) and
#' the same number of rows as the input `tof_tibble`. Each entry in the column
#' indicates the metacluster label assigned to the same row in `tof_tibble`.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr rename
#'
#' @importFrom rlang arg_match
#' @importFrom rlang check_installed
#' @importFrom rlang is_installed
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
#' tof_metacluster_consensus(tof_tibble = sim_data, cluster_col = cluster_id)
#'
#'
tof_metacluster_consensus <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    num_metaclusters = 10L,
    proportion_clusters = 0.9,
    proportion_features = 1,
    num_reps = 20L,
    clustering_algorithm = c("hierarchical", "pam", "kmeans"),
    distance_function = c("euclidean", "minkowski", "pearson", "spearman", "maximum", "binary", "canberra"),
    ... # optional additional arguments to add to CCP
  ) {

    # check for ConsensusClusterPlus package
    rlang::check_installed(pkg = "ConsensusClusterPlus")

    if (!rlang::is_installed(pkg = "ConsensusClusterPlus")) {
      stop("tof_metacluster_consensus requires the ConsensusClusterPlus package to be installed")
    }

    # check arguments
    distance_function <- rlang::arg_match(distance_function)
    clustering_algorithm <- rlang::arg_match(clustering_algorithm)

    if (missing(cluster_col)) {
      stop("cluster_col must be specified")
    }

    # convert clustering_algorithm string to something CCP will understand
    clustering_algorithm <-
      switch(
        clustering_algorithm,
        "hierarchical" = "hc",
        "pam" = "pam",
        "kmeans" = "km"
      )

    # find centroids of all input clusters in tof_tibble
    meta_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = {{metacluster_cols}},
        central_tendency_function = central_tendency_function
      )

    # data_matrix object
    data_matrix <-
      meta_tibble %>%
      dplyr::select(-{{cluster_col}}) %>%
      as.matrix() %>%
      t()
    colnames(data_matrix) <- dplyr::pull(meta_tibble, {{cluster_col}})

    # consensus clustering
    ## create a temporary file to dump the CCP plot result, which cannot
    ## be disabled, into
    temp_file <- tempdir()
    ccp_object <-
      suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(
          d = data_matrix,
          maxK = num_metaclusters,
          reps = num_reps,
          pItem = proportion_clusters,
          pFeature = proportion_features,
          clusterAlg = clustering_algorithm,
          title = temp_file,
          plot = "pdf",
          distance = distance_function,
          ...
        )
      )
    ccp_metaclusters <-
      ccp_object[[num_metaclusters]]$consensusClass %>%
      as.character()

    # return result
    result <-
      tof_tibble %>%
      tof_join_metacluster(
        cluster_col = {{cluster_col}},
        meta_tibble = meta_tibble,
        metacluster_vector = ccp_metaclusters
      ) %>%
      dplyr::rename(.consensus_metacluster = .data$.metacluster)
    return(result)
  }

# tof_metacluster_flowsom ------------------------------------------------------


#' Metacluster clustered CyTOF data using FlowSOM's built-in metaclustering algorithm
#'
#' This function performs metaclustering on a `tof_tbl` containing CyTOF data
#' using a user-specified selection of input variables/CyTOF measurements and
#' the number of desired metaclusters. It takes advantage of the FlowSOM package's
#' built-in functionality for automatically detecting the number of metaclusters
#' and can use several strategies as adapted by the FlowSOM team: consensus
#' metaclustering, hierarchical metaclustering, k-means metaclustering, or
#' metaclustering using the FlowSOM algorithm itself.
#' See \code{\link[FlowSOM]{MetaClustering}} for additional
#' details.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param num_metaclusters An integer indicating the maximum number of clusters
#' that should be returned. Defaults to 10. Note that for this function, the output
#' may provide a small number of metaclusters than requested. This is because
#' \code{\link[FlowSOM]{MetaClustering}} uses the "Elbow method" to automatically
#' detect the optimal number of metaclusters.
#'
#' @param clustering_algorithm A string indicating which clustering algorithm
#' \code{\link[FlowSOM]{MetaClustering}} should use to perform the metaclustering.
#' Options are "consensus" (the default), "hierarchical", "kmeans", and "som"
#' (i.e. self-organizing map; the FlowSOM algorithm itself).
#'
#'
#' @param ... Optional additional arguments to pass to
#' \code{\link[FlowSOM]{MetaClustering}}.
#'
#' @return A tibble with a single column (`.flowsom_metacluster`) and
#' the same number of rows as the input `tof_tibble`. Each entry in the column
#' indicates the metacluster label assigned to the same row in `tof_tibble`.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr select
#'
#' @importFrom rlang check_installed
#' @importFrom rlang is_installed
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
#' tof_metacluster_flowsom(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     clustering_algorithm = "consensus"
#' )
#'
#' tof_metacluster_flowsom(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     clustering_algorithm = "som"
#' )
#'
#'
tof_metacluster_flowsom <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    num_metaclusters = 10L,
    clustering_algorithm = c("consensus", "hierarchical", "kmeans", "som"),
    ...
  ) {
    # check that FlowSOM is installed
    rlang::check_installed(pkg = "FlowSOM")

    if (!rlang::is_installed(pkg = "FlowSOM")) {
      stop("tof_metacluster_flowsom() requires the FlowSOM package")
    }

    # check arguments
    clustering_algorithm <- rlang::arg_match(clustering_algorithm)

    if (missing(cluster_col)) {
      stop("cluster_col must be specified")
    }

    # convert clustering_algorithm string to something CCP will understand
    clustering_algorithm <-
      switch(
        clustering_algorithm,
        "consensus" = "metaClustering_consensus",
        "hierarchical" = "tof_metaClustering_hclust",
        "kmeans" = "tof_metaClustering_kmeans",
        "som" = "tof_metaClustering_som"
      )

    # find centroids of all input clusters in tof_tibble
    meta_tibble <-
      tof_tibble %>%
      tof_summarize_clusters(
        cluster_col = {{cluster_col}},
        metacluster_cols = {{metacluster_cols}},
        central_tendency_function = central_tendency_function
      )

    # data_matrix
    data_matrix <-
      meta_tibble %>%
      dplyr::select(-{{cluster_col}}) %>%
      as.matrix()

    # perform metaclustering
    flowsom_metaclusters <-
      MetaClustering(
        data = data_matrix,
        method = clustering_algorithm,
        max = num_metaclusters,
        ...
      ) %>%
      as.character()

    # return result
    result <-
      tof_tibble %>%
      tof_join_metacluster(
        meta_tibble = meta_tibble,
        cluster_col = {{cluster_col}},
        metacluster_vector = flowsom_metaclusters
      ) %>%
      dplyr::rename(.flowsom_metacluster = .data$.metacluster)
    return(result)
  }

# tof_metacluster --------------------------------------------------------------

#' Metacluster clustered CyTOF data.
#'
#' This function is a wrapper around {tidytof}'s tof_metacluster_* function family.
#' It performs metaclustering on CyTOF data using a user-specified method (of 5 choices)
#' and each method's corresponding input parameters.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param metacluster_cols Unquoted column names indicating which columns in
#' `tof_tibble` to use in computing the metaclusters.
#' Defaults to all numeric columns in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that should be used to
#' calculate the measurement of central tendency for each cluster before
#' metaclustering. This function will be used to compute a summary statistic for
#' each input cluster in `cluster_col` across all columns specified by
#' `metacluster_cols`, and the resulting vector (one for each cluster) will be
#' used as the input for metaclustering.
#' Defaults to \code{\link[stats]{median}}.
#'
#' @param ... Additional arguments to pass to the `tof_metacluster_*` function
#' family member corresponding to the chosen `method`.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' metacluster ids of each cell as a new column in `tof_tibble` (TRUE; the default) or if
#' a single-column tibble including only the metacluster ids should be returned (FALSE).
#'
#' @param method A string indicating which clustering method should be used. Valid
#' values include "consensus", "hierarchical", "kmeans", "phenograph", and "flowsom".
#'
#' @return A `tof_tbl` or `tibble` If augment = FALSE, it will have a single column encoding
#' the metacluster ids for each cell in `tof_tibble`. If augment = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the metacluster ids.
#'
#' @family metaclustering functions
#'
#' @export
#'
#' @importFrom dplyr bind_cols
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
#' tof_metacluster(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     clustering_algorithm = "consensus",
#'     method = "flowsom"
#' )
#'
#' tof_metacluster(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     method = "phenograph"
#' )
#'
tof_metacluster <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    ...,
    augment = TRUE,
    method = c("consensus", "hierarchical", "kmeans", "phenograph", "flowsom")
  ) {
    # check arguments
    method <- rlang::arg_match(method)

    if (method == "consensus") {
      metaclusters <-
        tof_metacluster_consensus(
          tof_tibble = tof_tibble,
          cluster_col = {{cluster_col}},
          metacluster_cols = {{metacluster_cols}},
          central_tendency_function = central_tendency_function,
          ...
        )
    } else if (method == "hierarchical") {
      metaclusters <-
        tof_metacluster_hierarchical(
          tof_tibble = tof_tibble,
          cluster_col = {{cluster_col}},
          metacluster_cols = {{metacluster_cols}},
          central_tendency_function = central_tendency_function,
          ...
        )
    } else if (method == "kmeans") {
      metaclusters <-
        tof_metacluster_kmeans(
          tof_tibble = tof_tibble,
          cluster_col = {{cluster_col}},
          metacluster_cols = {{metacluster_cols}},
          central_tendency_function = central_tendency_function,
          ...
        )
    } else if (method == "phenograph") {
      metaclusters <-
        tof_metacluster_phenograph(
          tof_tibble = tof_tibble,
          cluster_col = {{cluster_col}},
          metacluster_cols = {{metacluster_cols}},
          central_tendency_function = central_tendency_function,
          ...
        )
    } else if (method == "flowsom") {
      metaclusters <-
        tof_metacluster_flowsom(
          tof_tibble = tof_tibble,
          cluster_col = {{cluster_col}},
          metacluster_cols = {{metacluster_cols}},
          central_tendency_function = central_tendency_function,
          ...
        )
    } else {
      stop("Not a valid metaclustering method.")
    }

    # return result
    if (augment) {
      result <- dplyr::bind_cols(tof_tibble, metaclusters)
    } else {
      result <- metaclusters
    }

    return(result)
  }

