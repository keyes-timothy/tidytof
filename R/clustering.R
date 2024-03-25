# clustering.R
# This file contains functions relevant to performing single-cell clustering
# on tof_tibble objects containing high-dimensional cytometry data.

# tof_cluster_flowsom ----------------------------------------------------------

#' Perform FlowSOM clustering on high-dimensional cytometry data
#'
#' This function performs FlowSOM clustering on high-dimensional cytometry data using a user-specified
#' selection of input variables/high-dimensional cytometry measurements. It is mostly a convenient
#' wrapper around \code{\link[FlowSOM]{SOM}} and \code{\link[FlowSOM]{MetaClustering}}.
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
#' @param ... Optional additional parameters that can be passed to the \code{\link[FlowSOM]{BuildSOM}}
#' function.
#'
#' @return A tibble with one column named `.flowsom_cluster` or `.flowsom_metacluster`
#' depending on the value of `perform_metaclustering`. The column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the flowSOM cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @family clustering functions
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 200),
#'         cd38 = rnorm(n = 200),
#'         cd34 = rnorm(n = 200),
#'         cd19 = rnorm(n = 200)
#'     )
#'
#' tof_cluster_flowsom(tof_tibble = sim_data, cluster_cols = c(cd45, cd19))
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
        ...) {
        # check that flowSOM is installed
        has_flowsom <- requireNamespace(package = "FlowSOM", quietly = TRUE)
        if (!has_flowsom) {
          stop(
            "This function requires the {FlowSOM} package. Install it with this code:\n
            if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\n
                install.packages(\"BiocManager\")\n
            }\n
            BiocManager::install(\"FlowSOM\")\n"
          )
        }

        som_distance_function <-
            match.arg(
                arg = som_distance_function,
                choices = c("euclidean", "manhattan", "chebyshev", "cosine")
            )

        # extract string indicating which markers should be used for clustering
        clustering_markers <-
            tof_tibble |>
            dplyr::select({{ cluster_cols }}) |>
            colnames()

        # convert character distance function name to a number that SOM understands
        distf <-
            switch(som_distance_function,
                manhattan = 1,
                euclidean = 2,
                chebyshev = 3,
                cosine = 4
            )

        som <-
            suppressMessages(
                FlowSOM::SOM(
                    data = as.matrix(tof_tibble[, clustering_markers]),
                    xdim = som_xdim,
                    ydim = som_ydim,
                    distf = distf,
                    ...
                )
            )

        # if no metaclustering, return FlowSOM cluster labels
        if (!perform_metaclustering) {
            flowsom_clusters <- as.character(som$mapping[, 1])
            return(dplyr::tibble(.flowsom_cluster = flowsom_clusters))

            # otherwise, perform metaclustering
        } else {
            metacluster_result <-
                FlowSOM::metaClustering_consensus(
                    data = som$codes,
                    k = num_metaclusters
                )

            flowsom_metaclusters <-
                metacluster_result[som$mapping[, 1]] |>
                as.integer() |>
                as.character()

            return(dplyr::tibble(.flowsom_metacluster = flowsom_metaclusters))
        }
    }


# tof_cluster_phenograph -------------------------------------------------------

#' Perform PhenoGraph clustering on high-dimensional cytometry data.
#'
#' This function performs PhenoGraph clustering on high-dimensional cytometry data using a user-specified
#' selection of input variables/high-dimensional cytometry measurements.
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
#'
#'
#' @param ... Optional additional parameters that can be passed to
#' \code{\link{tof_find_knn}}.
#'
#' @return A tibble with one column named `.phenograph_cluster`. This column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the PhenoGraph cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @family clustering functions
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#' tof_cluster_phenograph(tof_tibble = sim_data)
#' tof_cluster_phenograph(tof_tibble = sim_data, cluster_cols = c(cd45, cd19))
#'
tof_cluster_phenograph <-
    function(
        tof_tibble,
        cluster_cols = where(tof_is_numeric),
        num_neighbors = 30,
        distance_function = c("euclidean", "cosine"),
        ...) {
        result <-
            phenograph_cluster(
                tof_tibble,
                cluster_cols = {{ cluster_cols }},
                num_neighbors = num_neighbors,
                distance_function = distance_function,
                ...
            )

        return(result)
    }


# tof_cluster_kmeans -----------------------------------------------------------

#' Perform k-means clustering on high-dimensional cytometry data.
#'
#' This function performs k-means clustering on high-dimensional cytometry data using a user-specified
#' selection of input variables/high-dimensional cytometry measurements. It is mostly a convenient
#' wrapper around \code{\link[stats]{kmeans}}.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the k-means clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param num_clusters An integer indicating the maximum number of clusters
#' that should be returned. Defaults to 20.
#'
#'
#' @param ... Optional additional arguments that can be passed to
#' \code{\link[stats]{kmeans}}.
#'
#' @return A tibble with one column named `.kmeans_cluster`. This column will contain an
#' integer vector of length `nrow(tof_tibble)` indicating the id of
#' the k-means cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#'
#' @family clustering functions
#'
#' @export
#'
#' @importFrom stats kmeans
#' @importFrom purrr pluck
#' @importFrom dplyr tibble
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#' tof_cluster_kmeans(tof_tibble = sim_data)
#' tof_cluster_kmeans(tof_tibble = sim_data, cluster_cols = c(cd45, cd19))
#'
tof_cluster_kmeans <-
    function(
        tof_tibble,
        cluster_cols = where(tof_is_numeric),
        num_clusters = 20,
        ...) {
        kmeans_clusters <-
            stats::kmeans(
                x = select(tof_tibble, {{ cluster_cols }}),
                centers = num_clusters,
                ...
            ) |>
            purrr::pluck("cluster")

        return(dplyr::tibble(.kmeans_cluster = as.character(kmeans_clusters)))
    }



# tof_cluster_ddpr -------------------------------------------------------------

#' Perform developmental clustering on high-dimensional cytometry data.
#'
#' This function performs distance-based clustering on high-dimensional cytometry data
#' by sorting cancer cells (passed into the function as `tof_tibble`) into
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
#' `.\{distance_function\}_cluster`, a character vector of length `nrow(tof_tibble)`
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
#' If `return_distances = FALSE`, a tibble with one column named `.\{distance_function\}_cluster`.
#' This column will contain an integer vector of length `nrow(tof_tibble)` indicating the id of
#' the developmental cluster to which each cell (i.e. each row) in `tof_tibble` was assigned.
#'
#' @family clustering functions
#'
#' @export
#'
#' @importFrom tidyselect all_of
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr rename_with
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#'
#' healthy_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 200),
#'         cd38 = rnorm(n = 200),
#'         cd34 = rnorm(n = 200),
#'         cd19 = rnorm(n = 200),
#'         cluster_id = c(rep("a", times = 100), rep("b", times = 100))
#'     )
#'
#' tof_cluster_ddpr(
#'     tof_tibble = sim_data,
#'     healthy_tibble = healthy_data,
#'     healthy_label_col = cluster_id
#' )
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
        verbose = FALSE) {
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
                healthy_tibble = dplyr::select(healthy_tibble, -{{ healthy_label_col }}),
                healthy_cell_labels = dplyr::pull(healthy_tibble, {{ healthy_label_col }}),
                classifier_markers = {{ cluster_cols }},
                verbose = verbose
            )

        # apply classifier
        if (missing(parallel_cols)) {
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
                    parallel_vars = {{ parallel_cols }}
                )
        }

        # return result
        result <-
            result |>
            dplyr::rename_with(.fn = function(x) paste0(".", x))

        if (!return_distances) {
            result <-
                result |>
                dplyr::select(tidyselect::all_of(paste0(".", distance_function, "_cluster")))
        }

        return(result)
    }



# tof_cluster ------------------------------------------------------------------

#' Cluster high-dimensional cytometry data.
#'
#' This function is a wrapper around tidytof's tof_cluster_* function family.
#' It performs clustering on high-dimensional cytometry data using a user-specified method (of 5 choices)
#' and each method's corresponding input parameters.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_cols Unquoted column names indicating which columns in `tof_tibble` to
#' use in computing the clusters. Defaults to all numeric columns
#' in `tof_tibble`. Supports tidyselect helpers.
#'
#' @param group_cols Optional. Unquoted column names indicating which columns
#' should be used to group cells before clustering. Clustering is then performed
#' on each group independently. Supports tidyselect helpers.
#'
#' @param ... Additional arguments to pass to the `tof_cluster_*`
#' function family member corresponding to the chosen method.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' cluster ids of each cell as a new column in `tof_tibble` (TRUE, the default) or if
#' a single-column tibble including only the cluster ids should be returned (FALSE).
#'
#' @param method A string indicating which clustering methods should be used. Valid
#' values include "flowsom", "phenograph", "kmeans", "ddpr", and "xshift".
#'
#' @return A `tof_tbl` or `tibble` If augment = FALSE, it will have a single column encoding
#' the cluster ids for each cell in `tof_tibble`. If augment = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the cluster ids.
#'
#' @family clustering functions
#'
#' @importFrom dplyr select
#' @importFrom dplyr select
#' @importFrom dplyr select
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 500),
#'         cd38 = rnorm(n = 500),
#'         cd34 = rnorm(n = 500),
#'         cd19 = rnorm(n = 500)
#'     )
#'
#' tof_cluster(tof_tibble = sim_data, method = "kmeans")
#' tof_cluster(tof_tibble = sim_data, method = "phenograph")
#'
tof_cluster <-
    function(
        tof_tibble,
        cluster_cols = where(tof_is_numeric),
        group_cols = NULL,
        ...,
        augment = TRUE,
        method) {
        # find grouping column names
        group_colnames <-
            tof_tibble |>
            dplyr::select({{ group_cols }}) |>
            colnames()

        # if groups are present
        if (length(group_colnames) > 0) {
            result <-
                tof_cluster_grouped(
                    tof_tibble = tof_tibble,
                    group_cols = {{ group_cols }},
                    cluster_cols = {{ cluster_cols }},
                    ...,
                    augment = augment,
                    method = method
                )

            # if groups are not present
        } else {
            result <-
                tof_cluster_tibble(
                    tof_tibble = tof_tibble,
                    cluster_cols = {{ cluster_cols }},
                    ...,
                    augment = augment,
                    method = method
                )
        }
        return(result)
    }

#' Cluster (ungrouped) high-dimensional cytometry data.
#'
#' This function is a wrapper around tidytof's tof_cluster_* function family and
#' provides a low-level API for clustering ungrouped data frames. It is a subroutine
#' of tof_cluster and shouldn't be called directly by users.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param ... Additional arguments to pass to the `tof_cluster_*`
#' function family member corresponding to the chosen method.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' cluster ids of each cell as a new column in `tof_tibble` (TRUE, the default) or if
#' a single-column tibble including only the cluster ids should be returned (FALSE).
#'
#' @param method A string indicating which clustering methods should be used. Valid
#' values include "flowsom", "phenograph", "kmeans", "ddpr", and "xshift".
#'
#' @return A `tof_tbl` or `tibble` If augment = FALSE, it will have a single column encoding
#' the cluster ids for each cell in `tof_tibble`. If augment = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the cluster ids.
#'
#' @importFrom dplyr bind_cols
#'
tof_cluster_tibble <-
    function(
        tof_tibble,
        ...,
        augment = TRUE,
        method) {
        if (method == "flowsom") {
            clusters <-
                tof_tibble |>
                tof_cluster_flowsom(...)
        } else if (method == "phenograph") {
            clusters <-
                tof_tibble |>
                tof_cluster_phenograph(...)
        } else if (method == "kmeans") {
            clusters <-
                tof_tibble |>
                tof_cluster_kmeans(...)
        } else if (method == "ddpr") {
            clusters <-
                tof_tibble |>
                tof_cluster_ddpr(...)
        } else if (method == "xshift") {
            stop("X-shift is not yet implemented.")
        } else {
            stop("Not a valid clustering method.")
        }

        if (augment) {
            result <-
                dplyr::bind_cols(tof_tibble, clusters)
        } else {
            result <- clusters
        }

        return(result)
    }


#' Cluster (grouped) high-dimensional cytometry data.
#'
#' This function is a wrapper around tidytof's tof_cluster_* function family and
#' provides a low-level API for clustering grouped data frames. It is a subroutine
#' of tof_cluster and shouldn't be called directly by users.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param group_cols An unquoted column name indicating which columns
#' should be used to group cells before clustering. Clustering is then performed
#' on each group independently.
#'
#' @param ... Additional arguments to pass to the `tof_cluster_*`
#' function family member corresponding to the chosen method.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' cluster ids of each cell as a new column in `tof_tibble` (TRUE, the default) or if
#' a single-column tibble including only the cluster ids should be returned (FALSE).
#'
#' @param method A string indicating which clustering methods should be used. Valid
#' values include "flowsom", "phenograph", "kmeans", "ddpr", and "xshift".
#'
#' @return A `tof_tbl` or `tibble` If augment = FALSE, it will have a single column encoding
#' the cluster ids for each cell in `tof_tibble`. If augment = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the cluster ids.
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
#' @importFrom purrr map
#' @importFrom purrr map2
#'
#' @importFrom tidyr unite
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
tof_cluster_grouped <-
    function(
        tof_tibble,
        group_cols,
        ...,
        augment = TRUE,
        method) {
        nested_tibble <-
            tof_tibble |>
            dplyr::group_by({{ group_cols }}) |>
            tidyr::nest() |>
            dplyr::ungroup()

        # a list of data frames containing the independently
        # clustered results (appended as columns to the rest of the input tibbles
        # if requested by augment)
        nested_clusters <-
            purrr::map(
                .x = nested_tibble$data,
                .f = tof_cluster_tibble,
                augment = FALSE,
                method = method,
                ...
            )

        # append each group's cluster labels to the group information itself
        # so that cluster labels are distinct among all groups
        result <-
            nested_tibble |>
            tidyr::unite(col = ".prefix", {{ group_cols }}, remove = FALSE) |>
            dplyr::mutate(
                clusters = nested_clusters,
                clusters =
                    purrr::map2(
                        .x = .data$.prefix,
                        .y = .data$clusters,
                        .f = function(.x, .y) {
                            colname <- colnames(.y)
                            new_clusters <-
                                dplyr::tibble(cluster = paste(.y[[1]], .x, sep = "@"))
                            colnames(new_clusters) <- colname
                            return(new_clusters)
                        }
                    )
            ) |>
            dplyr::select(-".prefix")

        if (augment) {
            result <-
                result |>
                #
                tidyr::unnest(cols = c("data", "clusters"))
        } else {
            result <-
                result |>
                dplyr::select("clusters") |>
                #
                tidyr::unnest(cols = "clusters")
        }

        return(result)
    }


# tof_annotate_clusters --------------------------------------------------------

#' Manually annotate tidytof-computed clusters using user-specified labels
#'
#' This function adds an additional column to a `tibble` or `tof_tbl` to allow
#' users to incorporate manual cell type labels for clusters identified using
#' unsupervised algorithms.
#'
#' @param tof_tibble `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' contains the ids of the unsupervised cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param annotations A data structure indicating how to annotate each cluster id
#' in `cluster_col`. `annotations` can be provided as a data.frame with two columns
#' (the first should have the same name as `cluster_col` and contain each unique
#' cluster id; the second can have any name and should contain a character vector
#' indicating which manual annotation should be matched with each cluster
#' id in the first column). `annotations` can also be provided as a named character vector;
#' in this case, each entry in `annotations` should be a unique cluster id, and the
#' names for each entry should be the corresponding manual cluster annotation. See
#' below for examples.
#'
#' @return A `tof_tbl` with the same number of rows as `tof_tibble` and one
#' additional column containing the manual cluster annotations for each cell
#' (as a character vector). If `annotations` was provided as a data.frame, the
#' new column will have the same name as the column containing the cluster annotations
#' in `annotations`. If `annotations` was provided as a named character vector,
#' the new column will be named `\{cluster_col\}_annotation`.
#'
#' @export
#'
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#'
#' @importFrom tibble enframe
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
#' # using named character vector
#' sim_data |>
#'     tof_annotate_clusters(
#'         cluster_col = cluster_id,
#'         annotations = c("macrophage" = "a", "dendritic cell" = "b")
#'     )
#'
#' # using two-column data.frame
#' annotation_data_frame <-
#'     data.frame(
#'         cluster_id = c("a", "b"),
#'         cluster_annotation = c("macrophage", "dendritic cell")
#'     )
#'
#' sim_data |>
#'     tof_annotate_clusters(
#'         cluster_col = cluster_id,
#'         annotations = annotation_data_frame
#'     )
#'
tof_annotate_clusters <- function(tof_tibble, cluster_col, annotations) {
    cluster_colname <-
        tof_tibble |>
        dplyr::select({{ cluster_col }}) |>
        colnames()

    # if annotations are provided as a named vector
    if (is.character(annotations)) {
        annotations <-
            tibble::enframe(
                x = annotations,
                name = paste0(cluster_colname, "_annotation"),
                value = cluster_colname
            )

        # if annotations are provided as a data.frame
    } else if (is.data.frame(annotations)) {
        if (!(cluster_colname %in% colnames(annotations))) {
            stop("One of the columns of `annotations` must have the same name as cluster_col.")
        } else if (ncol(annotations) != 2) {
            stop("`annotations` must have exactly 2 columns.")
        }
    }

    # compute and return result
    result <-
        tof_tibble |>
        dplyr::left_join(
            annotations,
            by = cluster_colname
        )

    return(result)
}
