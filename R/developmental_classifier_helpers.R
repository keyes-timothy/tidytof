# developmental_classifier_helpers.R
# This file contains helper functions for performing "developmental classification"
# of cancer cells into their most-similar healthy developmental subpopulation per
# Good et al., 2018 (Nature Medicine).
#
# All functions in this file are subroutines of the function tof_cluster_ddpr
# in clustering.R

#' Calculate centroids and covariance matrices for each cell subpopulation in
#' healthy CyTOF data.
#'
#' This function takes a `tibble` or `tof_tibble` storing healthy cell measurements
#' in each of its rows and a vector (`healthy_cell_labels`) representing the
#' cell subpopulation to which each cell belongs. It uses these values to calculate
#' several values required to perform "developmental classification" as described in
#' \href{https://pubmed.ncbi.nlm.nih.gov/29505032/}{this paper}.
#'
#' @param healthy_tibble A `tibble` or `tof_tibble` containing cells from only
#' healthy control samples (i.e. not disease samples).
#'
#' @param healthy_cell_labels A character or integer vector of length `nrow(healthy_tibble)`.
#' Each entry in this vector should represent the cell subpopulation label (or cluster id) for
#' the corresponding row in `healthy_tibble`.
#'
#' @param classifier_markers Unquoted column names indicating which columns in `healthy_tibble` to
#' use in the developmental classification. Defaults to all numeric columns
#' in `healthy_tibble`. Supports tidyselect helpers.
#'
#' @param verbose A boolean value indicating if updates should be printed to the
#' console during classification. Defaults to FALSE.
#'
#' @return A tibble with three columns:
#' \strong{population} (id of the healthy cell population),
#' \strong{centroid} (the centroid vector for that cell population), and
#' \strong{covariance_matrix} (the covariance matrix for that cell population)
#'
#'
tof_build_classifier <- function(
    healthy_tibble = NULL,
    healthy_cell_labels = NULL,
    classifier_markers = where(tof_is_numeric),
    verbose = FALSE) {
    if (verbose) {
        message("Building classifier from healthy cells...")

        message("Reshaping healthy data...")
    }

    healthy_tibble <-
        healthy_tibble |>
        dplyr::mutate(population = healthy_cell_labels) |>
        dplyr::select({{ classifier_markers }}, "population") |>
        dplyr::group_by(.data$population) |>
        tidyr::nest()

    # Calculate mean and covariance matrix for each population
    if (verbose) {
        message("Calculating mean and covariance matrix for all healthy populations")
    }
    classifier_fit <-
        healthy_tibble |>
        dplyr::transmute(
            .data$population,
            centroid =
                purrr::map(
                    data,
                    ~ dplyr::summarize(.x, dplyr::across(dplyr::everything(), mean)) |>
                        tidyr::pivot_longer(tidyr::everything()) |>
                        deframe()
                ),
            covariance_matrix = purrr::map(data, cov)
        )
    if (verbose) {
        message("Done! Returning classifier_fit object")
    }
    return(classifier_fit)
}




#' Classify each cell (i.e. each row) in a matrix of cancer cells into its most
#' similar healthy developmental subpopulation.
#'
#' This function uses a specified distance metric to classify each cell in a data.frame
#' or matrix (`cancer_data`) into one of `nrow(classifier_fit)` subpopulations
#' based on minimum distance, as described in \href{https://pubmed.ncbi.nlm.nih.gov/29505032/}{this paper}.
#'
#' @param classifier_fit A tibble produced by \code{\link{tof_build_classifier}}.
#'
#' @param cancer_data A matrix in which each row corresponds to a cell and each
#' column corresponds to a measured CyTOF antigen.
#'
#' @param distance_function A string indicating which of three distance functions should
#' be used to calculate the distances between each row of `cancer_data` and the
#' healthy developmental subpopulations corresponding to each row of `classifier_fit`.
#'
#' @return A data.frame in which each column represents the distance between
#' a cell in the input data and each healthy subpopulation cells are being
#' classified into.
#'
tof_classify_cells <-
    function(
        classifier_fit,
        cancer_data,
        distance_function = c("mahalanobis", "cosine", "pearson")) {
        distance_function <-
            match.arg(distance_function, c("mahalanobis", "cosine", "pearson"))

        # Calculate distance to each population
        if (distance_function == "mahalanobis") {
            classifications <-
                classifier_fit |>
                dplyr::transmute(
                    classification_data =
                        purrr::map2(
                            .x = .data$centroid,
                            .y = .data$covariance_matrix,
                            .f = ~ mahalanobis(cancer_data, center = .x, cov = .y)
                        )
                ) |>
                dplyr::pull(.data$classification_data) |>
                purrr::quietly(dplyr::bind_cols)() |>
                purrr::pluck("result")
        } else if (distance_function == "cosine") {
            classifications <-
                classifier_fit |>
                dplyr::transmute(
                    classification_data =
                        purrr::map(
                            .x = .data$centroid,
                            .f =
                                ~ tof_cosine_dist(
                                    matrix = as.matrix(cancer_data),
                                    vector = .x
                                )
                        )
                ) |>
                pull(.data$classification_data) |>
                purrr::quietly(dplyr::bind_cols)() |>
                purrr::pluck("result")
        } else if (distance_function == "pearson") {
            classifications <-
                classifier_fit |>
                dplyr::mutate(
                    classification_data =
                        purrr::map(
                            .x = .data$centroid,
                            .f = ~ 1 - cor(x = t(as.matrix(cancer_data)), y = .x)
                        )
                ) |>
                dplyr::pull(.data$classification_data) |>
                purrr::quietly(dplyr::bind_cols)() |>
                purrr::pluck("result")
        }

        colnames(classifications) <- as.character(classifier_fit$population)
        if (distance_function == "pearson") {
            classifications <- dplyr::as_tibble(as.matrix(classifications))
        }

        population_names <-
            classifications |>
            dplyr::mutate(..cell_id = seq_len(length.out = nrow(classifications))) |>
            tidyr::pivot_longer(
                cols = -"..cell_id",
                names_to = "cluster",
                values_to = "distance"
            ) |>
            dplyr::group_by(.data$..cell_id) |>
            dplyr::filter(.data$distance == min(.data$distance)) |>
            dplyr::select(-"distance")

        classifications <-
            classifications |>
            dplyr::mutate(..cell_id = seq_len(length.out = nrow(classifications))) |>
            dplyr::left_join(population_names, by = "..cell_id") |>
            dplyr::select(-"..cell_id")

        return(classifications)
    }


#' A function for finding the cosine distance between each of the rows of a numeric
#' matrix and a numeric vector.
#'
#' @param matrix A numeric matrix.
#'
#' @param vector A numeric vector.
#'
#' @return A numeric vector of distances of length `nrow(matrix)` in which the
#' ith entry represents the cosine distance between the ith row of `matrix` and
#' `vector`.
#'
#' @examples
#' NULL
tof_cosine_dist <-
    function(matrix, vector) {
        diag_matrix_matrix_t <- rep(0, nrow(matrix))
        for (i in seq_len(nrow(matrix))) {
            diag_matrix_matrix_t[[i]] <- crossprod(matrix[i, ])
        }
        distances <-
            (matrix %*% vector) /
                (as.vector(sqrt(crossprod(vector))) * sqrt(diag_matrix_matrix_t))

        distances <- as.numeric(1 - distances)

        return(distances)
    }


#' Perform developmental clustering on CyTOF data using a pre-fit classifier
#'
#' @param cancer_tibble A `tibble` or `tof_tibble` containing cells to be classified
#' into their nearest healthy subpopulation (generally cancer cells).
#'
#' @param classifier_fit A nested `tibble` produced by `tof_build_classifier` in which
#' each row represents a healthy cell subpopulation into which the cells in `cancer_tibble`
#' should be classified using minimum distance.
#'
#' @param distance_function A string indicating which distance function should
#' be used to perform the classification. Options are "mahalanobis" (the default),
#' "cosine", and "pearson".
#'
#' @param num_cores An integer indicating the number of CPU cores used to parallelize
#' the classification. Defaults to 1 (a single core).
#'
#' @param parallel_vars Unquoted column names indicating which columns in `cancer_tibble` to
#' use for breaking up the data in order to parallelize the classification.
#' Defaults to NULL. Supports tidyselect helpers.
#'
#' @return A tibble with `nrow(cancer_tibble)` rows and `nrow(classifier_fit) + 1`
#' columns. Each row represents a cell from `cancer_tibble`, and `nrow(classifier_fit)`
#' of the columns represent the distance between the cell and each of the healthy
#' subpopulations' cluster centroids. The final column represents the cluster id of
#' the healthy subpopulation with the minimum distance to the cell represented
#' by that row.
#'
#' @examples
#' NULL
#'
tof_apply_classifier <- function(
    cancer_tibble = NULL,
    classifier_fit = NULL,
    distance_function = c("mahalanobis", "cosine", "pearson"),
    num_cores = 1,
    parallel_vars) {
    # check distance function
    distance_function <-
        match.arg(distance_function, c("mahalanobis", "cosine", "pearson"))

    classifier_markers <- colnames(classifier_fit$covariance_matrix[[1]])

    # If no variable over which to parallelize was specified
    if (missing(parallel_vars)) {
        classification_data <-
            tof_classify_cells(
                classifier_fit = classifier_fit,
                cancer_data =
                    dplyr::select(cancer_tibble, dplyr::any_of(classifier_markers)),
                distance_function = distance_function
            ) |>
            dplyr::rename_with(function(x) stringr::str_c(distance_function, x, sep = "_"), .cols = everything())

        # otherwise,
    } else {
        # initialize values related to parallel processing
        `%my_do%` <- `%do%`
        my_cluster <- NULL

        # nest cancer data
        cancer_tibble <-
            cancer_tibble |>
            dplyr::select(
                {{ parallel_vars }},
                dplyr::any_of(classifier_markers)
            ) |>
            mutate(
                ..cell_id = seq_len(nrow(cancer_tibble))
            ) |>
            # check this
            tidyr::nest(data = -{{ parallel_vars }}, ..cell_ids = "..cell_id") |>
            dplyr::ungroup()

        # set up parallel back-end
        if (num_cores != 1) {
            my_cluster <- parallel::makeCluster(num_cores)
            doParallel::registerDoParallel(my_cluster)
            `%my_do%` <- `%dopar%`
        }

        # hack-y way to get around R CMD check note for expression matrix
        # not having a visible global binding due to NSE
        if (is.null(cancer_tibble)) {
            expression_matrix <- NULL
        }

        # cluster each value of parallel_vars on a difference core simultaneously
        classification_data <-
            foreach::foreach(
                expression_matrix = cancer_tibble$data,
                .combine = list,
                .packages = c("dplyr", "purrr", "tidyr"),
                .export = c("tof_classify_cells", "tof_cosine_dist"),
                .multicombine = TRUE,
                .maxcombine = nrow(cancer_tibble)
            ) %my_do%
            tof_classify_cells(
                classifier_fit = classifier_fit,
                cancer_data = dplyr::select(expression_matrix, -"..cell_id"),
                distance_function = distance_function
            )

        # stop cluster if it was set up
        if (!is.null(my_cluster)) {
            parallel::stopCluster(my_cluster)
        }

        # unnest final classification data and return it
        classification_data <-
            dplyr::tibble(
                classification_data = classification_data,
                ..cell_ids = cancer_tibble$..cell_ids
            ) |>
            tidyr::unnest(cols = c("classification_data", "..cell_ids")) |>
            dplyr::arrange(.data$..cell_id) |>
            dplyr::select(-"..cell_id") |>
            dplyr::rename_with(
                function(x) stringr::str_c(distance_function, x, sep = "_"),
                .cols = tidyselect::everything()
            )
    }

    return(classification_data)
}
