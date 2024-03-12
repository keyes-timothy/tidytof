# batch_correction.R

# This file contains functions relevant to performing simple batch correction procedures
# for high-dimensional cytometry data.


#' Batch-correct a tibble of high-dimensional cytometry data using quantile
#' normalization.
#'
#' This function performs quantile normalization on high-dimensional cytometry
#' data in tidy format using \code{\link[preprocessCore]{normalize.quantiles}}.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#'
#' @param augment A boolean value indicating if the output should replace the
#' `channel_cols` in `tof_tibble` with the new, batch corrected columns (TRUE, the default)
#' or if it should only return the batch-corrected columns (FALSE) with all other columns
#' omitted.
#'
#' @return If augment = TRUE, a tibble with the same number of rows and columns as
#' tof_tibble, with the columns specified by `channel_cols` batch-corrected. If
#' augment = FALSE, a tibble containing only the batch-corrected `channel_cols`.
#'
#' @importFrom dplyr all_of
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#'
#' @examples
#' NULL
#'
tof_batch_correct_quantile_tibble <-
    function(tof_tibble, channel_cols, augment = TRUE) {
        # check to see if preprocessCore is installed
        rlang::check_installed(pkg = "preprocessCore")

        if (!requireNamespace(package = "preprocessCore")) {
            stop("Quantile normalization requires the preprocessCore package to be installed from Bioconductor.")
        }

        col_order <- colnames(tof_tibble)
        channel_names <-
            tof_tibble |>
            dplyr::select({{ channel_cols }}) |>
            colnames()

        unchanged_cols <-
            tof_tibble |>
            dplyr::select(-{{ channel_cols }})

        result <-
            tof_tibble |>
            dplyr::select({{ channel_cols }}) |>
            as.matrix() |>
            t() |>
            preprocessCore::normalize.quantiles(
                copy = FALSE
            ) |>
            t() |>
            dplyr::as_tibble()

        colnames(result) <- channel_names

        if (augment) {
            result <-
                dplyr::bind_cols(result, unchanged_cols) |>
                dplyr::relocate(dplyr::all_of(col_order))
        }

        return(result)
    }




#' Batch-correct a tibble of high-dimensional cytometry data using quantile
#' normalization.
#'
#' This function performs quantile normalization on high-dimensional cytometry
#' data in tidy format using \code{\link[preprocessCore]{normalize.quantiles}}.
#' Optionally, groups can be specified and normalized separately.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#'
#' @param group_cols Optional. Unquoted column names indicating which columns
#' should be used to group cells before batch correction. Batch correction is then
#' performed independently within each group. Supports tidyselect helpers.
#'
#' @param augment A boolean value indicating if the output should replace the
#' `channel_cols` in `tof_tibble` with the new, batch corrected columns (TRUE, the default)
#' or if it should only return the batch-corrected columns (FALSE) with all other columns
#' omitted.
#'
#' @return If augment = TRUE, a tibble with the same number of rows and columns as
#' tof_tibble, with the columns specified by `channel_cols` batch-corrected. If
#' augment = FALSE, a tibble containing only the batch-corrected `channel_cols`.
#'
#' @export
#'
#' @importFrom dplyr mutate
#'
#' @importFrom purrr map
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @examples
#' NULL
#'
tof_batch_correct_quantile <-
    function(tof_tibble, channel_cols, group_cols, augment = TRUE) {
        # check to see if preprocessCore is installed
        rlang::check_installed(pkg = "preprocessCore")

        if (!requireNamespace(package = "preprocessCore")) {
            stop("Quantile normalization requires the preprocessCore package to be installed from Bioconductor.")
        }

        if (missing(group_cols)) {
            result <-
                tof_tibble |>
                tof_batch_correct_quantile_tibble(
                    channel_cols = {{ channel_cols }},
                    augment = augment
                )
        } else {
            result <-
                tof_tibble |>
                tidyr::nest(
                    data = -{{ group_cols }}
                ) |>
                dplyr::mutate(
                    data =
                        purrr::map(
                            .x = .data$data,
                            .f =
                                \(x) {
                                    tof_batch_correct_quantile_tibble(
                                        tof_tibble = x,
                                        channel_cols = {{ channel_cols }},
                                        augment = augment
                                    )
                                }
                        )
                ) |>
                tidyr::unnest(cols = .data$data)
        }

        return(result)
    }





#' Perform groupwise linear rescaling of high-dimensional cytometry measurements
#'
#' This function performs quantile normalization on high-dimensional cytometry
#' data in tidy format using linear rescaling. Each channel specified by
#' `channel_cols` is rescaled such that the maximum value is 1 and the minimum
#' value is 0. `group_cols` specifies the columns that should be used to break cells
#' into groups in which the rescaling should be performed separately.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#'
#' @param group_cols Optional. Unquoted column names indicating which columns
#' should be used to group cells before batch correction. Batch correction is then
#' performed independently within each group. Supports tidyselect helpers.
#'
#' @param augment A boolean value indicating if the output should replace the
#' `channel_cols` in `tof_tibble` with the new, batch corrected columns (TRUE, the default)
#' or if it should only return the batch-corrected columns (FALSE) with all other columns
#' omitted.
#'
#' @return If augment = TRUE, a tibble with the same number of rows and columns as
#' tof_tibble, with the columns specified by `channel_cols` batch-corrected. If
#' augment = FALSE, a tibble containing only the batch-corrected `channel_cols`.
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#'
#' @examples
#' NULL
#'
tof_batch_correct_rescale <-
    function(tof_tibble, channel_cols, group_cols, augment = TRUE) {
        result <-
            tof_tibble |>
            dplyr::group_by({{ group_cols }}) |>
            dplyr::mutate(
                dplyr::across(
                    .cols = {{ channel_cols }},
                    tof_rescale
                )
            ) |>
            dplyr::ungroup()

        if (!augment) {
            result <-
                result |>
                dplyr::select({{ group_cols }}, {{ channel_cols }})
        }

        return(result)
    }




#' Perform groupwise linear rescaling of high-dimensional cytometry measurements
#'
#' This function performs quantile normalization on high-dimensional cytometry
#' data in tidy format using either linear rescaling or quantile normalization.
#' Each channel specified by `channel_cols` is batch corrected, and `group_cols`
#' can be used to break cells
#' into groups for which the batch correction should be performed separately.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#'
#' @param group_cols Optional. Unquoted column names indicating which columns
#' should be used to group cells before batch correction. Batch correction is then
#' performed independently within each group. Supports tidyselect helpers.
#'
#' @param augment A boolean value indicating if the output should replace the
#' `channel_cols` in `tof_tibble` with the new, batch corrected columns (TRUE, the default)
#' or if it should only return the batch-corrected columns (FALSE) with all other columns
#' omitted.
#'
#' @param method A string indicating which batch correction method should be used.
#' Valid options are "rescale" for linear scaling (the default) and "quantile"
#' for quantile normalization using \code{\link[preprocessCore]{normalize.quantiles}}.
#'
#' @return If augment = TRUE, a tibble with the same number of rows and columns as
#' tof_tibble, with the columns specified by `channel_cols` batch-corrected. If
#' augment = FALSE, a tibble containing only the batch-corrected `channel_cols`.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_batch_correct <-
    function(
        tof_tibble,
        channel_cols,
        group_cols,
        augment = TRUE,
        method = c("rescale", "quantile")) {
        # check method

        # perform batch correction
        if (method == "rescale") {
            result <-
                tof_tibble |>
                tof_batch_correct_rescale(
                    channel_cols = {{ channel_cols }},
                    group_cols = {{ group_cols }},
                    augment = augment
                )
        } else if (method == "quantile") {
            result <-
                tof_tibble |>
                tof_batch_correct_quantile(
                    channel_cols = {{ channel_cols }},
                    group_cols = {{ group_cols }},
                    augment = augment
                )
        } else {
            stop("Not a valid method")
        }

        return(result)
    }
