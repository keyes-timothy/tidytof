# feature_extraction.R
# This file contains functions relevant to extracting patient- or sample-level
# features by aggregating single-cell data in tof_tibble objects.

# tof_extract_proportion -------------------------------------------------------

#' Extract the proportion of cells in each cluster in a `tof_tibble`.
#'
#' This feature extraction function allows you to calculate the proportion of
#' cells in each cluster in a `tof_tibble` - either overall or when broken down
#' into subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param format A string indicating if the data should be returned in "wide" format
#' (the default; each cluster proportion is given its own column) or in "long" format
#' (each cluster proportion is provided as its own row).
#'
#' @return A tibble.
#'
#' If format == "wide", the tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` and one column for each grouping variable
#' as well as one column for the proportion of cells in each cluster. The names of each
#' column containing cluster proportions is obtained using the following pattern:
#' "prop@\{cluster_id\}".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols` and each cluster id (i.e. level) in `cluster_col`. It will have one column for
#' each grouping variable, one column for the cluster ids, and one column (`prop`) containing the
#' cluster proportions.
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unite
#' @importFrom tidyr separate
#' @importFrom rlang arg_match
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract proportion of each cluster in each patient in wide format
#' tof_extract_proportion(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient
#' )
#'
#' # extract proportion of each cluster in each patient in long format
#' tof_extract_proportion(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     format = "long"
#' )
#'
tof_extract_proportion <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    format = c("wide", "long")
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    group_colnames <-
      tof_tibble |>
      dplyr::select({{group_cols}}) |>
      colnames()

      my_sep <- "________"

      if (length(group_colnames) != 0) {
        abundances <-
          tof_tibble |>
          tidyr::unite(
            col = ".group",
            {{group_cols}},
            sep = my_sep
          ) |>
          dplyr::mutate(.group = as.factor(.data$.group))

      } else {
        group_cols <- NULL

        abundances <-
          tof_tibble |>
          dplyr::mutate(.group = 1)
      }

      abundances <-
        abundances |>
        dplyr::mutate(
          group = as.factor(.data$.group),
          "{{cluster_col}}" := as.factor({{cluster_col}})
        ) |>
        dplyr::count(.data$.group, {{cluster_col}}, .drop = FALSE, name = "abundance") |>
        dplyr::group_by(.data$.group) |>
        dplyr::mutate(
          prop = .data$abundance / sum(.data$abundance),
          "{{cluster_col}}" := as.character({{cluster_col}})
        ) |>
        dplyr::ungroup()

      if (length(group_colnames) != 0) {
        abundances <-
          abundances |>
          tidyr::separate(col = .data$.group, into = group_colnames, sep = my_sep)
      } else {
        abundances <-
          abundances |>
          dplyr::select(-".group")
      }
      abundances <-
        abundances |>
        dplyr::select({{group_cols}}, {{cluster_col}}, "prop")

      if (format == "wide") {
      abundances <-
        abundances |>
        dplyr::filter(.data$prop != 0) |>
        tidyr::pivot_wider(
          names_from = {{cluster_col}},
          values_from = "prop",
          names_prefix = "prop@",
          # note that if a cluster is not present for a given group, it will
          # be filled in as having a relative abundance of 0.
          values_fill = 0
        )
      }

      return(dplyr::ungroup(abundances))
  }

# tof_extract_central_tendency -------------------------------------------------

#' Extract the central tendencies of CyTOF markers in each cluster in a `tof_tibble`.
#'
#' This feature extraction function calculates a user-specified measurement of central tendency
#' (i.e. median or mode) of the cells in each cluster in a `tof_tibble` across a
#' user-specified selection of CyTOF markers. These calculations can be done either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the feature extraction
#' calculation. Defaults to all numeric (integer or double) columns.
#' Supports tidyselection.
#'
#' @param stimulation_col Optional. An unquoted column name that indicates which
#' column in `tof_tibble` contains information about which stimulation condition each cell
#' was exposed to during data acquisition. If provided, the feature extraction will be
#' further broken down into subgroups by stimulation condition (and features from each stimulation
#' condition will be included as their own features in wide format).
#'
#' @param central_tendency_function The function that will be used to calculate
#' the measurement of central tendency for each cluster (to be used
#' as the dependent variable in the linear model). Defaults to \code{\link[stats]{median}}.
#'
#' @param format A string indicating if the data should be returned in "wide" format
#' (the default; each cluster feature is given its own column) or in "long" format
#' (each cluster feature is provided as its own row).
#'
#' @return A tibble.
#'
#' If format == "wide", the tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` and one column for each grouping variable,
#' one column for each extracted feature (the central tendency of a given marker in a given cluster).
#' The names of each column containing cluster features is obtained using the following pattern:
#' "\{marker_id\}@\{cluster_id\}_ct".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#' @importFrom stringr str_c
#' @importFrom stringr str_remove
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract proportion of each cluster in each patient in wide format
#' tof_extract_central_tendency(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient
#' )
#'
#' # extract proportion of each cluster in each patient in long format
#' tof_extract_central_tendency(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     format = "long"
#' )
#'
tof_extract_central_tendency <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    stimulation_col = NULL,
    central_tendency_function = stats::median,
    format = c("wide", "long")
  ) {

    # check format argument
    format <- rlang::arg_match(format)

    central_tendencies <-
      tof_tibble |>
      # if cluster_col isn't a character vector, make it one
      dplyr::mutate("{{cluster_col}}" := as.character({{cluster_col}})) |>
      dplyr::group_by(dplyr::across({{group_cols}}), {{cluster_col}}, {{stimulation_col}}) |>
      dplyr::summarize(dplyr::across({{marker_cols}}, central_tendency_function)) |>
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "channel",
        values_to = "values"
      )

    if (format == "wide") {
      central_tendencies <-
        central_tendencies |>
        dplyr::group_by(dplyr::across({{group_cols}})) |>
        dplyr::transmute(
          col_names =
            stringr::str_c(
              {{stimulation_col}},
              "_",
              .data$channel,
              "@",
              {{cluster_col}},
              "_ct"
            ) |>
            stringr::str_remove("^_"),
          .data$values
        ) |>
        tidyr::pivot_wider(
          names_from = "col_names",
          values_from = "values"
        )
    }

    return(dplyr::ungroup(central_tendencies))
  }




# tof_extract_threshold --------------------------------------------------------

#' Extract aggregated features from CyTOF data using a binary threshold
#'
#' This feature extraction function calculates the proportion of cells in a given cluster
#' that have a CyTOF antigen expression over a user-specified threshold across a
#' user-specified selection of CyTOF markers. These calculations can be done either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the feature extraction
#' calculation. Defaults to all numeric (integer or double) columns.
#' Supports tidyselect helpers.
#'
#' @param stimulation_col Optional. An unquoted column name that indicates which
#' column in `tof_tibble` contains information about which stimulation condition each cell
#' was exposed to during data acquisition. If provided, the feature extraction will be
#' further broken down into subgroups by stimulation condition (and features from each stimulation
#' condition will be included as their own features in wide format).
#'
#' @param threshold A double or integer of length 1 indicating what threshold should be used.
#'
#' @param format A string indicating if the data should be returned in "wide" format
#' (the default; each cluster feature is given its own column) or in "long" format
#' (each cluster feature is provided as its own row).
#'
#' @return A tibble.
#'
#' If format == "wide", the tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` and one column for each grouping variable,
#' one column for each extracted feature (the proportion of cells in a given cluster
#' over with marker expression values over `threshold`).
#' The names of each column containing cluster features is obtained using the following pattern:
#' "\{marker_id\}@\{cluster_id\}_threshold".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_c
#' @importFrom stringr str_remove
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract proportion of each cluster in each patient in wide format
#' tof_extract_threshold(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient
#' )
#'
#' # extract proportion of each cluster in each patient in long format
#' tof_extract_threshold(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     format = "long"
#' )
#'
tof_extract_threshold <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    stimulation_col = NULL,
    threshold = asinh(10 / 5),
    format = c("wide", "long")
  ) {

    # check format argument
    format <- rlang::arg_match(format)

    threshold_features <-
      tof_tibble |>
      # if the cluster column isn't a character vector, make it one
      dplyr::mutate("{{cluster_col}}" := as.character({{cluster_col}})) |>
      dplyr::group_by(
        dplyr::across({{group_cols}}), {{cluster_col}}, {{stimulation_col}}
      ) |>
      dplyr::summarize(
        dplyr::across({{marker_cols}}, ~ mean(.x > threshold))
      ) |>
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "channel",
        values_to = "values"
      )

    if (format == "wide") {
      threshold_features <-
        threshold_features |>
        dplyr::group_by(dplyr::across({{group_cols}})) |>
        dplyr::transmute(
          col_names =
            stringr::str_c(
              {{stimulation_col}},
              "_",
              .data$channel,
              "@",
              {{cluster_col}},
              "_threshold"
            ) |>
            stringr::str_remove("^_"),
          .data$values
        ) |>
        tidyr::pivot_wider(
          names_from = "col_names",
          values_from = "values"
        )
    }

    return(dplyr::ungroup(threshold_features))
  }



# tof_extract_emd --------------------------------------------------------------

#' Extract aggregated features from CyTOF data using earth-mover's distance (EMD)
#'
#' This feature extraction function calculates the earth-mover's distance (EMD) between
#' the stimulated and unstimulated ("basal") experimental conditions of samples in a
#' CyTOF experiment. This calculation is performed across a
#' user-specified selection of CyTOF antigens and can be performed either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the earth-mover's distance
#' calculation. Defaults to all numeric (integer or double) columns.
#' Supports tidyselect helpers.
#'
#' @param emd_col An unquoted column name that indicates which
#' column in `tof_tibble` should be used to group cells into different distributions
#' to be compared with one another during the EMD calculation. For example,
#' if you want to compare marker expression distributions across stimulation
#' conditions, `emd_col` should be the column in `tof_tibble` containing
#' information about which stimulation condition each cell
#' was exposed to during data acquisition.
#'
#' If provided, the feature extraction will be
#' further broken down into subgroups by stimulation condition (and features from each stimulation
#' condition will be included as their own features in wide format).
#'
#' @param reference_level A string indicating what the value in `emd_col`
#' corresponds to the "reference" value to which all other values in `emd_col`
#' should be compared. For example, if `emd_col` represents the stimulation
#' condition for a cell, reference_level might take the value of "basal" or
#' "unstimulated" if you want to compare each stimulation to the basal state.
#'
#' @param format A string indicating if the data should be returned in "wide" format
#' (the default; each cluster feature is given its own column) or in "long" format
#' (each cluster feature is provided as its own row).
#'
#' @param num_bins Optional. The number of bins to use in dividing one-dimensional
#' marker distributions into discrete segments for the EMD calculation. Defaults to 100.
#'
#' @return A tibble.
#'
#' If format == "wide", the tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` and one column for each grouping variable,
#' one column for each extracted feature (the EMD between the distribution of a
#' given marker in a given cluster in the basal condition and the distribution of
#' that marker in a given cluster in a stimulated condition).
#' The names of each column containing cluster features is obtained using the following pattern:
#' "\{stimulation_id\}_\{marker_id\}@\{cluster_id\}_emd".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr map_lgl
#' @importFrom purrr map2_dbl
#' @importFrom tidyselect all_of
#' @importFrom stringr str_c
#' @importFrom rlang arg_match
#' @importFrom rlang enquo
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract emd of each cluster in each patient (using the "basal" stim
#' # condition as a reference) in wide format
#' tof_extract_emd(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     emd_col = stim,
#'     reference_level = "basal"
#' )
#'
#' # extract emd of each cluster (using the "basal" stim
#' # condition as a reference) in long format
#' tof_extract_emd(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     emd_col = stim,
#'     reference_level = "basal",
#'     format = "long"
#' )
#'
tof_extract_emd <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    emd_col,
    reference_level,
    format = c("wide", "long"),
    num_bins = 100
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    # make sure that stimulation information is provided
    if (missing(emd_col)) {
      stop("`emd_col` must be provided.")
    } else if (missing(reference_level)) {
      stop("`reference_level` must be provided.")
    }

    # check that the emd_col is not one of the group_cols
    emd_colname <-
      tidyselect::eval_select(
        expr = rlang::enquo(emd_col),
        data = tof_tibble
      ) |>
      names()

    group_colnames <-
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) |>
      names()

    if (emd_colname %in% group_colnames) {
      stop("`emd_col` should not be one of the `group_cols`")
    }

    # extract the stimulation "levels" present in the original tof_tibble
    stim_levels <-
      tof_tibble |>
      dplyr::pull({{emd_col}}) |>
      unique()

    non_basal_stim_levels <-
      setdiff(stim_levels, reference_level)

    # nest data
    nested_data <-
      tof_tibble |>
      dplyr::select({{group_cols}}, {{cluster_col}}, {{emd_col}}, {{marker_cols}}) |>
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "marker",
        values_to = "expression"
      ) |>
      tidyr::nest(data = "expression") |>
      tidyr::pivot_wider(names_from = {{emd_col}}, values_from = data) |>
      # filter out any rows that don't have a basal stimulation condition,
      # as this means that emd to basal is undefined
      dplyr::filter(
        purrr::map_lgl(.x = .data[[reference_level]], .f = ~!is.null(.x))
      ) |>
      # convert each column corresponding to a stimulation condition into a vector
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(stim_levels),
          ~ purrr::map(
            .x = .x,
            .f = ~ pull_unless_null(.x, expression)
          )
        )
      )

    emd_tibble <-
      nested_data |>
      dplyr::transmute(
        dplyr::across({{group_cols}}),
        {{cluster_col}},
        .data$marker,
        dplyr::across(
          tidyselect::all_of(non_basal_stim_levels),
          .fns = ~ purrr::map2_dbl(
            .x = .x,
            .y = .data[[reference_level]],
            .f = ~ tof_find_emd(.y, .x, num_bins = num_bins)
          )
        )
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(non_basal_stim_levels),
        names_to = "stimulation",
        values_to = "emd"
      )

    if (format == "wide") {
      emd_tibble <-
        emd_tibble |>
        dplyr::mutate(
          col_name =
            stringr::str_c(
              .data$stimulation,
              "_",
              .data$marker,
              "@",
              {{cluster_col}},
              "_emd"
            )
        ) |>
        dplyr::select({{group_cols}}, "emd", "col_name") |>
        tidyr::pivot_wider(
          names_from = "col_name",
          values_from = "emd"
        )
    }

    return(emd_tibble)
  }



# tof_extract_jsd --------------------------------------------------------------

#' Extract aggregated features from CyTOF data using the Jensen-Shannon Distance (JSD)
#'
#' This feature extraction function calculates the Jensen-Shannon Distance (JSD) between
#' the stimulated and unstimulated ("basal") experimental conditions of samples in a
#' CyTOF experiment. This calculation is performed across a
#' user-specified selection of CyTOF antigens and can be performed either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the feature extraction
#' calculation. Defaults to all numeric (integer or double) columns.
#' Supports tidyselect helpers.
#'
#' @param jsd_col An unquoted column name that indicates which
#' column in `tof_tibble` contains information about which stimulation condition each cell
#' was exposed to during data acquisition.
#'
#' If provided, the feature extraction will be
#' further broken down into subgroups by stimulation condition (and features from each stimulation
#' condition will be included as their own features in wide format).
#'
#' @param reference_level A string indicating what the value in `jsd_col`
#' corresponds to the basal stimulation condition (i.e. "basal" or "unstimulated").
#'
#' @param format A string indicating if the data should be returned in "wide" format
#' (the default; each cluster feature is given its own column) or in "long" format
#' (each cluster feature is provided as its own row).
#'
#' @param num_bins Optional. The number of bins to use in dividing one-dimensional
#' marker distributions into discrete segments for the JSD calculation. Defaults to 100.
#'
#' @return A tibble.
#'
#' If format == "wide", the tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` and one column for each grouping variable,
#' one column for each extracted feature (the JSD between the distribution of a
#' given marker in a given cluster in the basal condition and the distribution of
#' that marker in the same cluster in a stimulated condition).
#' The names of each column containing cluster features is obtained using the following pattern:
#' "\{stimulation_id\}_\{marker_id\}@\{cluster_id\}_jsd".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @importFrom purrr map_lgl
#' @importFrom purrr map2_dbl
#' @importFrom tidyselect all_of
#' @importFrom stringr str_c
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract jsd of each cluster in each patient (using the "basal" stim
#' # condition as a reference) in wide format
#' tof_extract_jsd(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     jsd_col = stim,
#'     reference_level = "basal"
#' )
#'
#' # extract jsd of each cluster (using the "basal" stim
#' # condition as a reference) in long format
#' tof_extract_jsd(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     jsd_col = stim,
#'     reference_level = "basal",
#'     format = "long"
#' )
#'
tof_extract_jsd <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    jsd_col,
    reference_level,
    format = c("wide", "long"),
    num_bins = 100
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    # make sure that stimulation information is provided
    if (missing(jsd_col)) {
      stop("`jsd_col` must be provided.")
    } else if (missing(reference_level)) {
      stop("`reference_level` must be provided.")
    }

    # check that the jsd_col is not one of the group_cols
    jsd_colname <-
      tidyselect::eval_select(
        expr = rlang::enquo(jsd_col),
        data = tof_tibble
      ) |>
      names()

    group_colnames <-
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) |>
      names()

    if (jsd_colname %in% group_colnames) {
      stop("`jsd_col` should not be one of the `group_cols`")
    }

    # extract the stimulation "levels" present in the original tof_tibble
    stim_levels <-
      tof_tibble |>
      dplyr::pull({{jsd_col}}) |>
      unique()

    non_basal_stim_levels <-
      setdiff(stim_levels, reference_level)

    # nest data
    nested_data <-
      tof_tibble |>
      dplyr::select({{group_cols}}, {{cluster_col}}, {{jsd_col}}, {{marker_cols}}) |>
      tidyr::pivot_longer(cols = {{marker_cols}}, names_to = "marker", values_to = "expression") |>
      tidyr::nest(data = expression) |>
      tidyr::pivot_wider(names_from = {{jsd_col}}, values_from = data) |>
      # filter out any rows that don't have a basal stimulation condition,
      # as this means that emd to basal is undefined
      dplyr::filter(
        purrr::map_lgl(.x = .data[[reference_level]], .f = ~!is.null(.x))
      ) |>
      # convert each column corresponding to a stimulation condition into a vector
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(stim_levels),
          ~ purrr::map(
            .x = .x,
            .f = ~ pull_unless_null(.x, expression)
          )
        )
      )

    jsd_tibble <-
      nested_data |>
      dplyr::transmute(
        dplyr::across({{group_cols}}),
        {{cluster_col}},
        .data$marker,
        dplyr::across(
          tidyselect::all_of(non_basal_stim_levels),
          .fns = ~ purrr::map2_dbl(
            .x = .x,
            .y = .data[[reference_level]],
            .f = ~ tof_find_jsd(.y, .x, num_bins = num_bins)
          )
        )
      ) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(non_basal_stim_levels),
        names_to = "stimulation",
        values_to = "jsd"
      )

    if (format == "wide") {
      jsd_tibble <-
        jsd_tibble |>
        mutate(
          col_name =
            stringr::str_c(
              .data$stimulation,
              "_",
              .data$marker,
              "@",
              {{cluster_col}},
              "_jsd"
            )
        ) |>
        dplyr::select({{group_cols}}, "jsd", "col_name") |>
        tidyr::pivot_wider(
          names_from = "col_name",
          values_from = "jsd"
        )
    }

    return(jsd_tibble)
  }



# tof_extract_features ---------------------------------------------------------


#' Extract aggregated, sample-level features from CyTOF data.
#'
#' This function wraps other members of the `tof_extract_*` function family to extract
#' sample-level features from both lineage (i.e. cell surface antigen) CyTOF channels
#' assumed to be stable across stimulation conditions and signaling CyTOF channels
#' assumed to change across stimulation conditions. Features are extracted for
#' each cluster within each independent sample (as defined with the `group_cols` argument).
#'
#' Lineage channels are specified using the `lineage_cols` argument, and their
#' extracted features will be measurements of central tendency (as computed by the
#' user-supplied `central_tendency_function`).
#'
#' Signaling channels are specified
#' using the `signaling_cols` argument, and their extracted features will depend
#' on the user's chosen `signaling_method`. If `signaling method` == "threshold"
#' (the default), \code{\link{tof_extract_threshold}} will be used to calculate the proportion of
#' cells in each cluster with signaling marker expression over `threshold` in each
#' stimulation condition. If `signaling_method` == "emd" or `signaling_method` == "jsd",
#' \code{\link{tof_extract_emd}} or \code{\link{tof_extract_jsd}} will be used to calculate the earth-mover's
#' distance (EMD) or Jensen-Shannon Distance (JSD), respectively, between the basal
#' condition and each of the stimulated conditions in each cluster for each sample.
#' Finally, if none of these options are chosen, \code{\link{tof_extract_central_tendency}}
#' will be used to calculate measurements of central tendency.
#'
#' In addition, \code{\link{tof_extract_proportion}} will be used to extract
#' the proportion of cells in each cluster will be computed for each sample.
#'
#' These calculations can be performed either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param group_cols Unquoted column names representing which columns in `tof_tibble`
#' should be used to break the rows of `tof_tibble` into subgroups for the feature
#' extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).
#'
#' @param stimulation_col Optional. An unquoted column name that indicates which
#' column in `tof_tibble` contains information about which stimulation condition each cell
#' was exposed to during data acquisition. If provided, the feature extraction will be
#' further broken down into subgroups by stimulation condition (and features from each stimulation
#' condition will be included as their own features in wide format).
#'
#' @param lineage_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be considered lineage markers in the
#' feature extraction calculation. Supports tidyselect helpers.
#'
#' @param signaling_cols  Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be considered signaling markers in the
#' feature extraction calculation. Supports tidyselect helpers.
#'
#' @param central_tendency_function The function that will be used to calculate
#' the measurement of central tendency for each cluster (to be used
#' as the dependent variable in the linear model). Defaults to \code{\link[stats]{median}}.
#'
#' @param signaling_method A string indicating which feature extraction method to use
#' for signaling markers (as identified by the `signaling_cols` argument). Options are
#' "threshold" (the default), "emd", "jsd", and "central tendency".
#'
#' @param basal_level A string indicating what the value in `stimulation_col`
#' corresponds to the basal stimulation condition (i.e. "basal" or "unstimulated").
#'
#' @param ... Optional additional arguments to be passed to tof_extract_threshold,
#' \code{\link{tof_extract_emd}}, or \code{\link{tof_extract_jsd}}.
#'
#' @return A tibble.
#'
#' The output tibble will have 1 row for each combination of
#' the grouping variables provided in `group_cols` (thus, each row will represent
#' what is considered a single "sample" based on the grouping provided).
#' It will have one column for each grouping variable and one column for each
#' extracted feature ("wide" format).
#'
#' @family feature extraction functions
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000),
#'         cluster_id = sample(letters, size = 1000, replace = TRUE),
#'         patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
#'         stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
#'     )
#'
#' # extract the following features from each cluster in each
#' # patient/stimulation:
#' #    - proportion of each cluster
#' #    - central tendency (median) of cd45 and cd38 in each cluster
#' #    - the proportion of cells in each cluster with cd34 expression over
#' #      the default threshold (asinh(10 / 5))
#' tof_extract_features(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     lineage_cols = c(cd45, cd38),
#'     signaling_cols = cd34,
#'     stimulation_col = stim
#' )
#'
#' # extract the following features from each cluster in each
#' # patient/stimulation:
#' #    - proportion of each cluster
#' #    - central tendency (mean) of cd45 and cd38 in each cluster
#' #    - the earth mover's distance between each cluster's cd34 histogram in
#' #      the "basal" and "stim" conditions
#' tof_extract_features(
#'     tof_tibble = sim_data,
#'     cluster_col = cluster_id,
#'     group_cols = patient,
#'     lineage_cols = c(cd45, cd38),
#'     signaling_cols = cd34,
#'     central_tendency_function = mean,
#'     stimulation_col = stim,
#'     signaling_method = "emd",
#'     basal_level = "basal"
#' )
#'
tof_extract_features <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    stimulation_col = NULL,
    lineage_cols,
    signaling_cols,
    central_tendency_function = stats::median,
    signaling_method = c("threshold", "emd", "jsd", "central tendency"),
    basal_level = NULL,
    ... # Optional additional arguments to be passed to tof_extract_threshold, tof_extract_emd, or tof_extract_jsd
  ) {

    # check that lineage and signaling cols are provided
    if (missing(lineage_cols) | missing(signaling_cols)) {
      stop("Both lineage_cols and signaling_cols must be provided.\n")
    }

    # checking signaling method argument
    signaling_method <- rlang::arg_match(signaling_method)

    # extract stimulation column's name
    # will be an empty character vector if stimulation_col = NULL
    stimulation_colname <-
      tidyselect::eval_select(
        expr = rlang::enquo(stimulation_col),
        data = tof_tibble
      ) |>
      names()

    # extract grouping columns' names'
    # will be an empty character vector if stimulation_col = NULL
    group_colnames <-
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) |>
      names()

    # check stimulation_col argument
    if (signaling_method %in% c("emd", "jsd") & length(stimulation_colname) == 0) {
      stop("stimulation_col must be specified for the chosen signaling_method.\n")
    }

    # check basal_level argument
    if (signaling_method %in% c("emd", "jsd") & is.null(basal_level)) {
      stop("basal_level must be specified for the chosen signaling_method.
           It should be a character vector (a quoted string) indicating which
           value in stimulation_col represents the basal (unstimulated) state.\n")
    }

    # find cluster abundance features
    abundance_features <-
      tof_tibble |>
      tof_extract_proportion(
        cluster_col = {{cluster_col}},
        group_cols = {{group_cols}},
        format = "wide"
      )

    # find lineage features
    lineage_features <-
      tof_tibble |>
      tof_extract_central_tendency(
        cluster_col = {{cluster_col}},
        group_cols = {{group_cols}},
        marker_cols = {{lineage_cols}},
        central_tendency_function = central_tendency_function,
        format = "wide"
      )

    # if no stimulations,
    if (length(stimulation_colname) == 0) {
      if (signaling_method == "central tendency") {
        # find signaling features as above
        signaling_features <-
          tof_tibble |>
          tof_extract_central_tendency(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            central_tendency_function = central_tendency_function,
            format = "wide",
            ...
          )

      } else if (signaling_method == "threshold") {
        # find signaling features using the threshold method
        signaling_features <-
          tof_tibble |>
          tof_extract_threshold(
            cluster_col = {{cluster_col}},
            group_cols = c({{group_cols}}),
            marker_cols = {{signaling_cols}},
            stimulation_col = {{stimulation_col}},
            format = "wide",
            ...
          )
      } else {
        stop("stimulation_col must be specified for the chosen method.\n")
      }

    # if there *are* stimulations
    } else {
      if (signaling_method == "central tendency") {
        # find signaling features as above
        signaling_features <-
          tof_tibble |>
          tof_extract_central_tendency(
            cluster_col = {{cluster_col}},
            group_cols = c({{group_cols}}),
            marker_cols = {{signaling_cols}},
            central_tendency_function = central_tendency_function,
            format = "wide",
            ...
          )
        # change tof_extract_central_tendency function
        # so that it can accommodate a stimulation column

        # find threshold features
      } else if (signaling_method == "threshold") {
        signaling_features <-
          tof_tibble |>
          tof_extract_threshold(
            cluster_col = {{cluster_col}},
            group_cols = c({{group_cols}}),
            marker_cols = {{signaling_cols}},
            stimulation_col = {{stimulation_col}},
            format = "wide",
            ...
          )

        # find emd features
      } else if (signaling_method == "emd") {
        signaling_features <-
          tof_tibble |>
          tof_extract_emd(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            emd_col = {{stimulation_col}},
            reference_level = basal_level,
            format = "wide",
            ...
          )

      } else if (signaling_method == "jsd") {
        signaling_features <-
          tof_tibble |>
          tof_extract_jsd(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            jsd_col = {{stimulation_col}},
            reference_level = basal_level,
            format = "wide",
            ...
          )
      }

      if (signaling_method %in% c("emd", "jsd")) {
        # find basal stimulation features
        basal_stim_features <-
          tof_tibble |>
          dplyr::filter({{stimulation_col}} == basal_level) |>
          tof_extract_central_tendency(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            central_tendency_function = central_tendency_function,
            format = "wide"
          ) |>
          dplyr::rename_with(
            .fn = function(x) stringr::str_c(basal_level, "_", x),
            .cols = -{{group_cols}}
          )

        if (length(group_colnames) == 0) {
          signaling_features <-
            basal_stim_features |>
            dplyr::bind_cols(signaling_features)

        } else {
          signaling_features <-
            basal_stim_features |>
            dplyr::left_join(signaling_features, by = group_colnames)
        }
      }
    }

    # if there are no grouping variables
    if (length(group_colnames) == 0) {
      final_features <-
        abundance_features |>
        dplyr::bind_cols(lineage_features) |>
        dplyr::bind_cols(signaling_features)

    } else {
      # if there are grouping columns
      final_features <-
        abundance_features |>
        dplyr::left_join(lineage_features, by = group_colnames) |>
        dplyr::left_join(signaling_features, by = group_colnames)
    }

    return(final_features)

  }





