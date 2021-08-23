# feature_extraction.R
# This file contains functions relevant to extracting patient- or sample-level
# features by aggregating single-cell data in tof_tibble objects.


#' Extract the proportion of cells in each cluster in a `tof_tibble`.
#'
#' This feature extraction function allows you to calculate the proportion of
#' cells in each cluster in a `tof_tibble` - either overall or when broken down
#' into subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
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
#' "prop@{cluster_id}".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols` and each cluster id (i.e. level) in `cluster_col`. It will have one column for
#' each grouping variable, one column for the cluster ids, and one column (`prop`) containing the
#' cluster proportions.
#'
#' @export
#'
#' @examples
#' NULL
tof_extract_proportion <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    format = c("wide", "long")
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    abundances <-
      tof_tibble %>%
      dplyr::group_by(dplyr::across({{group_cols}}), {{cluster_col}}) %>%
      dplyr::summarize(abundance = dplyr::n()) %>%
      # two lines below can be compressed into a transmute if multidplyr is not being used
      dplyr::mutate(
        prop = abundance / sum(abundance)
        #"{{cluster_col}}" := str_c("prop", {{cluster_col}}, sep = "@")
      ) %>%
      dplyr::select({{group_cols}}, {{cluster_col}}, prop)

    if (format == "wide") {
      abundances <-
        abundances %>%
        tidyr::pivot_wider(
          names_from = {{cluster_col}},
          values_from = prop,
          names_prefix = "prop@",
          # note that if a cluster is not present for a given group, it will
          # be filled in as having a relative abundance of 0.
          values_fill = 0
        )
    }

    return(dplyr::ungroup(abundances))
  }


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
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
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
#' "{marker_id}@{cluster_id}".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_extract_central_tendency <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    format = c("wide", "long")
  ) {

    # check format argument
    format <- rlang::arg_match(format)

    central_tendencies <-
      tof_tibble %>%
      # if cluster_col is not a character, make it one
      dplyr::mutate("{{cluster_col}}" := as.character({{cluster_col}})) %>%
      dplyr::group_by({{cluster_col}}, dplyr::across({{group_cols}})) %>%
      dplyr::summarize(dplyr::across({{marker_cols}}, central_tendency_function)) %>%
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "channel",
        values_to = "value"
      )

    if (format == "wide") {
      central_tendencies <-
        central_tendencies %>%
        tidyr::pivot_wider(
          names_from = c(channel, {{cluster_col}}),
          values_from = value,
          names_sep = "@"
        )
    }

    return(dplyr::ungroup(central_tendencies))
  }


#' Extract aggregated features from CyTOF data using a binary threshold
#'
#' This feature extraction function calculates the proportion of cells in a given cluster
#' that have a CyTOF antigen expression over a user-specified threshold across a
#' user-specified selection of CyTOF markers. These calculations can be done either
#' overall (across all cells in the dataset) or after breaking down the cells into
#' subgroups using `group_cols`.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
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
#' one column for each extracted feature (the central tendency of a given marker in a given cluster).
#' The names of each column containing cluster features is obtained using the following pattern:
#' "{marker_id}@{cluster_id}".
#'
#' If format == "long", the tibble will have 1 row for each combination of the grouping variables
#' in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
#' It will have one column for each grouping variable, one column for the cluster ids, one
#' column for the CyTOF channel names, and one column (`value`) containing the features.
#'
#' @export
#'
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
      tof_tibble %>%
      # if the cluster column isn't a character vector, make it one
      dplyr::mutate("{{cluster_col}}" := as.character({{cluster_col}})) %>%
      dplyr::group_by(dplyr::across({{group_cols}}), {{cluster_col}}, {{stimulation_col}}) %>%
      dplyr::summarize(dplyr::across({{marker_cols}}, ~ mean(.x > threshold))) %>%
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "channel",
        values_to = "values"
      )

    if (format == "wide") {
      threshold_features <-
        threshold_features %>%
        dplyr::group_by(dplyr::across({{group_cols}})) %>%
        dplyr::transmute(
          col_names =
            stringr::str_c({{stimulation_col}}, "_", channel, "@", {{cluster_col}}) %>%
            stringr::str_remove("^_"),
          values
        ) %>%
        tidyr::pivot_wider(
          names_from = col_names,
          values_from = values
        )
    }

    return(dplyr::ungroup(threshold_features))
  }





#' Title
#'
#' @param tof_tibble TO DO
#'
#' @param cluster_col TO DO
#'
#' @param group_cols TO DO
#'
#' @param marker_cols TO DO
#'
#' @param stimulation_col TO DO
#'
#' @param basal_level TO DO
#'
#' @param format TO DO
#'
#' @param num_bins TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_extract_emd <-

  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    stimulation_col,
    basal_level, # a string indicating what the basal condition is called in stimulation_col,
    format = c("wide", "long"),
    num_bins = 100
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    # extract the stimulation "levels" present in the original tof_tibble
    stim_levels <-
      tof_tibble %>%
      dplyr::pull({{stimulation_col}}) %>%
      base::unique()

    non_basal_stim_levels <-
      generics::setdiff(stim_levels, basal_level)

    # nest data
    nested_data <-
      tof_tibble %>%
      dplyr::select({{group_cols}}, {{cluster_col}}, {{stimulation_col}}, {{marker_cols}}) %>%
      tidyr::pivot_longer(cols = {{marker_cols}}, names_to = "marker", values_to = "expression") %>%
      tidyr::nest(data = expression) %>%
      tidyr::pivot_wider(names_from = {{stimulation_col}}, values_from = data) %>%
      # filter out any rows that don't have a basal stimulation condition,
      # as this means that emd to basal is undefined
      dplyr::filter(
        purrr::map_lgl(.x = .data[[basal_level]], .f = ~!is.null(.x))
      ) %>%
      # convert each column corresponding to a stimulation condition into a vector
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(stim_levels),
          ~ purrr::map(
            .x = .x,
            .f = ~pull_unless_null(.x, expression)
          )
        )
      )

    emd_tibble <-
      nested_data %>%
      dplyr::transmute(
        {{group_cols}},
        {{cluster_col}},
        marker,
        dplyr::across(
          tidyselect::all_of(non_basal_stim_levels),
          .fns = ~purrr::map2_dbl(
            .x = .x,
            .y = .data[[basal_level]],
            .f = ~tof_find_emd(.y, .x, num_bins = num_bins)
          )
        )
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(non_basal_stim_levels),
        names_to = "stimulation",
        values_to = "emd"
      )

    if (format == "wide") {
      emd_tibble <-
        emd_tibble %>%
        mutate(
          col_name =
            stringr::str_c(stimulation, "_", marker, "@", {{cluster_col}}, "_emd")
        ) %>%
        dplyr::select({{group_cols}}, emd, col_name) %>%
        tidyr::pivot_wider(
          names_from = col_name,
          values_from = emd
        )
    }

    return(emd_tibble)
  }


#' Title
#'
#' @param tof_tibble TO DO
#'
#' @param cluster_col TO DO
#'
#' @param group_cols TO DO
#'
#' @param marker_cols TO DO
#'
#' @param stimulation_col TO DO
#'
#' @param basal_level TO DO
#'
#' @param format TO DO
#'
#' @param num_bins TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_extract_jsd <-

  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    marker_cols = where(tof_is_numeric),
    stimulation_col,
    basal_level, # a string indicating what the basal condition is called in stimulation_col,
    format = c("wide", "long"),
    num_bins = 100
  ) {
    # check format argument
    format <- rlang::arg_match(format)

    # extract the stimulation "levels" present in the original tof_tibble
    stim_levels <-
      tof_tibble %>%
      dplyr::pull({{stimulation_col}}) %>%
      base::unique()

    non_basal_stim_levels <-
      generics::setdiff(stim_levels, basal_level)

    # nest data
    nested_data <-
      tof_tibble %>%
      dplyr::select({{group_cols}}, {{cluster_col}}, {{stimulation_col}}, {{marker_cols}}) %>%
      tidyr::pivot_longer(cols = {{marker_cols}}, names_to = "marker", values_to = "expression") %>%
      tidyr::nest(data = expression) %>%
      tidyr::pivot_wider(names_from = {{stimulation_col}}, values_from = data) %>%
      # filter out any rows that don't have a basal stimulation condition,
      # as this means that emd to basal is undefined
      dplyr::filter(
        purrr::map_lgl(.x = .data[[basal_level]], .f = ~!is.null(.x))
      ) %>%
      # convert each column corresponding to a stimulation condition into a vector
      dplyr::mutate(
        dplyr::across(
          tidyselect::all_of(stim_levels),
          ~ purrr::map(
            .x = .x,
            .f = ~pull_unless_null(.x, expression)
          )
        )
      )

    jsd_tibble <-
      nested_data %>%
      dplyr::transmute(
        {{group_cols}},
        {{cluster_col}},
        marker,
        dplyr::across(
          tidyselect::all_of(non_basal_stim_levels),
          .fns = ~purrr::map2_dbl(
            .x = .x,
            .y = .data[[basal_level]],
            .f = ~tof_find_jsd(.y, .x, num_bins = num_bins)
          )
        )
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(non_basal_stim_levels),
        names_to = "stimulation",
        values_to = "jsd"
      )

    if (format == "wide") {
      jsd_tibble <-
        jsd_tibble %>%
        mutate(
          col_name =
            stringr::str_c(stimulation, "_", marker, "@", {{cluster_col}}, "_jsd")
        ) %>%
        dplyr::select({{group_cols}}, jsd, col_name) %>%
        tidyr::pivot_wider(
          names_from = col_name,
          values_from = jsd
        )
    }

    return(jsd_tibble)
  }



#' Title
#'
#' @param tof_tibble TO DO
#'
#' @param cluster_col TO DO
#'
#' @param group_cols TO DO
#'
#' @param stimulation_col TO DO
#'
#' @param lineage_cols TO DO
#'
#' @param signaling_cols TO DO
#'
#' @param central_tendency_function TO DO
#'
#' @param signaling_method TO DO
#'
#' @param basal_level TO DO
#'
#' @param ... TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_extract_features <-
  function(
    tof_tibble,
    cluster_col,
    group_cols = NULL,
    stimulation_col = NULL,
    lineage_cols,
    signaling_cols,
    central_tendency_function = mean,
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
      rlang::enquo(stimulation_col) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
      names()

    # extract grouping columns' names'
    # will be an empty character vector if stimulation_col = NULL
    group_colnames <-
      rlang::enquo(group_cols) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
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
      tof_tibble %>%
      tof_extract_proportion(
        cluster_col = {{cluster_col}},
        group_cols = {{group_cols}},
        format = "wide"
      )

    # find lineage features
    lineage_features <-
      tof_tibble %>%
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
          tof_tibble %>%
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
          tof_tibble %>%
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
          tof_tibble %>%
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
          tof_tibble %>%
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
          tof_tibble %>%
          tof_extract_emd(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            stimulation_col = {{stimulation_col}},
            basal_level = basal_level,
            format = "wide",
            ...
          )

      } else if (signaling_method == "jsd") {
        signaling_features <-
          tof_tibble %>%
          tof_extract_jsd(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            stimulation_col = {{stimulation_col}},
            basal_level = basal_level,
            format = "wide",
            ...
          )
      }

      if (signaling_method %in% c("emd", "jsd")) {
        # find basal stimulation features
        basal_stim_features <-
          tof_tibble %>%
          dplyr::filter({{stimulation_col}} == basal_level) %>%
          tof_extract_central_tendency(
            cluster_col = {{cluster_col}},
            group_cols = {{group_cols}},
            marker_cols = {{signaling_cols}},
            central_tendency_function = central_tendency_function,
            format = "wide"
          ) %>%
          dplyr::rename_with(
            .fn = function(x) stringr::str_c(basal_level, "_", x),
            .cols = -{{group_cols}}
          )

        if (length(group_colnames) == 0) {
          signaling_features <-
            basal_stim_features %>%
            dplyr::bind_cols(signaling_features)

        } else {
          signaling_features <-
            basal_stim_features %>%
            dplyr::left_join(signaling_features, by = group_colnames)
        }
      }

      # if there are no grouping variables
      if (length(group_colnames) == 0) {
        final_features <-
          abundance_features %>%
          dplyr::bind_cols(lineage_features) %>%
          dplyr::bind_cols(signaling_features)

      } else {
        # if there are grouping columns
        final_features <-
          abundance_features %>%
          dplyr::left_join(lineage_features, by = group_colnames) %>%
          dplyr::left_join(signaling_features, by = group_colnames)
      }

      return(final_features)
    }
  }





