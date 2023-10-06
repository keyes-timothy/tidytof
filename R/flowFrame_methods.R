# flowFrame_methods.R
# This file contains tidyverse-style methods for flowCore::flowFrame and
# flowCore::flowSet objects.


# dplyr methods ----------------------------------------------------------------


## mutate-------------

#' Create, modify, and delete columns.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Name-value pairs. The name (the left side of the equals sign)
#' gives the name of the column in the output. The right side of the equation
#' performs computations using the names of each channel according to
#'  \code{\link[flowCore]{featureNames}}. Supports tidyselection.
#'
#' @returns A \code{\link[flowCore]{flowFrame}}. The output has the following properties:
#' * Columns from .data will be preserved according to the .keep argument.
#' * Existing columns that are modified by ... will always be returned in their original location.
#' * New columns created through ... will be placed according to the .before and .after arguments.
#' * The number of rows is not affected.
#' * Columns given the value NULL will be removed.
#'
#' @importFrom dplyr mutate
#'
#' @importFrom flowCore identifier
#'
#' @export
#'
mutate.flowFrame <- function(.data, ...) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::mutate(...) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Create, modify, and delete columns.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Name-value pairs. The name (the left side of the equals sign)
#' gives the name of the column in the output. The right side of the equation
#' performs computations using the names of each channel according to
#'  \code{\link[flowCore]{featureNames}}. Supports tidyselection.
#'
#'  @returns A \code{\link[flowCore]{flowSet}}. The output has the following properties:
#' * Columns from .data will be preserved according to the .keep argument.
#' * Existing columns that are modified by ... will always be returned in their original location.
#' * New columns created through ... will be placed according to the .before and .after arguments.
#' * The number of rows is not affected.
#' * Columns given the value NULL will be removed.
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
mutate.flowSet <- function(.data, ...) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::mutate(.data[[i]], ...)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## select -----------

#' Keep or drop columns using their names and types.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... One or more unquoted expressions separated by commas. Variables names
#' (as specified by \code{\link[flowCore]{featureNames}}) can be used as if they
#' were positions in the \code{\link[flowCore]{flowFrame}}). Supports tidyselection.
#'
#' @returns A \code{\link[flowCore]{flowFrame}}. The output has the following properties:
#' * Rows are not affected.
#' * Output columns are a subset of input columns, potentially with a different order. Columns will be renamed if new_name = old_name form is used.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @export
select.flowFrame <- function(.data, ...) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::select(...) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Keep or drop columns using their names and types.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... One or more unquoted expressions separated by commas. Variables names
#' (as specified by the \code{\link[flowCore]{featureNames}} of the component flowFrames
#' that make up the flowSet) can be used as if they
#' were positions in the \code{\link[flowCore]{flowSet}}). Supports tidyselection.
#'
#' @returns A \code{\link[flowCore]{flowSet}}. The output has the following properties:
#' * Rows are not affected.
#' * Output columns are a subset of input columns, potentially with a different order. Columns will be renamed if new_name = old_name form is used.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @export
select.flowSet <- function(.data, ...) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::select(.data[[i]], ...)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## rename -----------

#' Rename columns in a \code{\link[flowCore]{flowFrame}}
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Unquoted name-value pairs (as specified by \code{\link[flowCore]{featureNames}}).
#' Use new_name = old_name to rename selected columns
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Rows are not affected.
#' * Column names are changed; column order is preserved.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @importFrom dplyr rename
#'
#' @importFrom flowCore identifier
#'
#' @export
#'
#'
rename.flowFrame <- function(.data, ...) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::rename(...) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Rename columns in a \code{\link[flowCore]{flowSet}}
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Unquoted name-value pairs (as specified by the \code{\link[flowCore]{featureNames}} of
#' the \code{\link[flowCore]{flowFrame}}s making up the \code{\link[flowCore]{flowSet}}).
#' Use new_name = old_name to rename selected columns
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Rows are not affected.
#' * Column names are changed; column order is preserved.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @importFrom dplyr rename
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
#'
rename.flowSet <- function(.data, ...) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::rename(.data[[i]], ...)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## rename_with -----------

#' Rename columns in a \code{\link[flowCore]{flowFrame}}
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param .fn A function used to transform the selected .cols.
#' Should return a character vector the same length as the input.
#'
#' @param .cols Unquoted column names indicating which columns to rename
#' (as specified by \code{\link[flowCore]{featureNames}}).
#'
#' @param ... Additional arguments passed onto .fn.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Rows are not affected.
#' * Column names are changed; column order is preserved.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @importFrom dplyr rename_with
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
rename_with.flowFrame <- function(.data, .fn, .cols = dplyr::everything(), ...) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::rename_with(.fn = .fn, .cols = {{ .cols }}, ...) |>
    as_flowFrame()
  flowCore::identifier(result) <- flowCore::identifier(.data)
  return(result)
}

#' Rename columns in a \code{\link[flowCore]{flowSet}}
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param .fn A function used to transform the selected .cols. Should return a character vector the same length as the input.
#'
#' @param .cols Unquoted column names indicating which columns to rename (as specified by the \code{\link[flowCore]{featureNames}} of
#' the \code{\link[flowCore]{flowFrame}}s making up the \code{\link[flowCore]{flowSet}}).
#'
#' @param ... Additional arguments passed onto .fn.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Rows are not affected.
#' * Column names are changed; column order is preserved.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @importFrom dplyr rename_with
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
rename_with.flowSet <- function(.data, .fn, .cols = dplyr::everything(), ...) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::rename_with(.data[[i]], .fn = .fn, .cols = {{ .cols }}, ...)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}


## filter ----------

#' Keep rows that match a condition.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Expressions that return a logical value, and are defined in terms
#' of the variables in the \code{\link[flowCore]{featureNames}} of .data.
#' If multiple expressions are included, they are combined with the & operator.
#' Only rows for which all conditions evaluate to TRUE are kept.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .preserve Unused.
#'
#' @returns An object of the same type as .data. The output has the
#' following properties:
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @importFrom dplyr filter
#'
#' @importFrom flowCore identifier
#'
#' @export
#'
#'
filter.flowFrame <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::filter(..., .by = {{ .by }}, .preserve = .preserve) |>
    as_flowFrame()
  flowCore::identifier(result) <- flowCore::identifier(.data)
  return(result)
}

#' Keep rows that match a condition.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Expressions that return a logical value, and are defined in terms
#' of the variables in the \code{\link[flowCore]{featureNames}} of the
#' \code{\link[flowCore]{flowFrame}}s in .data.
#' If multiple expressions are included, they are combined with the & operator.
#' Only rows for which all conditions evaluate to TRUE are kept.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .preserve Unused.
#'
#' @returns An object of the same type as .data. The output has the
#' following properties:
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @importFrom dplyr filter
#'
#' @importFrom flowCore pData
#' @importFrom flowCore identifier
#'
#' @export
#'
filter.flowSet <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::filter(.data[[i]], ..., .by = {{ .by }}, .preserve = .preserve)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## arrange -----------

#' Order rows using column values
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Variables, or functions of variables, to arrange by.
#'
#' @param .by_group Unused.
#'
#' @returns An object of the same type as .data. The output has the following
#' properties:
#' * All rows appear in the output, but (usually) in a different place.
#' * Columns are not modified.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @importFrom dplyr arrange
#'
#' @importFrom flowCore identifier
#'
#' @export
arrange.flowFrame <- function(.data, ..., .by_group = FALSE) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::arrange(..., .by_group = .by_group) |>
    as_flowFrame()
  flowCore::identifier(result) <- flowCore::identifier(.data)
  return(result)
}

#' Order rows using column values
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Variables, or functions of variables, to arrange by.
#'
#' @param .by_group Unused.
#'
#' @returns An object of the same type as .data. The output has the following
#' properties:
#' * All rows appear in the output, but (usually) in a different place.
#' * Columns are not modified.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @importFrom dplyr arrange
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
arrange.flowSet <- function(.data, ..., .by_group = FALSE) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::arrange(.data[[i]], ..., .by_group = .by_group)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## transmute -----------

#' Create, modify, and delete columns.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Name-value pairs. The name (the left side of the equals sign)
#' gives the name of the column in the output. The right side of the equation
#' performs computations using the names of each channel according to
#'  \code{\link[flowCore]{featureNames}}. Supports tidyselection.
#'
#' @returns A \code{\link[flowCore]{flowFrame}}. The output has the following properties:
#' * Columns created or modified through ... will be returned in the order specified by ....
#' * The number of rows is not affected.
#' * Columns given the value NULL will be removed.
#' * The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
#'
#' @importFrom dplyr transmute
#'
#' @importFrom flowCore identifier
#'
#' @export
#'
transmute.flowFrame <- function(.data, ...) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::transmute(...) |>
    as_flowFrame()
  flowCore::identifier(result) <- flowCore::identifier(.data)
  return(result)
}

#' Create, modify, and delete columns.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Name-value pairs. The name (the left side of the equals sign)
#' gives the name of the column in the output. The right side of the equation
#' performs computations using the names of each channel according to the
#'  \code{\link[flowCore]{featureNames}} of .data's constituent
#'  \code{\link[flowCore]{flowFrame}}s. Supports tidyselection.
#'
#' @returns A \code{\link[flowCore]{flowSet}}. The output has the following properties:
#' * Columns created or modified through ... will be returned in the order specified by ....
#' * The number of rows is not affected.
#' * Columns given the value NULL will be removed.
#' * The \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} will be preserved.
#'
#' @importFrom dplyr transmute
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
transmute.flowSet <- function(.data, ...) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::transmute(.data[[i]], ...)
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

## summarize -----------

#' Summarize a flowFrame.
#'
#' @param .data .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Name-value pairs of summary functions.
#' The name will be the name of the variable in the result.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .groups Grouping structure of the result.
#' * "drop_last": dropping the last level of grouping.
#' * "drop": All levels of grouping are dropped.
#' * "keep": Same grouping structure as .data.
#' * "rowwise": Each row is its own group.
#'
#' @returns A data.frame containing the summarized result.
#'
#' @importFrom dplyr summarize
#'
#' @export
#'
summarise.flowFrame <- function(.data, ..., .by = NULL, .groups = NULL) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::summarize(..., .by = {{ .by }}, .groups = NULL)
  return(result)
}

#' Summarize a flowFrame.
#'
#' @param .data .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Name-value pairs of summary functions.
#' The name will be the name of the variable in the result.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .groups Grouping structure of the result.
#' * "drop_last": dropping the last level of grouping.
#' * "drop": All levels of grouping are dropped.
#' * "keep": Same grouping structure as .data.
#' * "rowwise": Each row is its own group.
#'
#' @returns A data.frame containing the summarized result.
#'
#' @importFrom dplyr summarise
#'
#' @export
summarize.flowFrame <- function(.data, ..., .by = NULL, .groups = NULL) {
  result <-
    dplyr::summarise(.data = .data, ..., .by = {{ .by }}, .groups = .groups)
  return(result)
}

#' Summarize a flowSet.
#'
#' @param .data .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Name-value pairs of summary functions.
#' The name will be the name of the variable in the result.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .groups Grouping structure of the result.
#' * "drop_last": dropping the last level of grouping.
#' * "drop": All levels of grouping are dropped.
#' * "keep": Same grouping structure as .data.
#' * "rowwise": Each row is its own group.
#'
#' @returns A data.frame containing the summarized result.
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
summarise.flowSet <- function(.data, ..., .by = NULL, .groups = NULL) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::summarize(.data[[i]], ..., .by = {{ .by }}, .groups = .groups) |>
      dplyr::mutate(.flowframe_identifier = identifier)
  }
  result <-
    dplyr::bind_rows(result_list) |>
    dplyr::left_join(
      flowCore::pData(.data),
      by = c(".flowframe_identifier" = "name")
    ) |>
    dplyr::select(-dplyr::any_of(c(".tidytof_unique_identifier")))
  return(result)
}

#' Summarize a flowSet.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Name-value pairs of summary functions.
#' The name will be the name of the variable in the result.
#'
#' @param .by Optionally, a selection of columns to group by for just this
#' operation, functioning as an alternative to group_by().
#'
#' @param .groups Grouping structure of the result.
#' * "drop_last": dropping the last level of grouping.
#' * "drop": All levels of grouping are dropped.
#' * "keep": Same grouping structure as .data.
#' * "rowwise": Each row is its own group.
#'
#' @returns A data.frame containing the summarized result.
#'
#' @importFrom dplyr summarise
#'
#' @export
summarize.flowSet <- function(.data, ..., .by = NULL, .groups = NULL) {
  result <-
    dplyr::summarise(.data = .data, ..., .by = {{ .by }}, .groups = .groups)
  return(result)
}

## group_by -----------

#' Group a flowFrame into a flowSet using one or more variables.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Unquoted variables or columns to group by according to .data's
#'  \code{\link[flowCore]{featureNames}}.
#'
#' @param .add Unused.
#'
#' @param .drop Unused.
#'
#' @returns A \code{\link[flowCore]{flowSet}} containing one
#' \code{\link[flowFrame]{flowFrame}} for each of the unique combinations of columns
#' selected in .... Metadata about grouping columns will be stored in the output
#' \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}}.
#'
#' @importFrom dplyr select
#' @importFrom dplyr any_of
#' @importFrom dplyr group_by_drop_default
#'
#' @export
group_by.flowFrame <-
  function(.data, ..., .add = FALSE, .drop = dplyr::group_by_drop_default(.data)) {
    tof_tibble <-
      .data |>
      as_tof_tbl(.name_method = "featureNames")

    group_colnames <-
      tof_tibble |>
      dplyr::select(...) |>
      colnames()

    result <-
      tof_tibble |>
      as_flowSet(group_cols = dplyr::any_of(group_colnames))

    return(result)
  }

## ungroup ---------

# TODO: Throws an error when name is present

#' Ungroup a flowSet
#'
#' @param x A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Variables/columns in \code{\link[flowCore]{pData}} to remove
#' from the grouping. The column ".tidytof_unique_identifier" is used internally
#' and will not have any effect on the ungrouping.
#'
#' @returns A \code{\link[flowCore]{flowFrame}} or
#' \code{\link[flowCore]{flowSet}} depending on the degree of ungrouping.
#' Note that unnest-ing and ungrouping a \code{\link[flowCore]{flowSet}} are
#' equivalent.
#'
#' @importFrom dplyr across
#' @importFrom dplyr any_of
#' @importFrom dplyr mutate
#' @importFrom dplyr where
#' @importFrom dplyr select
#'
#' @importFrom flowCore pData
#'
#' @importFrom rlang enquos
#' @importFrom rlang !!!
#'
#' @export
#'
ungroup.flowSet <- function(x, ...) {
  # Note that .tidytof_unique_identifier is not meant to be directly accessed
  # by the user and will not have any effect if chosen for the ungrouping
  .cols <- rlang::enquos(...)

  metadata_colnames <-
    flowCore::pData(x) |>
    dplyr::select(-dplyr::any_of(".tidytof_unique_identifier")) |>
    colnames()

  ungroup_colnames <-
    flowCore::pData(x) |>
    dplyr::select(!!! .cols) |>
    colnames()

  # gives expected behavior such that calling ungroup without any arguments
  # will return a fully ungrouped dataset
  if(length(ungroup_colnames) == 0) {
    ungroup_colnames <- metadata_colnames
  }

  group_colnames <- setdiff(metadata_colnames, ungroup_colnames)

  # convert any character or factor columns to be ungrouped into numeric values
  # ignore the name column, which is special in flowCore and shouldn't be changed
  character_ungroup_colnames <-
    flowCore::pData(x) |>
    dplyr::select(dplyr::any_of(ungroup_colnames)) |>
    dplyr::select(dplyr::where(\(.x) !is.numeric(.x))) |>
    colnames()

  # non-numeric values cannot be "ungrouped" and turned into a column in a flowFrame
  # so we convert any non-numeric columns that are being ungrouped into numerics
  # character vectors are converted in alphabetical order; factors in level order
  flowCore::pData(x) <-
    # causes a bug when "name" is one of the metadata variables being ungrouped
    # because it causes a key mismatch in #tof_tbl.R's as_tof_tbl function
    flowCore::pData(x) |>
    # experimental
    dplyr::mutate(.tidytof_unique_identifier = name) |>
    dplyr::mutate(dplyr::across(dplyr::any_of(character_ungroup_colnames), \(.x) as.numeric(as.factor(.x))))

  if (length(group_colnames) != 0) {
    result <-
      x |>
      as_tof_tbl(.name_method = "featureNames", include_metadata = TRUE) |>
      as_flowSet(group_cols = dplyr::any_of(group_colnames))
  } else {
    result <-
      x |>
      as_tof_tbl(.name_method = "featureNames", include_metadata = TRUE) |>
      as_flowSet()
  }

  return(result)
}

## count -----------

#' Count the observations in each group.
#'
#' @param x A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Variables to group by, named according to
#' \code{\link[flowCore]{featureNames}}
#'
#' @param wt If NULL (the default), counts the number of rows in each group.
#' If a variable, computes sum(wt) for each group.
#'
#' @param sort If TRUE, will show the largest groups at the top.
#'
#' @param name If omitted, it will default to n. If there's already a column
#' called n, it will use nn. If there's a column called n and nn, it'll use nnn,
#' and so on, adding ns until it gets a new name.
#'
#' @returns A data.frame containing the groupwise counts.
#'
#' @importFrom dplyr count
#'
#' @export
#'
count.flowFrame <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  tof_tibble <-
    x |>
    as_tof_tbl(.name_method = "featureNames")

    result <-
      dplyr::count(x = tof_tibble, ..., wt = {{ wt }}, sort = sort, name = name)
  return(result)
}

#' Count the observations in each group.
#'
#' @param x A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Variables to group by, named according to
#' \code{\link[flowCore]{featureNames}}
#'
#' @param wt If NULL (the default), counts the number of rows in each group.
#' If a variable, computes sum(wt) for each group.
#'
#' @param sort If TRUE, will show the largest groups at the top.
#'
#' @param name If omitted, it will default to n. If there's already a column
#' called n, it will use nn. If there's a column called n and nn, it'll use nnn,
#' and so on, adding ns until it gets a new name.
#'
#' @returns A data.frame containing the groupwise counts.
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr count
#' @importFrom dplyr mutate
#'
#' @importFrom flowCore identifier
#'
#' @export
count.flowSet <- function(x, ..., wt = NULL, sort = FALSE, name = NULL) {
  # result_list <- list()
  # for (i in 1:length(x)) {
  #   identifier <- flowCore::identifier(x[[i]])
  #   result_list[[identifier]] <-
  #     dplyr::count(x[[i]], ..., wt = wt, sort = sort, name = name)
  #
  #   if (ncol(result_list[[identifier]]) == 1) {
  #     result_list[[identifier]] <-
  #       result_list[[identifier]] |>
  #       dplyr::mutate(.flowframe_identifier = identifier)
  #   }
  #
  # }
  # result <- dplyr::bind_rows(result_list)

  result <-
    x |>
    as_tof_tbl(.name_method = "featureNames") |>
    dplyr::count(..., wt = {{ wt }}, sort = sort, name = name)

  return(result)
}

## pull -----------

#' Extract a single column.
#'
#' pull() is similar to $. It's mostly useful because it looks a little nicer
#' in pipes.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}.
#'
#' @param var A variable specified as:
#' * a literal variable name
#' * a positive integer, giving the position counting from the left
#' * a negative integer, giving the position counting from the right.
#'
#' @param name An optional parameter that specifies the column to be used as
#' names for a named vector. Specified in a similar manner as var.
#'
#' @param ... For use by methods.
#'
#' @returns A vector the same size as .data.
#'
#' @importFrom dplyr pull
#'
#' @export
pull.flowFrame <- function(.data, var = -1, name = NULL, ...) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::pull(var = {{ var }}, name = {{ name }}, ...)
  return(result)
}

#' Extract a single column.
#'
#' pull() is similar to $. It's mostly useful because it looks a little nicer
#' in pipes.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}.
#'
#' @param var A variable specified as:
#' * a literal variable name
#' * a positive integer, giving the position counting from the left
#' * a negative integer, giving the position counting from the right.
#'
#' @param name An optional parameter that specifies the column to be used as
#' names for a named vector. Specified in a similar manner as var.
#'
#' @param ... For use by methods.
#'
#' @returns A vector the same size as .data.
#'
#' @importFrom dplyr pull
#'
#' @export
pull.flowSet <- function(.data, var = -1, name = NULL, ...) {
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::pull(var = {{ var }}, name = {{ name }}, ...)
  return(result)
}


## slice methods -----------

### slice---------

#' Subset rows using their positions
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Integer row values (to keep).
#'
#' @param .by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param .preserve Currently unused.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice
#'
#' @importFrom flowCore identifier
#'
#' @export
#'
slice.flowFrame <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::slice(..., .by = {{ .by }}, .preserve = .preserve) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}


#' Subset rows using their positions
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Integer row values (to keep).
#'
#' @param .by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param .preserve Currently unused.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
slice.flowSet <- function(.data, ..., .by = NULL, .preserve = FALSE) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::slice(
        .data[[i]],
        ...,
        .by = {{ .by }},
        .preserve = .preserve
      )
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}


### slice_sample ------------

#' Subset rows randomly
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by  Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param weight_by Sampling weights. This must evaluate to a vector of
#' non-negative numbers the same length as the input.
#' Weights are automatically standardized to sum to 1.
#'
#' @param replace Should sampling be performed with (TRUE) or without
#' (FALSE, the default) replacement.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} is preserved.
#'
#' @importFrom dplyr slice_sample
#'
#' @importFrom flowCore identifier
#'
#' @export
slice_sample.flowFrame <- function(.data, ..., n, prop, by = NULL, weight_by = NULL, replace = FALSE) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::slice_sample(
      ...,
      n = n,
      prop = prop,
      by = {{ by }},
      weight_by = weight_by,
      replace = replace
    ) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Subset rows randomly
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param weight_by Sampling weights. This must evaluate to a vector of
#' non-negative numbers the same length as the input.
#' Weights are automatically standardized to sum to 1.
#'
#' @param replace Should sampling be performed with (TRUE) or without
#' (FALSE, the default) replacement.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#'
#' @importFrom dplyr slice_sample
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
slice_sample.flowSet <- function(.data, ..., n, prop,  by = NULL, weight_by = NULL, replace = FALSE) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::slice_sample(
        .data[[i]],
        ...,
        n = n,
        prop = prop,
        by = {{ by }},
        weight_by = weight_by,
        replace = replace
      )
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}


### slice_head --------

#' Subset rows at the head of a data structure.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} is preserved.
#'
#' @importFrom dplyr slice_head
#'
#' @importFrom flowCore identifier
#'
#' @export
slice_head.flowFrame <- function(.data, ..., n, prop, by = NULL) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::slice_head(
      ...,
      n = n,
      prop = prop,
      by = {{ by }}
    ) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Subset rows at the head of a data structure.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice_head
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
slice_head.flowSet <- function(.data, ..., n, prop, by = NULL) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::slice_sample(
        .data[[i]],
        ...,
        n = n,
        prop = prop,
        by = {{ by }}
      )
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

### slice_tail --------

#' Subset rows at the tail of a data structure.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} is preserved.
#'
#' @importFrom dplyr slice_tail
#'
#' @importFrom flowCore identifier
#'
#' @export
slice_tail.flowFrame <- function(.data, ..., n, prop, by = NULL) {
  identifier <- flowCore::identifier(.data)
  tof_tibble <-
    .data |>
    as_tof_tbl(.name_method = "featureNames")
  result <-
    tof_tibble |>
    dplyr::slice_tail(
      ...,
      n = n,
      prop = prop,
      by = {{ by }}
    ) |>
    as_flowFrame()
  flowCore::identifier(result) <- identifier
  return(result)
}

#' Subset rows at the tail of a data structure.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice_tail
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
#'
slice_tail.flowSet <- function(.data, ..., n, prop, by = NULL) {
  result_list <- list()
  for (i in 1:length(.data)) {
    identifier <- flowCore::identifier(.data[[i]])
    result_list[[identifier]] <-
      dplyr::slice_tail(
        .data[[i]],
        ...,
        n = n,
        prop = prop,
        by = {{ by }}
      )
  }
  result <- as(result_list, "flowSet")
  flowCore::pData(result) <- flowCore::pData(.data)
  return(result)
}

### slice_max ------------

#' Subset rows of a data structure in order.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param order_by Variable or function of variables to order by.
#' To order by multiple variables, wrap them in a data frame or tibble.
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param with_ties Should ties be kept together? The default, TRUE, may return
#' more rows than you request. Use FALSE to ignore ties, and return the first
#' n rows.
#'
#' @param na_rm Should missing values in order_by be removed from the result?
#' If FALSE, NA values are sorted to the end so they will only be included if
#' there are insufficient non-missing values to reach n/prop.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} is preserved.
#'
#' @importFrom dplyr slice_max
#'
#' @importFrom flowCore identifier
#'
#' @export
slice_max.flowFrame <-
  function(
    .data,
    order_by,
    ...,
    n,
    prop,
    by = NULL,
    with_ties = TRUE,
    na_rm = FALSE
  ) {
    identifier <- flowCore::identifier(.data)
    tof_tibble <-
      .data |>
      as_tof_tbl(.name_method = "featureNames")
    result <-
      tof_tibble |>
      dplyr::slice_max(
        order_by = {{ order_by }},
        ...,
        n = n,
        prop = prop,
        by = {{ by }},
        with_ties = with_ties,
        na_rm = na_rm
      ) |>
      as_flowFrame()
    flowCore::identifier(result) <- identifier
    return(result)
  }

#' Subset rows of a data structure in order.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param order_by Variable or function of variables to order by.
#' To order by multiple variables, wrap them in a data frame or tibble.
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param with_ties Should ties be kept together? The default, TRUE, may return
#' more rows than you request. Use FALSE to ignore ties, and return the first
#' n rows.
#'
#' @param na_rm Should missing values in order_by be removed from the result?
#' If FALSE, NA values are sorted to the end so they will only be included if
#' there are insufficient non-missing values to reach n/prop.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice_max
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
slice_max.flowSet <-
  function(
    .data,
    order_by,
    ...,
    n,
    prop,
    by = NULL,
    with_ties = TRUE,
    na_rm = FALSE
  ){
    result_list <- list()
    for (i in 1:length(.data)) {
      identifier <- flowCore::identifier(.data[[i]])
      result_list[[identifier]] <-
        dplyr::slice_max(
          .data[[i]],
          order_by = {{ order_by }},
          ...,
          n = n,
          prop = prop,
          by = {{ by }},
          with_ties = with_ties,
          na_rm = na_rm
        )
    }
    result <- as(result_list, "flowSet")
    flowCore::pData(result) <- flowCore::pData(.data)
    return(result)
  }


### slice_min --------


#' Subset rows of a data structure in order.
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param order_by Variable or function of variables to order by.
#' To order by multiple variables, wrap them in a data frame or tibble.
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param with_ties Should ties be kept together? The default, TRUE, may return
#' more rows than you request. Use FALSE to ignore ties, and return the first
#' n rows.
#'
#' @param na_rm Should missing values in order_by be removed from the result?
#' If FALSE, NA values are sorted to the end so they will only be included if
#' there are insufficient non-missing values to reach n/prop.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} is preserved.
#'
#' @importFrom dplyr slice_min
#'
#' @importFrom flowCore identifier
#'
#' @export
slice_min.flowFrame <-
  function(
    .data,
    order_by,
    ...,
    n,
    prop,
    by = NULL,
    with_ties = TRUE,
    na_rm = FALSE
  ) {
    identifier <- flowCore::identifier(.data)
    tof_tibble <-
      .data |>
      as_tof_tbl()
    result <-
      tof_tibble |>
      dplyr::slice_min(
        order_by = {{ order_by }},
        ...,
        n = n,
        prop = prop,
        by = {{ by }},
        with_ties = with_ties,
        na_rm = na_rm
      ) |>
      as_flowFrame()
    flowCore::identifier(result) <- identifier
    return(result)
  }

#' Subset rows of a data structure in order.
#'
#' @param .data A \code{\link[flowCore]{flowSet}}
#'
#' @param order_by Variable or function of variables to order by.
#' To order by multiple variables, wrap them in a data frame or tibble.
#'
#' @param ... Unused.
#'
#' @param n,prop Provide either n, the number of rows, or prop, the proportion
#' of rows to select. If neither are supplied, n = 1 will be used.
#' If n is greater than the number of rows in the group (or prop > 1), the
#' result will be silently truncated to the group size. prop will be rounded
#' towards zero to generate an integer number of rows.
#'
#' A negative value of n or prop will be subtracted from the group size.
#' For example, n = -2 with a group of 5 rows will select 5 - 2 = 3 rows;
#' prop = -0.25 with 8 rows will select 8 * (1 - 0.25) = 6 rows.
#'
#' @param by Optionally, an unquoted selection of columns to group by for just this operation.
#' An alternative to group_by.
#'
#' @param with_ties Should ties be kept together? The default, TRUE, may return
#' more rows than you request. Use FALSE to ignore ties, and return the first
#' n rows.
#'
#' @param na_rm Should missing values in order_by be removed from the result?
#' If FALSE, NA values are sorted to the end so they will only be included if
#' there are insufficient non-missing values to reach n/prop.
#'
#' @returns An object of the same type as .data. The output has the following properties:
#' * Each row may appear 0, 1, or many times in the output.
#' * Columns are not modified.
#' * Groups are not modified.
#' * A \code{\link[flowCore]{flowSet}}'s \code{\link[flowCore]{pData}} is preserved.
#'
#' @importFrom dplyr slice_min
#'
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @export
slice_min.flowSet <-
  function(
    .data,
    order_by,
    ...,
    n,
    prop,
    by = NULL,
    with_ties = TRUE,
    na_rm = FALSE
  ) {
    result_list <- list()
    for (i in 1:length(.data)) {
      identifier <- flowCore::identifier(.data[[i]])
      result_list[[identifier]] <-
        dplyr::slice_min(
          .data[[i]],
          order_by = {{ order_by }},
          ...,
          n = n,
          prop = prop,
          by = {{ by }},
          with_ties = with_ties,
          na_rm = na_rm
        )
    }
    result <- as(result_list, "flowSet")
    flowCore::pData(result) <- flowCore::pData(.data)
    return(result)
  }


# tidyr methods ----------------------------------------------------------------

#' Nest a \code{\link[flowCore]{flowFrame}} into a \code{\link[flowCore]{flowSet}}
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Columns to nest; these will appear in the inner \code{\link[flowCore]{flowFrame}}s
#' comprising the output \code{\link[flowCore]{flowSet}}.
#' Specified using name-variable pairs of the form new_col = c(col1, col2, col3).
#' The right hand side can be any valid tidyselect expression.
#' If not supplied, then ... is derived as all columns not selected by .by.
#'
#' @param .by Columns to nest by; these will be stored in the
#' \code{\link[flowCore]{pData}} of the output \code{\link[flowCore]{flowSet}}.
#' .by can be used in place of or in conjunction with columns supplied through ....
#' If not supplied, then .by is derived as all columns not selected by ....
#'
#' @param .key Unused.
#'
#' @param .names_sep Unused.
#'
#' @returns A \code{\link[flowCore]{flowSet}} wherein cells are grouped into
#' constituent \code{\link[flowCore]{flowFrame}}s based on which columns are used
#' to nest.
#'
#' @importFrom rlang enquos
#' @importFrom rlang !!
#' @importFrom rlang !!!
#'
#' @export
#'
nest.flowFrame <- function(.data, ..., .by = NULL, .key = NULL, .names_sep = NULL) {

  .cols <- rlang::enquos(...)
  if (length(.cols) == 0) {
    result <-
      .data |>
      group_by.flowFrame({{ .by }})
  } else if (length(.cols) == 1) {
    names(.cols) <- NULL
    .cols <- .cols[[1]]
    result <-
      .data |>
      group_by.flowFrame(!! .cols)
  } else {
    names(.cols) <- NULL
    result <-
      .data |>
      group_by.flowFrame(!!! .cols)
  }

  return(result)
}

#' Unnest a \code{\link[flowCore]{flowSet}} into a single
#' \code{\link[flowCore]{flowFrame}}
#'
#' @param data A \code{\link[flowCore]{flowSet}}
#'
#' @param cols Columns in \code{\link[flowCore]{pData}} to unnest.
#'
#' @param ... Unused.
#'
#' @param keep_empty Unused.
#'
#' @param ptype Unused.
#'
#' @param names_sep Unused.
#'
#' @param names_repair Unused.
#'
#' @returns A \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}
#' depending on the degree of unnest-ing. Note that unnest-ing and
#' ungrouping a \code{\link[flowCore]{flowSet}} are equivalent.
#'
#' @export
unnest.flowSet <-
  function (
    data,
    cols,
    ...,
    keep_empty = FALSE,
    ptype = NULL,
    names_sep = NULL,
    names_repair = "check_unique"
  ) {
    result <-
      data |>
      ungroup.flowSet({{ cols }})

    return(result)
  }


# ggplot2 methods --------------------------------------------------------------

#' Create a new ggplot.
#'
#' @param data Default dataset to use for plot in the form of a
#' \code{\link[flowCore]{flowFrame}}. If not specified, must be supplied in each
#'  layer added to the plot.
#'
#' @param mapping Default list of aesthetic mappings to use for plot.
#' If not specified, must be supplied in each layer added to the plot. Note that
#' variable names used for aesthetic mappings come from the
#' \code{\link[flowCore]{featureNames}} of the input \code{\link[flowCore]{flowFrame}}.
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @param environment Deprecated. Used prior to tidy evaluation.
#'
#' @returns A \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#'
#' @export
#'
ggplot.flowFrame <-
  function(data = NULL, mapping = ggplot2::aes(), ..., environment = parent.frame()) {
    tof_tibble <-
      as_tof_tbl(data, .name_method = "featureNames")
    result <-
      ggplot2::ggplot(data = tof_tibble, mapping = mapping, ..., environment = environment)
    return(result)
  }

#' Create a new ggplot.
#'
#' @param data Default dataset to use for plot in the form of a
#' \code{\link[flowCore]{flowSet}}. If not specified, must be supplied in each
#'  layer added to the plot.
#'
#' @param mapping Default list of aesthetic mappings to use for plot.
#' If not specified, must be supplied in each layer added to the plot. Note that
#' variable names used for aesthetic mappings come from the
#' \code{\link[flowCore]{featureNames}} of the input
#' \code{\link[flowCore]{flowSet}}'s constituent \code{\link[flowCore]{flowFrame}}s.
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @param environment Deprecated. Used prior to tidy evaluation.
#'
#' @returns A \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#'
#' @export
ggplot.flowSet <-
  function(data = NULL, mapping = ggplot2::aes(), ..., environment = parent.frame()) {
    tof_tibble <-
      as_tof_tbl(data, .name_method = "featureNames", include_metadata = TRUE)
    result <-
      ggplot2::ggplot(data = tof_tibble, mapping = mapping, ..., environment = environment)
    return(result)
  }

