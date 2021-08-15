# io.R
# This file contains functions relevant to reading input CyTOF data from files
# and writing output data to files.


# tof_find_panel_info ----------------------------------------------------------

#' Use tidytof's opinionated heuristic for extracted a CyTOF panel's metal-antigen pairs
#'
#' Using the character vectors obtained from the `name` and `desc` columns of
#' the parameters of the data of a flowFrame, figure out the CyTOF panel used
#' to collect the data and return it as a tidy tibble.
#'
#' @param input_flowFrame a raw flowFrame (just read from an .fcs file) from which
#' a CyTOF panel should be extracted
#'
#' @return A tibble with 2 columns (`metals` and `antigens`) that correspond to the
#' metals and antigens of the CyTOF panel used during data acquisition.
#'
#' @examples
#' NULL
#'
tof_find_panel_info <- function(input_flowFrame) {

  # extract "desc" and "names" information from the flowFrame
  data_desc <- as.character(input_flowFrame@parameters@data$desc)
  data_names <- as.character(input_flowFrame@parameters@data$name)

  # find metals
  metals <-
    dplyr::if_else(
      stringr::str_detect(
        data_names,
        pattern = str_c(metal_masterlist, collapse = "|")
      ),
      stringr::str_extract(
        data_names,
        pattern = str_c(metal_masterlist, collapse = "|")
      ),
      stringr::str_extract(
        data_desc,
        pattern = str_c(metal_masterlist, collapse = "|")
      )
    )

  # if no metal could be detected, just throw whatever was in the names
  # slot (to be as informative as possible)
  metals <-
    dplyr::if_else(
      is.na(metals),
      data_names,
      metals
    )

  # find antigens --------------------------------------------------------------

  # first, look in the description slot and remove any metal patterns. What
  # remains (minus any punctuation) is a candidate antigen name.
  antigens <-
    dplyr::if_else(
      stringr::str_detect(data_desc, pattern = str_c(metal_masterlist, collapse = "|")),
      stringr::str_remove(data_desc, pattern = str_c(metal_masterlist, collapse = "|")),
      data_desc
    ) %>%
    stringr::str_remove("^[:punct:]|[:punct:]$") %>%
    stringr::str_remove_all("\\(|\\)|Di")

  # if a given antigen name is empty after the first round of candidates is
  # explored, check the desc slot. Remove any metal patterns (and punctuation)
  # and what remains should be the antigen name.
  antigens <-
    dplyr::if_else(
      antigens == "",
      stringr::str_remove(data_names, pattern = str_c(metal_masterlist, collapse = "|")),
      antigens
    ) %>%
    stringr::str_remove("^[:punct:]|[:punct:]$") %>%
    stringr::str_remove_all("\\(|\\)|Di") %>%
    # if the antigen name of any given channel is still empty, just put the
    # word "empty"
    dplyr::if_else(. == "", "empty", .)

  # return result
  result <-
    tibble::tibble(
      metals,
      antigens
    )

  return(result)
}



# tof_read_fcs -----------------------------------------------------------------

#' Read CyTOF data from an .fcs file into a tidy tibble.
#'
#' This function reads CyTOF data from a single .fcs file into a tidy data
#' structure called a "tof_tibble." tof_tibbles are identical to normal
#' tibbles except for an additional attribute ("panel") that stores information
#' about the CyTOF panel used during data acquisition.
#'
#' @param file_path A file path to a single .fcs file.
#'
#' @param sep A string to use to separate the antigen name and its associated
#' metal in the column names of the output tibble. Defaults to "|".
#'
#' @return a `tof_tibble` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tibble` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition).
#'
#' @examples
#' NULL
#'
tof_read_fcs <-
  function(file_path = NULL, sep = "|") {

    # read flowFrame from file
    tof_flowFrame <-
      file_path %>%
      flowCore::read.FCS(transformation = FALSE, truncate_max_range = FALSE)

    # extract panel information from inner parameters of the flowFrame
    panel_info <- tof_find_panel_info(input_flowFrame = tof_flowFrame)

    # derive and set the column names to use for the output tof_tibble
    col_names <-
      stringr::str_c(panel_info$antigens, panel_info$metals, sep = sep)

    tof_tibble <-
      tof_flowFrame %>%
      {
        setNames(
          object = as_tibble(flowCore::exprs(.)),
          nm = col_names
        )
      }

    tof_tibble <-
      new_tof_tibble(
        x = tof_tibble,
        panel = panel_info
      )

    return(tof_tibble)
  }

# tof_read_csv ------------------

#' Read CyTOF data from a .csv file into a tidy tibble.
#'
#' @param file_path A file path to a single .csv file.
#'
#' @param panel_info Optional. A tibble or data.frame containing information about the
#' panel used during CyTOF data acquisition. Two columns are required:
#' "metals" and "antigens".
#'
#' @return A `tof_tibble` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tibble` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition). Because panel information isn't
#' obvious from data read as a .csv file, this information must be provided
#' manually from the user (unlike in `tof_read_fcs`).
#'
#' @export
#'
#' @examples
#' NULL
tof_read_csv <-
  function(file_path = NULL, panel_info = tibble::tibble()) {

    tof_tibble <-
      file_path %>%
      readr::read_csv(file_path)

    # check that panel_info typing is correct
    if (is.data.frame(panel_info)) {
      panel_info <- tibble::as_tibble(panel_info)
    } else if (!is.null(panel-info)) {
      stop("panel_info must be a tibble, a data.frame, or NULL")
    }

    tof_tibble <-
      tibble::new_tibble(
        x = tof_tibble,
        panel = panel_info,
        nrow = nrow(tof_tibble),
        class = "tof_tbl"
      )

    return(tof_tibble)
  }


#' Read CyTOF data from a single .fcs or .csv file into a tidy tibble.
#'
#' @param file_path A file path to a single .fcs or .csv file.
#'
#' @param sep A string to use to separate the antigen name and its associated
#' metal in the column names of the output tibble. Defaults to "|". Only used
#' if the input file is an .fcs file.
#'
#' @param panel_info Optional. A tibble or data.frame containing information about the
#' panel used during CyTOF data acquisition. Two columns are required:
#' "metals" and "antigens". Only used if the input file is a .csv file.
#'
#' @return A `tof_tibble` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tibble` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition). Because panel information isn't
#' obvious from data read as a .csv file, this information must be provided
#' manually from the user.
#'
#' @examples
#' NULL
#'
tof_read_file <- function(file_path = NULL, sep = "|", panel_info = NULL) {
  if (get_extension(file_path) == "fcs") {
    tof_tibble <-
      file_path %>%
      tof_read_fcs(sep = sep)
  } else if (get_extension(file_path) == "csv") {
    tof_tibble <-
      file_path %>%
      tof_read_csv(panel_info = panel_info)
  }
  return(tof_tibble)
}

# tof_read_data ----------------------------------------------------------------
# Notes: Will be buggy in the event that a direcory has a combination of .fcs
# and .csv files, as their "panel" attribute may differ (and currently we do
# not do anything to check for this or fix it).


#' Read data from an .fcs/.csv file or a directory of .fcs/.csv files.
#'
#' @param path A file path to a single file or to a directory of files.
#' The only valid file types are .fcs files or .csv files
#' containing CyTOF data.
#'
#' @param sep Optional. A string to use to separate the antigen name and its associated
#' metal in the column names of the output tibble. Defaults to "|". Only used if
#' the input file is an .fcs file.
#'
#' @param panel_info Optional. A tibble or data.frame containing information about the
#' panel used during CyTOF data acquisition. Two columns are required:
#' "metals" and "antigens". Only used if the input file is a .csv file.
#'
#' @return An [c by m+1] tibble in which each row represents a single cell (of c
#' total in the dataset) and each column represents a CyTOF measurement
#' (of m total in the dataset). If more than one .fcs is read at once,
#' the last column of the tibble (`file_name`) will represent the file name
#' of the .fcs file from which each cell was read.
#'
#' @export
#'
#' @examples
#' NULL
#'
#'
tof_read_data <- function(path = NULL, sep = "|", panel_info = tibble::tibble()) {

  # if the path leads to a single file
  if (file_test(op = "-f", x = path)) {
    tof_tibble <- tof_read_file(path, sep = sep, panel_info = panel_info)
    return(tof_tibble)

    # if the path leads to a directory
  } else if (file_test(op = "-d", x = path)) {

    file_names <-
      list.files(path, full.names = FALSE) %>%
      stringr::str_remove_all("\\.fcs|\\.csv")

    # read files into a nested tibble
    tof_tibble <-
      tibble::tibble(
        file_name = file_names,
        data =
          map(
            .x = list.files(path, full.names = TRUE),
            .f = tof_read_file,
            sep = sep,
            panel_info = panel_info
          )
      )

    # check how many unique panels are present across all files being read
    panels <-
      map(.x = tof_tibble$data, ~attr(x = .x, which = "panel"))

    num_panels <-
      panels %>%
      unique() %>%
      length()

    if (num_panels > 1) {
      # group by panel
      tof_tibble <-
        tof_tibble %>%
        mutate(panel = panels) %>%
        nest(data = -panel) %>%
        mutate(
          data =
            map2(
              .x = data,
              .y = panel,
              .f = ~new_tof_tibble(x = .x, panel = .y)
            )
        ) %>%
        mutate(
          data =
            map2(
              .x = data,
              .y = panel,
              .f = ~
                new_tof_tibble(x = unnest(.x, cols = data), panel = .y)
            )
        )

    } else {
      # put everything together
      tof_tibble <-
        tof_tibble %>%
        unnest(cols = data)
    }

  }

  return(tof_tibble)
}

# tof_write_csv ----------------------------------------------------------------

#' Write a series of .csv files from a tof_tibble
#'
#' This function takes a given `tof_tibble` and writes the single-cell data
#' it contains into .csv files within the directory located at `out_path`. The
#' `group_vars` argument specifies how the rows of the `tof_tibble` (each cell)
#' should be broken into separate .csv files
#'
#' @param tof_tibble A `tof_tibble`.
#' @param group_vars Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults
#' to selecting all non-numeric (i.e. non-integer and non-double) columns.
#' @param out_path A system path indicating the directory where the output .csv
#' files should be saved. If the directory doesn't exist, it will be created.
#' @param sep Delimiter that should be used between each of the values of `group_vars`
#' to create the output .csv file names. Defaults to "_".
#'
#' @return This function does not return anything. Instead, it has the side-effect
#' of saving .csv files to `out_path`.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_write_csv <-
  function(
    tof_tibble,
    group_vars = where(~ !(purrr::is_integer(.x) || purrr::is_double(.x))),
    out_path,
    sep = "_"
  ) {

    dir.create(path = out_path, showWarnings = FALSE, recursive = TRUE)

    tof_tibble <-
      tof_tibble %>%
      group_by(across({{group_vars}})) %>%
      nest() %>%
      ungroup() %>%
      unite(col = "prefix", -data, sep = sep)

    walk2(
      .x = tof_tibble$prefix,
      .y = tof_tibble$data,
      .f = ~
        readr::write_csv(
          x = .y,
          file = file.path(out_path, stringr::str_c(.x, ".csv"))
        )
    )
  }



# tof_write_fcs ----------------------------------------------------------------

#' Write a series of .fcs files from a tof_tibble
#'
#' This function takes a given `tof_tibble` and writes the single-cell data
#' it contains into .fcs files within the directory located at `out_path`. The
#' `group_vars` argument specifies how the rows of the `tof_tibble` (each cell)
#' should be broken into separate .fcs files
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param group_vars Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults
#' to selecting all non-numeric (i.e. non-integer and non-double) columns.
#'
#' @param out_path A system path indicating the directory where the output .csv
#' files should be saved. If the directory doesn't exist, it will be created.
#'
#' @param sep Delimiter that should be used between each of the values of `group_vars`
#' to create the output .fcs file names. Defaults to "_".
#'
#' @return This function does not return anything. Instead, it has the side-effect
#' of saving .fcs files to `out_path`.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_write_fcs <-
  function(
    tof_tibble,
    group_vars = where(~ !(tof_is_numeric(.x))),
    out_path,
    sep = "_"
  ) {

    # create out_path
    dir.create(path = out_path, showWarnings = FALSE, recursive = TRUE)

    # eliminate all non-grouping and non-numeric columns from tof_tibble
    tof_tibble <-
      tof_tibble %>%
      dplyr::select({{group_vars}}, where(tof_is_numeric))

    # find max and min values for all non-grouping columns in tof_tibble
    maxes_and_mins <-
      tof_tibble %>%
      dplyr::summarize(
        dplyr::across(
          -{{group_vars}},
          .fns = list(max = ~max(.x, na.rm = TRUE), min= ~min(.x, na.rm = TRUE)),
          # use the many underscores because it's unlikely this will come up
          # in column names on their own...maybe make more rigorous?
          .names = "{.col}_____{.fn}"
        )
      ) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = c("antigen", "value_type"),
        values_to = "value",
        names_sep = "_____"
      )  %>%
      tidyr::pivot_wider(
        names_from = value_type,
        values_from = value
      )

    # extract the names of all non-grouping columns to be saved to the .fcs file
    data_cols <- maxes_and_mins$antigen

    # nest tof_tibble
    tof_tibble <-
      tof_tibble %>%
      dplyr::group_by(across({{group_vars}})) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      tidyr::unite(col = "prefix", -data, sep = sep)

    # make components of parameters AnnotatedDataFrame
    fcs_varMetadata <-
      data.frame(
        labelDescription =
          c(
            "Name of Parameter",
            "Description of Parameter",
            "Range of Parameter",
            "Minimum Parameter Value after Transformation",
            "Maximum Parameter Value after Transformation"
          )
      )

    fcs_data <-
      maxes_and_mins %>%
      dplyr::transmute(
        # have to change any instances of "|" in column names to another
        # separator, as "|" has special meaning as an .fcs file delimiter
        name = sstringr::tr_replace(antigen, "\\|", "_"),
        desc = string::str_replace(antigen, "\\|", "_"),
        range = max - min,
        minRange = min,
        maxRange = max
      ) %>%
      as.data.frame()

    row.names(fcs_data) <- stringr::str_c("$", "P", 1:nrow(fcs_data))

    # make the AnnotatedDataFrame
    parameters <-
      Biobase::AnnotatedDataFrame(data = fcs_data, varMetadata = fcs_varMetadata)

    # make flowFrames for each row of tof_tibble
    tof_tibble <-
      tof_tibble %>%
      dplyr::transmute(
        prefix,
        flowFrames =
          purrr:::map(
            .x = data,
            ~ flowCore::flowFrame(
              exprs =
                # have to change any instances of "|" in column names to another
                # separator, as "|" has special meaning as an .fcs file delimiter
                as.matrix(
                  dplyr::rename_with(
                    .x,
                    stringr::str_replace,
                    pattern = "\\|",
                    replacement = "_"
                  )
                ),
              parameters = parameters
            )
          )
      )

    # write out final .fcs files
    purrr::walk2(
      .x = tof_tibble$prefix,
      .y = tof_tibble$flowFrames,
      .f = ~
        flowCore::write.FCS(
          x = .y,
          filename = file.path(out_path, stringr::str_c(.x, ".fcs"))
        )
    )

  }


# tof_write_data ---------------------------------------------------------------

#' Write cytof data to a file or to a directory of files
#'
#' Write data (in the form of a tof_tibble) into either a .csv or an .fcs file for storage.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param group_vars group_vars Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults
#' to selecting all non-numeric (i.e. non-integer and non-double) columns.
#'
#' @param out_path Path to the directory where output files should be saved.
#'
#' @param format format for the files being written. Currently supports .csv and .fcs files
#'
#' @param sep Delimiter that should be used between each of the values of `group_vars`
#' to create the output .csv/.fcs file names. Defaults to "_".
#'
#' @return This function does not explicitly return any values. Instead,
#' it writes .csv or .fcs files to the specified `out_path`.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_write_data <-
  function(
    tof_tibble = NULL,
    group_vars = where(~ !(tof_is_numeric(.x))),
    out_path = NULL,
    format = c("csv", "fcs"),
    sep = "_"
  ) {
    # check that the format argument is correctly specified
    format <- rlang::arg_match(arg = format)

    # if .csv file is asked for
    if ("csv" %in% format) {
      tof_write_csv(
        tof_tibble = tof_tibble,
        group_vars = {{group_vars}},
        out_path = out_path,
        sep = sep
      )
    }

    # if .fcs file is asked for
    if ("fcs" %in% format) {
      tof_write_fcs(
        tof_tibble = tof_tibble,
        group_vars = {{group_vars}},
        out_path = out_path,
        sep = sep
      )
    }
  }

