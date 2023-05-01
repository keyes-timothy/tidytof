# io.R
# This file contains functions relevant to reading input CyTOF data from files
# and writing output data to files.

# tof_find_panel_info ----------------------------------------------------------

#' Use tidytof's opinionated heuristic for extracted a CyTOF panel's metal-antigen pairs
#' from a flowFrame (read from a .fcs file.)
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
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom stringr str_remove
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_c
#'
#' @importFrom dplyr if_else
#' @importFrom dplyr tibble
#'
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
        pattern = stringr::str_c(metal_masterlist, collapse = "|")
      ),
      stringr::str_c(
        stringr::str_extract(
          data_names,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) %>%
          stringr::str_extract("[:alpha:]+"),
        stringr::str_extract(
          data_names,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) %>%
          stringr::str_extract("[:digit:]+")
      ),
      stringr::str_c(
        stringr::str_extract(
          data_desc,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) %>%
          stringr::str_extract("[:alpha:]+"),
        stringr::str_extract(
          data_desc,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) %>%
          stringr::str_extract("[:digit:]+")
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
      stringr::str_detect(data_desc, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      stringr::str_remove(data_desc, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      data_desc
    ) %>%
    stringr::str_remove("^[:punct:]|[:punct:]$") %>%
    stringr::str_remove_all("\\(|\\)|Di")

  # if a given antigen name is empty after the first round of candidates is
  # explored, check the names slot. Remove any metal patterns (and punctuation)
  # and what remains should be the antigen name.
  antigens <-
    dplyr::if_else(
      antigens == "" | is.na(antigens),
      stringr::str_remove(data_names, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      antigens
    ) %>%
    stringr::str_remove("^[:punct:]|[:punct:]$") %>%
    stringr::str_remove_all("\\(|\\)|Di")

  # if the antigen name of any given channel is still empty (or NA), just put
  # the word "empty"
  antigens <-
    dplyr::if_else((antigens == "" | is.na(antigens)), "empty", antigens)

  # return result
  result <-
    dplyr::tibble(
      metals,
      antigens
    )

  return(result)
}



# tof_read_fcs -----------------------------------------------------------------

#' Read CyTOF data from an .fcs file into a tidy tibble.
#'
#' This function reads CyTOF data from a single .fcs file into a tidy data
#' structure called a `tof_tbl` ("tof_tibble"). tof_tibbles are identical to normal
#' tibbles except for an additional attribute ("panel") that stores information
#' about the CyTOF panel used during data acquisition.
#'
#' @param file_path A file path to a single .fcs file.
#'
#' @param sep A string to use to separate the antigen name and its associated
#' metal in the column names of the output tibble. Defaults to "|".
#'
#' @return a `tof_tbl` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tbl` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition).
#'
#' @importFrom dplyr if_else
#' @importFrom dplyr as_tibble
#'
#' @importFrom flowCore read.FCS
#'
#' @importFrom utils capture.output
#'
tof_read_fcs <-
  function(file_path = NULL, sep = "|") {

    # read flowFrame from file
    invisible(
      utils::capture.output(
        tof_flowFrame <-
          file_path %>%
          flowCore::read.FCS(transformation = FALSE, truncate_max_range = FALSE)
      )
    )

    # extract panel information from inner parameters of the flowFrame
    panel_info <- tof_find_panel_info(input_flowFrame = tof_flowFrame)

    # derive and set the column names to use for the output tof_tibble
    col_names <-
      # if the antigen and the metal name are the same,
      # use a simplified column name
      dplyr::if_else(
        panel_info$antigens == panel_info$metals,
        panel_info$antigens,
        base::paste(panel_info$antigens, panel_info$metals, sep = sep)
      )

    tof_tibble <-
      setNames(
        object = dplyr::as_tibble(flowCore::exprs(tof_flowFrame)),
        nm = col_names
      )

    tof_tibble <-
      new_tof_tibble(
        x = tof_tibble,
        panel = panel_info
      )

    return(tof_tibble)
  }

# tof_read_csv -----------------------------------------------------------------

#' Read CyTOF data from a .csv file into a tidy tibble.
#'
#' @param file_path A file path to a single .csv file.
#'
#' @param panel_info Optional. A tibble or data.frame containing information about the
#' panel used during CyTOF data acquisition. Two columns are required:
#' "metals" and "antigens".
#'
#' @return A `tof_tbl` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tbl` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition). Because panel information isn't
#' obvious from data read as a .csv file, this information must be provided
#' manually from the user (unlike in `tof_read_fcs`).
#'
#' @importFrom dplyr as_tibble
#' @importFrom dplyr tibble
#' @importFrom dplyr select
#'
#' @importFrom readr read_csv
#' @importFrom readr cols
#'
#'
tof_read_csv <-
  function(file_path = NULL, panel_info = dplyr::tibble()) {

    tof_tibble <-
      file_path %>%
      readr::read_csv(col_types = readr::cols(), progress = FALSE)

    # check that panel_info typing is correct
    if (is.data.frame(panel_info)) {
      panel_info <- dplyr::as_tibble(panel_info)

      panel_names <- colnames(panel_info)

      if(all(c("antigens", "metals") %in% panel_names)) {
        panel_info <-
          panel_info %>%
          dplyr::select("antigens", "metals")
      } else if (!identical(panel_info, dplyr::tibble())) {
        stop("panel_info must contain an \"antigens\" and a \"metals\" column")
      }

    } else if (!is.null(panel_info)) {
      stop("panel_info must be a tibble, a data.frame, or NULL")
    }

    tof_tibble <-
      new_tof_tibble(
        x = tof_tibble,
        panel = panel_info
      )

    return(tof_tibble)
  }


# tof_read_file ----------------------------------------------------------------

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
#' @return A `tof_tbl` in which each row represents a single cell and each
#' column represents a CyTOF antigen channel.
#'
#' A `tof_tbl` is an S3 class that extends the "tibble" class by storing
#' one additional attribute: "panel" (a tibble storing information about the
#' panel used during data acquisition). Because panel information isn't
#' obvious from data read as a .csv file, this information must be provided
#' manually by the user.
#'
#'
tof_read_file <- function(file_path = NULL, sep = "|", panel_info = dplyr::tibble()) {
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
#' @family input/output functions
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr pluck
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @importFrom stringr str_remove_all
#'
#' @examples
#'
#' input_file <- dir(tidytof_example_data("aml"), full.names = TRUE)[[1]]
#' tof_read_data(input_file)
#'
tof_read_data <- function(path = NULL, sep = "|", panel_info = dplyr::tibble()) {

  # if path gives multiple paths to multiple files
  if (length(path) > 1) {
    if (all(file_test(op = "-f", x = path))) {
      file_paths <- path
      file_names <- basename(file_paths)
    } else {
      stop("If length(path) > 1, each entry must point to a single file.")
    }

    # if path only gives a single path
  } else if (length(path == 1)) {

    # if the path leads to a single file
    if (file_test(op = "-f", x = path)) {
      tof_tibble <- tof_read_file(path, sep = sep, panel_info = panel_info)
      return(tof_tibble)

      # if the path leads to a directory
    } else if (file_test(op = "-d", x = path)) {
      file_paths <-
        list.files(path, full.names = TRUE)
      file_names <- basename(file_paths)
    } else {
      stop(
        "path must be a vector of paths to data files or a
        single path to a data file or directory."
      )
    }
  }

  # read files into a nested tibble
  tof_tibble <-
    dplyr::tibble(
      file_name = file_names,
      data =
        purrr::map(
          .x = file_paths,
          .f = tof_read_file,
          sep = sep,
          panel_info = panel_info
        )
    )

  # check how many unique panels are present across all files being read
  panels <-
    purrr::map(.x = tof_tibble$data, ~attr(x = .x, which = "panel"))

  num_panels <-
    panels %>%
    unique() %>%
    length()

  if (num_panels > 1) {
    # group by panel
    tof_tibble <-
      tof_tibble %>%
      dplyr::mutate(panel = panels) %>%
      tidyr::nest(data = -panel) %>%
      dplyr::mutate(
        data =
          purrr::map2(
            .x = data,
            .y = panel,
            .f = ~new_tof_tibble(x = .x, panel = .y)
          )
      ) %>%
      mutate(
        data =
          purrr::map2(
            .x = data,
            .y = panel,
            .f = ~
              new_tof_tibble(x = tidyr::unnest(.x, cols = data), panel = .y)
          )
      )

  } else {
    # put everything together
    tof_tibble <-
      tof_tibble %>%
      tidyr::unnest(cols = data)

    panel <-
      panels %>%
      unique()
    panel <- panel[[1]]

    tof_tibble <-
      new_tof_tibble(x = tof_tibble, panel = panel)
  }


  return(tof_tibble)
}

# tof_write_csv ----------------------------------------------------------------

#' Write a series of .csv files from a tof_tbl
#'
#' This function takes a given `tof_tbl` and writes the single-cell data
#' it contains into .csv files within the directory located at `out_path`. The
#' `group_cols` argument specifies how the rows of the `tof_tbl` (each cell)
#' should be broken into separate .csv files
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#' @param group_cols Optional. Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults to
#' NULL (all cells are written into a single file).
#' @param out_path A system path indicating the directory where the output .csv
#' files should be saved. If the directory doesn't exist, it will be created.
#' @param sep Delimiter that should be used between each of the values of `group_cols`
#' to create the output .csv file names. Defaults to "_".
#' @param file_name If `group_cols` isn't specified, the name (without an extension)
#' that should be used for the saved .csv file.
#'
#' @return This function does not return anything. Instead, it has the side-effect
#' of saving .csv files to `out_path`.
#'
#' @family input/output functions
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unite
#' @importFrom purrr walk2
#' @importFrom readr write_csv
#' @importFrom stringr str_c
#' @importFrom stringr str_replace
#'
tof_write_csv <-
  function(
    tof_tibble,
    group_cols,
    out_path,
    sep = "_",
    file_name
  ) {

    # create the output directory if it doesn't already exist
    dir.create(path = out_path, showWarnings = FALSE, recursive = TRUE)

    # if groups_cols is NULL, make sure there's a filename
    if (missing(group_cols) & missing(file_name)) {
      stop("if `group_cols` are not provided, you must specify a `file_name.`")
    } else if (missing(group_cols)) {
      file_name <- stringr::str_replace(file_name, "\\.csv$", "")
    }

    # nest cells using group_cols if they are provided, otherwise nest
    # all cells into a single tof_tibble
    if (missing(group_cols)) {
      tof_tibble <-
        suppressWarnings(
          tof_tibble %>%
            tidyr::nest() %>%
            dplyr::ungroup()
        ) %>%
        dplyr::mutate(prefix = file_name)

    } else {
      tof_tibble <-
        suppressWarnings(
          tof_tibble %>%
            dplyr::group_by(dplyr::across({{group_cols}})) %>%
            tidyr::nest() %>%
            dplyr::ungroup()
        ) %>%
        tidyr::unite(col = "prefix", -data, sep = sep)
    }

    # save each tof_tbl in the nested tibble as its own file
    purrr::walk2(
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

#' Write a series of .fcs files from a tof_tbl
#'
#' This function takes a given `tof_tbl` and writes the single-cell data
#' it contains into .fcs files within the directory located at `out_path`. The
#' `group_cols` argument specifies how the rows of the `tof_tbl` (each cell)
#' should be broken into separate .fcs files
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults to
#' NULL (all cells are written into a single file).
#'
#' @param out_path A system path indicating the directory where the output .csv
#' files should be saved. If the directory doesn't exist, it will be created.
#'
#' @param sep Delimiter that should be used between each of the values of `group_cols`
#' to create the output .fcs file names. Defaults to "_".
#'
#' @param file_name If `group_cols` isn't specified, the name (without an extension)
#' that should be used for the saved .csv file.
#'
#' @return This function does not return anything. Instead, it has the side-effect
#' of saving .fcs files to `out_path`.
#'
#' @family input/output functions
#'
#' @export
#'
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#'
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unite
#'
#' @importFrom tidyselect everything
#'
#' @importFrom stringr str_c
#' @importFrom stringr str_replace
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowCore write.FCS
#'
#' @importFrom purrr walk2
#'
#' @importFrom methods new
#'
tof_write_fcs <-
  function(
    tof_tibble,
    group_cols,
    out_path,
    sep = "_",
    file_name
  ) {

    # create out_path
    dir.create(path = out_path, showWarnings = FALSE, recursive = TRUE)

    # if groups_cols is NULL, make sure there's a filename
    if (missing(group_cols) & missing(file_name)) {
      stop("if `group_cols` are not provided, you must specify a `file_name.`")
    } else if (missing(group_cols)) {
      file_name <- stringr::str_replace(file_name, "\\.fcs$", "")
    }

    # eliminate all non-grouping and non-numeric columns from tof_tibble
    if (!missing(group_cols)) {
    tof_tibble <-
      tof_tibble %>%
      dplyr::select({{group_cols}}, where(tof_is_numeric))
    } else {
      tof_tibble <-
        tof_tibble %>%
        dplyr::select(where(tof_is_numeric))
    }

    # find max and min values for all non-grouping columns in tof_tibble
    if (!missing(group_cols)) {
      maxes_and_mins <-
        tof_tibble %>%
        dplyr::summarize(
          dplyr::across(
            -{{group_cols}},
            .fns =
              list(max = ~ max(.x, na.rm = TRUE), min = ~ min(.x, na.rm = TRUE)),
            # use the many underscores because it's unlikely this will come up
            # in column names on their own
            .names = "{.col}_____{.fn}"
          )
        )
    } else {
      maxes_and_mins <-
        tof_tibble %>%
        dplyr::summarize(
          dplyr::across(
            dplyr::everything(),
            .fns =
              list(max = ~ max(.x, na.rm = TRUE), min= ~ min(.x, na.rm = TRUE)),
            # use the many underscores because it's unlikely this will come up
            # in column names on their own
            .names = "{.col}_____{.fn}"
          )
        )
    }

    maxes_and_mins <-
      maxes_and_mins %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = c("antigen", "value_type"),
        values_to = "value",
        names_sep = "_____"
      )  %>%
      tidyr::pivot_wider(
        names_from = "value_type",
        values_from = "value"
      )

    # extract the names of all non-grouping columns to be saved to the .fcs file
    data_cols <- maxes_and_mins$antigen

    # nest cells using group_cols if they are provided, otherwise nest
    # all cells into a single tof_tibble
    if (missing(group_cols)) {
      tof_tibble <-
        suppressWarnings(
          tof_tibble %>%
            tidyr::nest() %>%
            dplyr::ungroup()
        ) %>%
        dplyr::mutate(prefix = file_name)

    } else {
      tof_tibble <-
        suppressWarnings(
          tof_tibble %>%
            dplyr::group_by(dplyr::across({{group_cols}})) %>%
            tidyr::nest() %>%
            dplyr::ungroup()
        ) %>%
        tidyr::unite(col = "prefix", -data, sep = sep)
    }

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
        name = stringr::str_replace(.data$antigen, "\\|", "_"),
        desc = stringr::str_replace(.data$antigen, "\\|", "_"),
        range = max - min,
        minRange = min,
        maxRange = max
      ) %>%
      as.data.frame()

    row.names(fcs_data) <- stringr::str_c("$", "P", 1:nrow(fcs_data))

    # make the AnnotatedDataFrame
    parameters <-
      methods::new(
        "AnnotatedDataFrame",
        data = fcs_data,
        varMetadata = fcs_varMetadata
      )

    # make flowFrames for each row of tof_tibble
    tof_tibble <-
      tof_tibble %>%
      dplyr::transmute(
        .data$prefix,
        flowFrames =
          purrr::map(
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
#' Write data (in the form of a `tof_tbl`) into either a .csv or an .fcs file for storage.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param group_cols Optional. Unquoted names of the columns in `tof_tibble` that should
#' be used to group cells into separate files. Supports tidyselect helpers. Defaults to
#' no grouping (all cells are written into a single file).
#'
#' @param out_path Path to the directory where output files should be saved.
#'
#' @param format format for the files being written. Currently supports .csv and .fcs files
#'
#' @param sep Delimiter that should be used between each of the values of `group_cols`
#' to create the output .csv/.fcs file names. Defaults to "_".
#'
#' @param file_name If `group_cols` isn't specified, the name (without an extension)
#' that should be used for the saved file.
#'
#' @return This function does not explicitly return any values. Instead,
#' it writes .csv and/or .fcs files to the specified `out_path`.
#'
#' @family input/output functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#'
tof_write_data <-
  function(
    tof_tibble = NULL,
    group_cols,
    out_path = NULL,
    format = c("fcs", "csv"),
    sep = "_",
    file_name
  ) {
    # check that the format argument is correctly specified
    format <- rlang::arg_match(arg = format)

    # if .csv file is requested
    if ("csv" %in% format) {
      if (!missing(group_cols)) {
        tof_write_csv(
          tof_tibble = tof_tibble,
          group_cols = {{group_cols}},
          out_path = out_path,
          sep = sep,
          file_name = file_name
        )
      } else {
        tof_write_csv(
          tof_tibble = tof_tibble,
          out_path = out_path,
          sep = sep,
          file_name = file_name
        )
      }
    }

    # if .fcs file is requested
    if ("fcs" %in% format) {
      if (!missing(group_cols)) {
        tof_write_fcs(
          tof_tibble = tof_tibble,
          group_cols = {{group_cols}},
          out_path = out_path,
          sep = sep,
          file_name = file_name
        )
      } else {
        tof_write_fcs(
          tof_tibble = tof_tibble,
          out_path = out_path,
          sep = sep,
          file_name = file_name
        )
      }
    }
  }

