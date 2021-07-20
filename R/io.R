# io.R
# This file contains functions relevant to reading input CyTOF data from files
# and writing output data to files.


# tof_find_panel_info ------------------

#' Use tkeyes's opinionated heuristic for extracted a CyTOF panel's metal-antigen pairs
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
    if_else(
      str_detect(data_names, pattern = str_c(metal_masterlist, collapse = "|")),
      str_extract(data_names, pattern = str_c(metal_masterlist, collapse = "|")),
      data_names
    )

  # find antigens
  antigens <-
    if_else(
      str_detect(data_desc, pattern = str_c(metal_masterlist, collapse = "|")),
      str_remove(data_desc, pattern = str_c(metal_masterlist, collapse = "|")),
      data_desc
    ) %>%
    str_remove("^[:punct:]|[:punct:]$") %>%
    str_remove_all("\\(|\\)|Di") %>%
    if_else(. == "", "empty", .)

  # return result
  result <-
    tibble(
      metals,
      antigens
    )

  return(result)
}





# tof_read_fcs ------------------

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
    col_names <- str_c(panel_info$antigens, panel_info$metals, sep = sep)

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
        panel = panel_info#,
        #nrow = nrow(tof_tibble),
        #class = "tof_tbl"
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
  function(file_path = NULL, panel_info = NULL) {

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

# tof_read_data ------------------

#' Read data from an .fcs/.csv file or a directory of .fcs/.csv files.
#'
#' @param path A file path to a single file or to a directory of files.
#' The only valid file types are .fcs files or .csv files
#' containing CyTOF data.
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
tof_read_data <- function(path = NULL, sep = "|") {

  # if the path leads to a single file
  if (file_test(op = "-f", x = path)) {
    if (get_extension(path) == "fcs") {
      tof_tibble <-
        path %>%
        tof_read_fcs()
    } else if (get_extension(path) == "csv") {
        tof_tibble <-
          path %>%
          tof_read_csv
    }
    return(tof_tibble)

    # if the path leads to a directory
  } else if (file_test(op = "-d", x = path)) {

    file_names <- list.files(path, full.names = FALSE)

    # read files into a nested tibble
    tof_tibble <-
      tibble::tibble(
        file_name = file_names,
        data = map(.x = list.files(path, full.names = TRUE), .f = tof_read_fcs)
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

# tof_write_csv ------------------

tof_write_csv <- function(tof_tibble, group_vars, out_path) {
  NULL
}



# tof_write_fcs ------------------

tof_write_fcs <- function(tof_tibble, group_vars, out_path) {
  NULL
}


# tof_write_out ------------------

#' Write cytof data to a file or directory of files
#'
#' Write data (in the form of a tibble) into either a .csv or an .fcs file for storage.
#'
#' @param tof_tibble A `tof_tibble`.
#'
#' @param group_vars Unquoted variable names representing columns to be used to
#' group cells into individual files. Defaults to `file_name` (as produced by
#' `tof_read_fcs`).
#'
#' @param out_path Path to the directory where output files should be saved.
#'
#' @param format format for the files being written. Currently supports .csv and .fcs files
#'
#' @return This function does not explicitly return any values. Instead,
#' it writes .csv or .fcs files to the specified `out_path`.
#'
#' @export
#'
#' @examples
#' NULL
tof_write_out <-
  function(
    tof_tibble = NULL,
    group_vars = file_name,
    out_path = NULL,
    format = c("csv", "fcs")
  ) {
    # check that the format argument is correctly specified
    format <- match.arg(arg = format, choices = c("csv", "fcs"), several.ok = TRUE)

    if ("csv" %in% format) {
      tof_tibble <-
        tof_tibble %>%
        group_by(across({{group_vars}})) %>%
        nest(data = everything()) %>%
        ungroup() %>%
        mutate(prefix = str_c(setdiff(colnames(.), "data")))

      #check this for bugs
      pwalk(
        .l = tof_tibble,
        .f =
          ~
          write_csv(
            x = ..1,
            path = file.path(out_path, str_c(!!! group_vars, ".csv"))
          )
      )
    }
    # has a bug
    if ("fcs" %in% format) {
      annotated_df <-
        Biobase::AnnotatedDataFrame(
          data = NULL,
          varMetadata = NULL
        )

      tof_tibble %>%
        group_by({group_vars}) %>%
        nest(data = everything()) %>%
        map(
          .x = data,
          .f = flowCore::write.fcs,
          path = file.path(out_path, str_c({group_vars}, ".fcs"))
        )
    }
  }


# tof_save_figures --------------------------

#' Save each entry in a list-column of ggplot figures in a tibble.
#'
#' Saves a list of plots contained in a plot_tibble (useful for nested tibbles).
#'
#' @param plot_tibble A tibble containing a column of ggplot objects
#' @param out_path File path to the directory in which the plots should be saved
#' @param label_column An unquoted column name representing a column in `plot_tibble`
#' containing strings that can be used to uniquely identify each plot in `plot_tibble`.
#' @param plot_column An unquoted column name indicating which column in `plot_tibble`
#' contains the plots to be saved.
#' @param label_prefixan Optional string to concatenate to the beginning of each file name.
#' @param label_suffix Optional string to concatenate to the end of each file name.
#' @param width The width for the plots being saved.
#' @param height The height for the plots being saved.
#' @param device "jpg", "tiff", or "pdf" file format.
#' @param ... Additional arguments to pass to `ggsave()`
#'
#' @return NULL
#'
#' @export
#'
#' @examples
tof_save_figures <-

  function(
    plot_tibble = NULL,
    out_path = NULL,
    label_column = NULL,
    plot_column = plots,
    label_prefix = "",
    label_suffix = "",
    width,
    height,
    device = "pdf",
    ...
  ) {
    walk2(
      .x = pull(plot_tibble, {{plot_column}}),
      .y = pull(plot_tibble, {{label_column}}),
      .f =
        ~ggsave(
          filename = str_c(label_prefix, .y, label_suffix, sep = "-"),
          plot = .x,
          device = device,
          path = out_path,
          width = width,
          height = height,
          units = "in"
        )
    )
  }
