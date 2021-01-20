############################    tof_read_fcs

# Description:
### Reads an .fcs file (or a director of .fcs files) into a tibble.
#
# Inputs:
#     - file_path = file path to a single .fcs file to be read into the R session.
#
# Outputs:
#     - tof_tibble = an [m x (n+1)] tibble in which each row represents a single cell
#                    (of m total in the dataset) and each column represents a metal
#                    measurement (of n total in the dataset). In addition, the last column
#                    of the tibble will represent the filename from which the cell was read.
#
# Dependencies:
#     - flowCore library
#     - tidyverse library
#     - data.table library

#optimization:
#              -add error if .fcs files have different panels
#

#' Read data from an .fcs file or a directory of .fcs files.
#'
#' @param file_path A file path. It can lead to either a single .fcs file or
#' to a directory filled with .fcs files.
#'
#' @return An [m by n+1] tibble in which each row represents a single cell (of m
#' total in the dataset) and each column represents a metal measurement
#' (of n total in the dataset). In addition, the last column of the tibble will
#' represent the filename from which the cell was read.
#'
#' @export
#'
#' @examples
#' NULL
tof_read_fcs <- function(file_path = NULL) {

  if (file_test(op = "-f", file_path)) {
    # if the path leads to a single file
    tof_tibble <-
      file_path %>%
      flowCore::read.FCS(transformation = FALSE, truncate_max_range = FALSE)

    col_names <-
      if_else(
        rlang::are_na(as.character(tof_tibble@parameters@data$desc)),
        tof_tibble@parameters@data$name,
        tof_tibble@parameters@data$desc
      ) %>%
      as.character()

    tof_tibble <-
      tof_tibble %>%
      {
        setNames(
          object = as_tibble(flowCore::exprs(.)),
          nm = col_names
        )
      }

  } else {
    # if the path leads to a directory
    file_names <- list.files(file_path, full.names = FALSE)

    tof_tibble <-
      file_path %>%
      list.files(full.names = TRUE) %>%
      map(read.FCS, transformation = FALSE, truncate_max_range = FALSE)

    col_names <-
      if_else(
        rlang::are_na(as.character(tof_tibble[[1]]@parameters@data$desc)),
        tof_tibble[[1]]@parameters@data$name,
        tof_tibble[[1]]@parameters@data$desc
      ) %>%
      as.character()

    tof_tibble <-
      tof_tibble %>%
      map(
        .f =
          ~ setNames(
            object = as_tibble((flowCore::exprs(.x))),
            nm = col_names
          )
      ) %>%
      map2(.y = file_names, .f = ~ mutate(.x, file_name = .y)) %>%
      #use data.table for the row binding, which can be quite slow otherwise
      data.table::rbindlist() %>%
      as_tibble()

  }
  return(tof_tibble)
}



############################    tof_write_out

#' Write data to a file
#'
#' Write data (in the form of a tibble) into either a .csv or an .fcs file for storage.
#'
#' @param tof_tibble a tibble, data.frame, or something that can be coerced into either
#' @param group_vars variables that should be used to group cells into individual files.
#' The default is to group cells based on the column "file_names" (as produced by tof_read_fcs).
#' @param out_folder_path file path to the folder where the output files should be saved
#' @param format format for the files being written. Currently supports .csv and .fcs files
#'
#' @return Nothing
#'
#' @export
#'
#' @examples
tof_write_out <-
  function(
    tof_tibble = NULL,
    group_vars = "file_names",
    out_folder_path = NULL,
    format = NULL
  ) {
    if ("csv" %in% format) {
      tof_tibble <-
        tof_tibble %>%
        group_by(across(one_of(group_vars))) %>%
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
            path = file.path(out_folder_path, str_c(!!! group_vars, ".csv"))
          )
      )
    }
    # has a bug
    if ("fcs" %in% format) {
      tof_tibble %>%
        group_by({group_vars}) %>%
        nest(data = everything()) %>%
        map(
          .x = data,
          .f = flowCore::write.fcs,
          path = file.path(out_folder_path, str_c({group_vars}, ".fcs"))
        )
    }
  }


############################    tof_save_figures

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
