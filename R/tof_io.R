###    tof_read_fcs

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
    tof_tibble <-
      file_path %>%
      read.FCS(transformation = FALSE, truncate_max_range = FALSE)

    col_names <-
      if_else(
        are_na(as.character(tof_tibble@parameters@data$desc)),
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
    file_names <- list.files(file_path, full.names = FALSE)

    tof_tibble <-
      file_path %>%
      list.files(full.names = TRUE) %>%
      map(read.FCS, transformation = FALSE, truncate_max_range = FALSE)

    col_names <-
      if_else(
        are_na(as.character(tof_tibble[[1]]@parameters@data$desc)),
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
      map2(.y = file_names, .f = ~ mutate(.x, file_names = .y)) %>%
      #use data.table for the row binding, which can be quite slow otherwise
      data.table::rbindlist() %>%
      as_tibble()

  }
  return(tof_tibble)
}



###    tof_write_out

# Description:
### Reads an .fcs file (or a directory of .fcs files) into a tibble.
#
# Inputs:
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either
#       of these
#     - group_vars = variables that should be used to break cells into individual files
#     - out_folder_path = file path to the folder where the output files should be saved.
#     - format = format for the files being written. Currently supports .csv and .fcs files.
#     - Others?
#
# Outputs:
#     - none
#
# Side-effects:
#       - Save .csv or .fcs files in the designated folder.
#
# Dependencies:
#     - flowCore library
#     - tidyverse library
#     - data.table library

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
      tof_tibble %>%
        group_by(across(one_of(group_vars))) %>%
        nest(data = everything()) %>%
        ungroup() %>%
        mutate(prefix = str_c(setdiff(colnames(.), "data"))) %>%
        pmap(
          .l = .,
          .f =
            ~
            write_csv(
              x = ..1,
              path = file.path(out_folder_path, str_c(!!! group_vars, ".csv"))
            )
        )
    }
    if ("fcs" %in% format) {
      tof_tibble %>%
        group_by({group_vars}) %>%
        nest(data = everything()) %>%
        map(
          .x = data,
          .f = write.fcs,
          path = file.path(out_folder_path, str_c({group_vars}, ".fcs"))
        )
    }
  }

