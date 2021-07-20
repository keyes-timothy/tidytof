
# tidytof_example_data ------------------

#' Get paths to tidytof example data
#'
#' tidytof comes bundled with a number of sample .fcs files in its
#' inst/extdata directory. This function makes them easy to access.
#'
#' @param dataset_name Name of the dataset you want to access. If NULL,
#' the names of the datasets (each of which is from a different study)
#' will be listed.
#'
#' @return A character vector of file paths where the requested .fcs
#' files are located. If `dataset_name` is NULL, a character vector of
#' dataset names (that can be used as values for `dataset_name`) is
#' returned instead.
#'
#' @export
#'
#' @examples
#' tidytof_example_data()
#' tidytof_example_data(dataset_name = "phenograph")
#'
tidytof_example_data <-
  function(dataset_name = NULL) {

    if (is.null(dataset_name)) {
      dir(system.file("extdata", package = "tidytof"))
    }
    else {
      dir(
        system.file("extdata", dataset_name, package = "tidytof", mustWork = TRUE),
        full.names = TRUE
      )
    }
  }


# get_extension ------------------

#' Find the extension for a file
#'
#' @param filename A string representing the name of a file in its local directory
#'
#' @return The the file extension of `filename`
#'
#' @examples
#' \dontrun{
#' # example file name
#' my_filename <- "my_file.txt"
#'
#' # find and print the extension
#' my_extension <- getExtension(my_filename)
#' print(my_extension)
#' }
get_extension <- function(filename) {
  ex <- strsplit(basename(filename), split="\\.")[[1]]
  return(ex[[length(ex)]])
}

# new_tof_tibble ------------------

#' Create a new tof_tibble.
#'
#' @param x
#'
#' @param panel
#'
#' @return
#'
#' @examples
#'
new_tof_tibble <- function(x = tibble::tibble(), panel = tibble::tibble()) {

  stopifnot(tibble::is_tibble(x))
  stopifnot(tibble::is_tibble(panel))

  # is this hack-y? I want tof_tbl's to inherit behavior from both
  # data.frame's, tbl's, and tbl_df's.
  # TO DO: implement certain methods (like unnest) that will
  # deal with panels in a way that makes sense
  vctrs::new_data_frame(
    x,
    panel = panel,
    n = nrow(x),
    class = c("tof_tbl", "tbl_df", "tbl")
  )
}

# tof_get_panel ------------------

#' Get panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @return A tibble containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in `tof_tibble`.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_get_panel <- function(tof_tibble) {
  panel <-
    tof_tibble %>%
    attr(which = "panel")

  return(panel)
}


# tof_set_panel -------------------

#' Set panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @param panel A tibble containing two columns representing the
#' information about a panel
#'
#' @return A `tof_tibble` containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in the input `tof_tibble`. Two columns are required: "metals" and "antigens".
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_set_panel <- function(tof_tibble, panel) {
  attr(tof_tibble, which = "panel") <- panel
  return(tof_tibble)
}


# tof_tbl methods --------------------------------------------

#' @export
nest.tof_tbl <- function(.data, ..., .names_sep = NULL, .key = deprecated()) {
  panel <- tof_get_panel(.data)
  #new_data <- tidyr::nest(tibble::as_tibble(.data), ...)
  #return(new_tof_tibble(x = new_data, panel = panel))
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
unnest.tof_tbl <-
  function(
    data,
    ...
  ) {
    start_panel <- tof_get_panel(data)
    return(new_tof_tibble(x = NextMethod(), panel = start_panel))
  }

#' @export
pivot_longer.tof_tbl <-
  function(
    data,
    ...
  ) {
    panel <- tof_get_panel(data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }

#' @export
pivot_wider.tof_tbl <-
  function(
    data,
    ...
  ) {
    panel <- tof_get_panel(data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }


