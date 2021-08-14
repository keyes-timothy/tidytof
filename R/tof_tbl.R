
# new_tof_tibble ------------------

#' Constructor for tof_tibble.
#'
#' @param x A data.frame or tibble containing single-cell mass cytometry data
#' such that rows are cells and columns are CyTOF measurements.
#'
#' @param panel A data.frame or tibble containing information about the panel
#' for the mass cytometry data in x.
#'
#' @return A `tof_tbl`, an tibble extension that tracks a few other attributes
#' that are useful for CyTOF data analysis.
#'
#' @examples
#'
new_tof_tibble <- function(x = tibble::tibble(), panel = tibble::tibble()) {

  stopifnot(tibble::is_tibble(x))
  stopifnot(tibble::is_tibble(panel))

  if("grouped_df" %in% class(x)) {
    subclasses <- c("grouped_tof_tbl", "grouped_df", "tof_tbl")
  } else {
    subclasses <- "tof_tbl"
  }

  tibble::new_tibble(
    x,
    panel = panel,
    nrow = nrow(x),
    class = subclasses
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


# tof_set_panel ---------------------------

#' Set panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @param panel A tibble containing two columns (`metals` and `antigens`) representing
#' the information about a panel
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

# tidyr methods

#' @export
nest.tof_tbl <- function(.data, ..., .names_sep = NULL, .key = deprecated()) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

# needs to be optimized
#' @export
unnest.tof_tbl <- function(data, ...) {
  start_panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = start_panel))
}


#' @export
pivot_longer.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
pivot_wider.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @importFrom dplyr group_by_drop_default
#' @export
group_by.tof_tbl <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}



# dplyr methods

#' @export
mutate.tof_tbl <- function(.data, ...) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}



# grouped_tof_tbl methods ---------------------------------

# tidyr methods

#' @export
nest.grouped_tof_tbl <- nest.tof_tbl

#' @export
ungroup.grouped_tof_tbl <- function(x, ...) {
  panel <- tof_get_panel(x)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}
