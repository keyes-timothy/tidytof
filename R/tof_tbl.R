
# new_tof_tibble ---------------------------------------------------------------

#' Constructor for a tof_tibble.
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
#' @family tof_tbl utilities
#'
new_tof_tibble <- function(x = dplyr::tibble(), panel = dplyr::tibble()) {

  stopifnot(inherits(x, "tbl_df"))
  stopifnot(inherits(panel, "tbl_df"))

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

# tof_get_panel ----------------------------------------------------------------

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
#' @family tof_tbl utilities
#'
#'
tof_get_panel <- function(tof_tibble) {
  panel <-
    tof_tibble %>%
    attr(which = "panel")

  return(panel)
}


# tof_set_panel ----------------------------------------------------------------

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
#' @family tof_tbl utilities
#'
#' @export
#'
#'
tof_set_panel <- function(tof_tibble, panel) {
  attr(tof_tibble, which = "panel") <- panel
  return(tof_tibble)
}


# tof_tbl methods --------------------------------------------------------------

## tidyr methods

#' @export
nest.tof_tbl <- function(.data, ..., .names_sep = NULL) {#, .key = deprecated()) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

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

#' @export
#'
#' @importFrom dplyr group_by_drop_default
group_by.tof_tbl <-
  function(.data, ..., .add = FALSE, .drop = dplyr::group_by_drop_default(.data)) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }


## dplyr methods

#' @export
mutate.tof_tbl <- function(.data, ...) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
slice_sample.tof_tbl <-
  function(.data, ..., n, prop, weight_by = NULL, replace = FALSE) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }

# grouped_tof_tbl methods ------------------------------------------------------

## tidyr methods

#' @export
nest.grouped_tof_tbl <- nest.tof_tbl

#' @export
ungroup.grouped_tof_tbl <- function(x, ...) {
  panel <- tof_get_panel(x)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

# for interoperability with flowCore -------------------------------------------

#' @export
#'
#' @importFrom flowCore fsApply
#' @importFrom flowCore exprs
#'
#' @importFrom dplyr as_tibble
as_tof_tbl.flowSet <- function(flow_data, sep = "|") {
  # check if flowset is empty
  if (length(flow_data) < 1) {
    stop("This flowSet is empty.")
  }
  panel_info <-
    flow_data[[1]] %>%
    tof_find_panel_info()

  flowset_exprs <-
    flow_data %>%
    flowCore::fsApply(FUN = flowCore::exprs) %>%
    dplyr::as_tibble()

  col_names <-
    base::paste(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  colnames(flowset_exprs) <- col_names

  result <- new_tof_tibble(x = flowset_exprs, panel = panel_info)

  return(result)
}

#' @export
#'
#' @importFrom dplyr as_tibble
#' @importFrom flowCore exprs
#'
as_tof_tbl.flowFrame <- function(flow_data, sep = "|") {
  panel_info <-
    flow_data %>%
    tof_find_panel_info()

  col_names <-
    #stringr::str_c(panel_info$antigens, panel_info$metals, sep = sep)
    base::paste(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  flowframe_exprs <-
    setNames(
      object = dplyr::as_tibble(flowCore::exprs(flow_data)),
      nm = col_names
    )

  result <-
    new_tof_tibble(
      x = flowframe_exprs,
      panel = panel_info
    )

  return(result)

}

#' Coerce flowFrames or flowSets into tof_tbl's.
#'
#' @param flow_data A flowFrame or flowSet
#'
#' @param sep A string indicating which symbol should be used to separate
#' antigen names and metal names in the columns of the output tof_tbl.
#'
#'
#' @export
as_tof_tbl <- function(flow_data, sep = "|") {
  UseMethod("as_tof_tbl")
}

