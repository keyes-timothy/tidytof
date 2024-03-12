#' @import Rcpp
#'
#' @importFrom methods as
#'
#' @importFrom stats cov
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats setNames
#'
#' @importFrom utils data
#' @importFrom utils file_test
#'
#' @importFrom rlang are_na
#'
#' @importFrom foreach `%dopar%`
#' @importFrom foreach `%do%`
#'
NULL

# dplyr reexports --------------------------------------------------------------

#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`


# rlang reexports --------------------------------------------------------------

#' @importFrom rlang `:=`
#' @export
rlang::`:=`

#' @importFrom rlang `.data`
#' @export
rlang::`.data`

# tidyselect reexports ---------------------------------------------------------

#' Select variables with a function
#'
#' This is a copy of \code{\link[tidyselect]{where}}, a selection helper that
#' selects the variables for which a predicate function returns TRUE. See
#'  \code{\link[tidyselect]{language}} for more details about tidyselection.
#'
#' This help file was replicated verbatim from \code{\link[tidyselect]{tidyselect-package}}.
#'
#' @param fn A function that returns TRUE or FALSE (technically, a predicate function).
#' Can also be a purrr-like formula.
#'
#' @importFrom tidyselect vars_select_helpers
#'
#' @export
#'
#' @return A predicate that can be used to select columns from a data.frame.
#'
#' @references  Lionel Henry and Hadley Wickham (2021). tidyselect:
#' Select from a Set of Strings. R package version 1.1.1.
#' https://CRAN.R-project.org/package=tidyselect
#'
#' @examples
#' NULL
#'
where <- tidyselect::vars_select_helpers$where

# Alias required for help links in downstream packages

#' @aliases select_helpers
#' @importFrom tidyselect contains
#' @export
tidyselect::contains

#' @importFrom tidyselect ends_with
#' @export
tidyselect::ends_with

#' @importFrom tidyselect everything
#' @export
tidyselect::everything

#' @importFrom tidyselect matches
#' @export
tidyselect::matches

#' @importFrom tidyselect num_range
#' @export
tidyselect::num_range

#' @importFrom tidyselect starts_with
#' @export
tidyselect::starts_with

#' @importFrom tidyselect last_col
#' @export
tidyselect::last_col

#' @importFrom tidyselect any_of
#' @export
tidyselect::any_of

#' @importFrom tidyselect all_of
#' @export
tidyselect::all_of
