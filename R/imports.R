#' @import dplyr
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
NULL

# dplyr reexports --------------------------------------------------------------

#' @importFrom dplyr %>%
#' @export
dplyr::`%>%`


# rlang reexports --------------------------------------------------------------

#' @importFrom rlang `:=`
#' @export
rlang::`:=`

#' @importFrom rlang `:=`
#' @export
rlang::.data

# tidyselect reexports ---------------------------------------------------------

#' @importFrom tidyselect vars_select_helpers
#' @export
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


