# preprocessing.R
# This file contains functions relevant to performing several common preprocessing
# (and postprocessing) steps in CyTOF data analysis, including hyperbolic
# arcsine transformation.


# tof_preprocess ---------------------------------------------------------------

#' Preprocess raw CyTOF data.
#'
#' This function transforms a `tof_tbl` of raw ion counts directly measured on
#' a mass cytometer using a user-provided function. It can be used to perform
#' standard pre-processing steps before CyTOF data analysis.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param channel_cols A vector of non-quoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#' If nothing is specified, the default is to transform all numeric columns.
#'
#' @param undo_noise A boolean value indicating whether to remove the uniform noise that
#' Fluidigm software adds to each CyTOF measurement for aesthetic
#' and visualization purposes. See \href{https://pubmed.ncbi.nlm.nih.gov/30277658/}{this paper}.
#' Defaults to TRUE.
#'
#' @param transform_fun A vectorized function to apply to each protein value for
#' variance stabilization. Defaults to \code{\link[base]{asinh}} transformation
#' (with a co-factor of 5).
#'
#' @return A `tof_tbl` with identical dimensions to the input `tof_tibble`, with all
#' columns specified in channel_cols transformed using `transform_fun` (with noise
#' removed or not removed depending on `undo_noise`).
#'
#' @seealso [tof_postprocess()]
#'
#' @export
#'
#'
#'
tof_preprocess <-
  function(
    tof_tibble = NULL,
    channel_cols = where(tof_is_numeric),
    undo_noise = TRUE,
    transform_fun = function(x) asinh(x/5)
  ) {
    # first remove noise if specified
    if (undo_noise) {
      tof_tibble <-
        tof_tibble %>%
        dplyr::mutate(across({{channel_cols}}, ~ floor(.x) + 1))
    }

    # then apply transformation to all channel_cols
    tof_tibble <-
      tof_tibble %>%
      dplyr::mutate(across({{channel_cols}}, transform_fun))

    return(tof_tibble)
  }



# tof_postprocess --------------------------------------------------------------

#' Post-process transformed CyTOF data.
#'
#' This function transforms a `tof_tibble` of transformed ion counts from a mass
#' cytometer back into something that looks more like an .fcs file that Fluidigm
#' software generates.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble`.
#'
#' @param channel_cols A vector of non-quoted column names indicating which columns
#' in `tof_tibble` contain protein measurements. Supports tidyselect helpers.
#' If nothing is specified, the default is to transform all numeric columns.
#'
#' @param redo_noise A boolean value indicating whether to add  uniform noise that
#' to each CyTOF measurement for aesthetic and visualization purposes. See \href{https://pubmed.ncbi.nlm.nih.gov/30277658/}{this paper}.
#' Defaults to FALSE
#'
#' @param transform_fun A vectorized function to apply to each column specified by
#' `channel_cols` for post-processing. Defaults to \code{\link{rev_asinh}} transformation
#' (with a cofactor of 5).
#'
#' @return A `tof_tbl` with identical dimensions to the input `tof_tibble`, with all
#' columns specified in channel_cols transformed using `transform_fun` (with noise
#' added or not removed depending on `redo_noise`).
#'
#' @seealso [tof_preprocess()]
#'
#' @export
#'
#'
tof_postprocess <-
  function(
    tof_tibble = NULL,
    channel_cols = where(tof_is_numeric),
    redo_noise = FALSE,
    transform_fun = function(x) rev_asinh(x, shift_factor = 0, scale_factor = 0.2)
  ) {

    # first apply transformation function to all channel_cols
    tof_tibble <-
      tof_tibble %>%
      dplyr::mutate(across({{channel_cols}}, transform_fun))

    # then remove noise if specified
    if (redo_noise) {
      tof_tibble <-
        tof_tibble %>%
        dplyr::mutate(across({{channel_cols}}, ~ .x - stats::runif(n = length(.x))))
    }

    return(tof_tibble)
  }





