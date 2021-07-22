# preprocessing.R
# This file contains functions relevant to performing several common preprocessing steps
# in CyTOF data analysis, including hyperbolic arcsine transformation.


# tof_preprocess ----------------------------

#' Preprocess raw CyTOF data.
#'
#' This function transforms a `tof_tibble` of raw ion counts directly measured on
#' a mass cytometer using a user-provided function.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble`.
#'
#' @param metadata_vars Unquoted column names of columns in `tof_tibble` that should
#' not be computed over, i.e. file names, patient names, stimulation names, etc.
#' Supports tidy selection using tidy_select helpers.
#' Not currently used.
#'
#' @param channel_vars A vector of non-quoted variables representing columns that contain
#' single-cell protein measurements. Anything that works in the first argument of
#' dplyr::across will work. See ?across. Supports tidy selection using tidyselect
#' helpers. If nothing is specified, the default is to transform all numeric columns.
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
#' @return A `tof_tibble` with identical dimensions to the input `tof_tibble`, with all
#' columns specified in channel_vars transformed using `transform_fun` (with noise
#' removed or not removed depending on `undo_noise`).
#'
#' @export
#'
#' @examples
#' NULL
#'
#' Note: This can be optimized by keeping track of all tranformations/analyses
#' that have been applied to a tof_tibble (giving it memory)...then creating a
#' function that can automatically undo it? Or it can at least be referenced.
#'
#' Note: Convert this to a dtplyr backend and so some speed testing?
#'
#'
tof_preprocess <-

  function(
    tof_tibble = NULL,
    metadata_vars = NULL,
    channel_vars = where(is.numeric),
    undo_noise = TRUE,
    transform_fun = function(x) asinh(x/5)
  ) {

    # first remove noise if specified
    if (undo_noise) {
      tof_tibble <-
        tof_tibble %>%
        mutate(across({{channel_vars}}, ~ floor(.x) + 1))
    }

    # then apply transformation to all channel_vars
    tof_tibble <-
      tof_tibble %>%
      mutate(across({{channel_vars}}, transform_fun))

    return(tof_tibble)
  }


#' Preprocess transformed CyTOF data.
#'
#' This function transforms a `tof_tibble` of transformed ion counts from a mass
#' cytometer back into something that looks more like an .fcs file that Fluidigm
#' software generates.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble`.
#'
#' @param metadata_vars Unquoted column names of columns in `tof_tibble` that should
#' not be computed over, i.e. file names, patient names, stimulation names, etc.
#' Supports tidy selection using tidy_select helpers.
#' Not currently used.
#'
#' @param channel_vars A vector of non-quoted variables representing columns that contain
#' single-cell protein measurements. Anything that works in the first argument of
#' dplyr::across will work. See ?across. Supports tidy selection using tidyselect
#' helpers. If nothing is specified, the default is to transform all numeric columns.
#'
#' @param redo_noise A boolean value indicating whether to add  uniform noise that
#' to each CyTOF measurement for aesthetic and visualization purposes. See \href{https://pubmed.ncbi.nlm.nih.gov/30277658/}{this paper}.
#' Defaults to FALSE
#'
#' @param transform_fun A vectorized function to apply to each column specified by
#' `channel_vars` for post-processing. Defaults to \code{\link{rev_asinh}} transformation
#' (with a cofactor of 5).
#'
#' @return A `tof_tibble` with identical dimensions to the input `tof_tibble`, with all
#' columns specified in channel_vars transformed using `transform_fun` (with noise
#' added or not removed depending on `redo_noise`).
#'
#' @export
#'
#' @examples
tof_postprocess <-
  function(
    tof_tibble = NULL,
    metadata_vars = NULL,
    channel_vars = where(is.numeric),
    redo_noise = FALSE,
    transform_fun = function(x) rev_asinh(x, shift_factor = 0, scale_factor = 0.2)
  ) {

    # first apply transformation function to all channel_vars
    tof_tibble <-
      tof_tibble %>%
      mutate(across({{channel_vars}}, transform_fun))

    # then remove noise if specified
    if (redo_noise) {
      tof_tibble <-
        tof_tibble %>%
        mutate(across({{channel_vars}}, ~ .x - stats::runif(n = length(.x))))
    }

    return(tof_tibble)
  }





