
######################   tof_preprocess

#' Preprocess raw CyTOF data.
#'
#' tof_preprocess() transforms a tibble of raw ion counts directly measured on a mass cytometer.
#'
#' @param tof_tibble a tibble, data.frame, or something that can be coerced into either
#' @param metadata_vars variables that contain information about variables that should not
#' computed over, i.e. file names, patient names, stimulation names, etc.
#' Supports tidy selection using tidy_select helpers.
#' Not currently used.
#' @param channel_vars A vector of non-quoted variables representing columns that contain
#' single-cell protein measurements. Anything that works in the first argument of
#' dplyr::across will work. See ?across. Supports tidy selection using tidy_select
#' helpers. If nothing is specified, the default is to transform all numeric columns.
#' @param undo_noise logical indicating if you'd like to remove the uniform noise that
#' Fluidigm software adds to each protein measurement for aesthetic
#' and visualization purposes. See _____. Default = TRUE.
#' @param transform_fun vectorized function to apply to each protein value for
#' variance stabilization. Default is asinh transformation (with a co-factor of 5).
#'
#' @return a tibble with identical dimensions to the input tof_tibble, with all
#' columns specified in channel_vars transformed using transform_fun (with noise
#' removed or not removed depending on undo_noise).
#'
#' @export
#'
#' @examples
tof_preprocess <-

  function(
    tof_tibble = NULL,
    metadata_vars = NULL,
    channel_vars = where(is.numeric),
    undo_noise = TRUE,
    transform_fun = function(x) asinh(x/5)
  ) {
    #channel_vars <- enexprs(channel_vars)
    if (undo_noise) {
      tof_tibble <-
        tof_tibble %>%
        mutate(across({{channel_vars}}, ~ floor(.x) + 1))
    }
    tof_tibble <-
      tof_tibble %>%
      mutate(across({{channel_vars}}, transform_fun))

    return(tof_tibble)
  }


