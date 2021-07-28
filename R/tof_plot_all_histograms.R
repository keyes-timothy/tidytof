################   tof_plot_all_histograms.R
### Description:
### Plots a series of overlain histograms for all CyTOF channels in your dataset broken up
### by a given grouping variable.
#
# Inputs:
#     - tof_tibble = a tibble, data.frame, or somehing that can be coerced into either
#       of these
#     - group_var = an unquoted variable name indicating which variable should be used to
#                   break the cells into different histograms.
#     - out_path = file path to the directory in which the plots should be saved
#     - label_column = an unquoted column name representing a column in `plot_tibble`
#                      containing strings that can be used to uniquely identify each
#                      plot in `plot_tibble`.
#     - plot_column = an unquoted column name indicating which column in `plot_tibble`
#                      contains the plots to be saved.
#     - label_prefix = an optional string to concatenate to the beginning of each file name.
#     - label_suffix = an optional string to concatenate to the beginning of each file name.
#     - width = the width for the plots being saved.
#     - height = the height of the plots being saved.
#     - device = "jpg", "tif", or "pdf" file format.
#     - ... = additional arguments to pass to `tof_histogram()`
#
# Outputs:
#     - None
#
# Side-effects:
#     - Save all density plots to the file path specified by `out_path`.
#
# Dependencies:
#     - tidyverse
#     - rlang

tof_plot_all_histograms <-
  function(
    tof_tibble,
    group_var,
    out_path = NULL,
    label_prefix = "",
    label_suffix = "",
    width,
    height,
    device = "pdf",
    ...
  ) {
    my_channels <- tof_tibble %>%
      select_if(is.numeric) %>%
      colnames()

      my_channels %>%
      syms() %>%
      walk(
        .f =
          function(x) {
            tof_histogram(
              tof_tibble = tof_tibble,
              channel_var = !!x,
              group_var = {{group_var}},
              ...
            ) %>%
              ggsave(
                filename = str_c(label_prefix, x, label_suffix, sep = "-"),
                plot = .,
                device = device,
                path = out_path,
                width = width,
                height = height,
                units = "in"
              )
          }
      )
  }
