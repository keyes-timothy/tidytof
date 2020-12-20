################   tof_save_figures.R 
### Description: 
### Saves a series of plots contained in a plot_tibble. 
#
# Inputs: 
#     - plot_tibble = a tibble containing a column of ggplot objects
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
#     - ... = additional arguments to pass to `ggsave()`
#
# Outputs: None
#
# Side-effects:
#     - Saves each plot in the column called `plot_column` in `plot_tibble` 
#       using `ggsave()`.
#
# Dependencies: 
#     - tidyverse
#     - rlang


tof_save_figures <- 
  function(
    plot_tibble = NULL,
    out_path = NULL, 
    label_column = NULL,
    plot_column = plots,
    label_prefix = "",
    label_suffix = "", 
    width, 
    height,
    device = "pdf", 
    ...
  ) { 
    walk2(
      .x = pull(plot_tibble, {{plot_column}}), 
      .y = pull(plot_tibble, {{label_column}}), 
      .f = 
        ~ggsave(
          filename = str_c(label_prefix, .y, label_suffix, sep = "-"),
          plot = .x, 
          device = device,
          path = out_path,
          width = width, 
          height = height, 
          units = "in"
        )
    )
  }
