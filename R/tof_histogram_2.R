################   tof_histogram_2.R 
### Description: 
### Plots a series of overlain histograms for a given CyTOF channel broken up 
### by a given grouping variable and a given splitting variable.
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#       of these
#     - channel_var = an unquoted variable name indicating which channel should be plotted 
#     - group_var = an unquoted variable name indicating which variable should be used to 
#                   break the cells into different histograms along the overlay.
#     - split_var = an unquoted variable name indicating which variable should be used to break 
#                   the cells into different histograms along the axis of `channel_var`.
#     - ordered = a logical indicating the the histograms along the y axis should be ordered by decreasing 
#                 median channel value.
#     - color_option = color palette name to be passed to scale_color_viridis as argument `option`. 
#                      Choices include "inferno", "plasma", "viridis", "magma", and "cividis". 
#     - lower_quantile = lowest quantile in `channel_var` to include in the density plot. ALl cells with 
#                        `channel_var` values lower than this will be filtered out of the plot. 
#     - upper_quantile = highest quantile in `channel_var` to include in the density plot. ALl cells with 
#                        `channel_var` values higher than this will be filtered out of the plot. 
#     - ... = additional arguments to be passed to `geom_density_ridges()`. 
#
# Outputs: 
#     - a ggplot object 
#
# Dependencies: 
#     - tidyverse
#     - rlang
#     - ggridges

tof_histogram_2 <- 
  function(
    tof_tibble = NULL,   
    channel_var = NULL, 
    group_var = NULL,
    split_var = NULL, 
    ordered = FALSE, 
    lower_quantile = 0.01,
    upper_quantile = 0.99, 
    ... # to be passed to `geom_density_ridges()`
  ) {
    if (ordered) { 
      tof_tibble <-
        tof_tibble %>% 
        mutate(
          "{{group_var}}" := 
            fct_reorder({{group_var}}, {{channel_var}}, median, .desc = TRUE)
        )
    }
    
    tof_tibble <- 
      tof_tibble %>% 
      select({{channel_var}}, {{group_var}}, {{split_var}}) %>% 
      filter(
        {{channel_var}} <= quantile({{channel_var}}, upper_quantile),
        {{channel_var}} >= quantile({{channel_var}}, lower_quantile),
      )
    
    tof_tibble %>% 
      ggplot(aes(x = {{channel_var}}, y = {{group_var}}, fill = {{split_var}})) + 
      geom_density_ridges(...) + 
      scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) + 
      scale_fill_tableau()
  }
