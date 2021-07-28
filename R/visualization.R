######################   tof_histogram

#' Generate univariate histograms of CyTOF data from a `tof_tibble`.
#'
#' Plots a series of overlaid histograms for a given CyTOF channel broken up
#' by a given grouping variable.
#'
#' @param tof_tibble A tibble, data.frame, or something that can be coerced into either
#'
#' @param channel_var An unquoted variable name indicating which channel should be plotted
#'
#' @param group_var An unquoted variable name indicating which variable should be used to
#' break the cells into different histograms along the overlay.
#'
#' @param ordered A logical indicating the the histograms along the y axis should
#' be ordered by decreasing median channel value.
#'
#' @param color_option Color palette name to be passed to scale_color_viridis as argument `option`.
#' Choices include "inferno", "plasma", "viridis", "magma", and "cividis".
#'
#' @param lower_quantile lowest quantile in `channel_var` to include in the density plot. All cells with
#' `channel_var` values lower than this will be filtered out of the plot.
#'
#' @param upper_quantile highest quantile in `channel_var` to include in the density plot. All cells with
#' `channel_var` values higher than this will be filtered out of the plot.
#'
#' @param ... Optional additional arguments to be passed to `geom_tof_ridges()`
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_histogram <- function(
  tof_tibble = NULL,
  channel_var = NULL,
  group_var = NULL,
  ordered = FALSE,
  color_option = "plasma",
  lower_quantile = 0.01,
  upper_quantile = 0.99,
  ...
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
    select({{channel_var}}, {{group_var}}) %>%
    filter(
      {{channel_var}} <= quantile({{channel_var}}, upper_quantile),
      {{channel_var}} >= quantile({{channel_var}}, lower_quantile),
    )

  tof_tibble %>%
    ggplot(aes(x = {{channel_var}}, y = {{group_var}}, fill = stat(x))) +
    geom_vline(
      xintercept = mean(pull(tof_tibble, {{channel_var}})),
      linetype = "dotted",
      size = 0.5
    ) +
    geom_tof_ridges(...) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_viridis_c(option = color_option) +
    theme_ridges() +
    theme(legend.position = "bottom") +
    guides(
      fill =
        guide_colorbar(
          title.position = "top",
          title.theme = element_text(size = 12),
          title.hjust = 0.5,
          barwidth = 10,
          barheight = 0.5
        )
    ) +
    labs(
      subtitle = NULL,
      x = NULL,
      y = as_name(enexpr(group_var)),
      fill = str_c(as_name(enexpr(channel_var)), " expression")
    )
}





#' Generate grouped univariate histograms of CyTOF data from a `tof_tibble`.
#'
#' @param tof_tibble A tibble, data.frame, or something that can be coerced into either
#'
#' @param channel_var An unquoted variable name indicating which channel should be plotted
#'
#' @param group_var An unquoted variable name indicating which variable should be used to
#' break the cells into different histograms along the overlay.
#'
#' @param split_var an unquoted variable name indicating which variable should be used to break
#' the cells into different histograms along the axis of `channel_var`.
#'
#' @param ordered A logical indicating the the histograms along the y axis should
#' be ordered by decreasing median channel value.
#'
#' @param lower_quantile lowest quantile in `channel_var` to include in the density plot. All cells with
#' `channel_var` values lower than this will be filtered out of the plot.
#'
#' @param upper_quantile highest quantile in `channel_var` to include in the density plot. All cells with
#' `channel_var` values higher than this will be filtered out of the plot.
#'
#' @param ... Optional additional arguments to be passed to `geom_tof_ridges()`
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_histogram_2 <-
  function(
    tof_tibble = NULL,
    channel_var = NULL,
    group_var = NULL,
    split_var = NULL,
    ordered = FALSE,
    lower_quantile = 0.01,
    upper_quantile = 0.99,
    ...
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




#' Title
#'
#' Plots a series of overlain histograms for all CyTOF channels in your dataset broken up
#' by a given grouping variable.
#'
#' @param tof_tibble a tibble, data.frame, or somehing that can be coerced into either
#' of these
#'
#' @param group_var an unquoted variable name indicating which variable should be used to
#' break the cells into different histograms.
#'
#' @param out_path file path to the directory in which the plots should be saved
#'
#' @param label_prefix An optional string to concatenate to the beginning of each file name.
#'
#' @param label_suffix An optional string to concatenate to the beginning of each file name.
#'
#' @param width The width for the plots being saved.
#'
#' @param height The height of the plots being saved.
#'
#' @param device "jpg", "tif", or "pdf" file format
#'
#' @param ... Optional additional arguments to pass to `tof_histogram()`
#'
#' @return Save all density plots to the file path specified by `out_path`.
#'
#' @export
#'
#' @examples
#' NULL
#'
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

