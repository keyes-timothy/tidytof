# library(dplyr)
# library(tidyr)
# library(tidytof)
# library(scales)
# library(magrittr)
# library(ggplot2)
#
# sample_list <-
#   map(
#     .x = seq(1, 20, length.out = 100),
#     .f = \(x)
#     rnorm(n = 40, mean = x, sd = 1.5) + rnorm(n = 40, mean = 0, sd = 0.3) |>
#       as.matrix(nrow = 40)
#   )
#
# samples <- purrr::reduce(.x = sample_list, .f = cbind)
#
# temp_tibble <-
#   dplyr::tibble(
#     channel_1 = c(1, 2, 3, 4, 5, 6),
#     group = c("a", "a", "a", "b", "b", "b")
#   )
#
# sample_tibble <- dplyr::as_tibble(t(samples))
# colnames(sample_tibble) <- paste0("channel_", 1:ncol(sample_tibble))
# sample_tibble$random <-
#   sample(x = c("a", "b"), size = nrow(sample_tibble), replace = TRUE)
#
#
# # apply quantile normalization without groups
# temp <-
#   sample_tibble |>
#   tof_batch_correct_quantile(channel_cols = -random)
#
# sample_tibble |>
#   select(-random) |>
#   t() |>
#   as_tibble() |>
#   tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = values, group = sample)) +
#   geom_density(linewidth = 0.1) +
#   labs(subtitle = "pre")
#
# temp |>
#   select(-random) |>
#   t() |>
#   as_tibble() |>
#   tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = values, group = sample)) +
#   geom_density() +
#   scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
#   labs(subtitle = "post")
#
# # apply quantile normalization with blocks
#
# temp_2 <-
#   sample_tibble |>
#   tof_batch_correct_quantile(
#     channel_cols = starts_with("channel_"),
#     group_cols = random
#   )
#
# sample_tibble |>
#   select(-random) |>
#   t() |>
#   as_tibble() |>
#   tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = values, group = sample)) +
#   geom_density(linewidth = 0.1) +
#   labs(subtitle = "pre")
#
# temp_2 |>
#   select(-random) |>
#   t() |>
#   as_tibble() |>
#   tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = values, group = sample)) +
#   geom_density() +
#   scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
#   labs(subtitle = "post")
#
#
# # rescaling
#
# temp_3 <-
#   sample_tibble |>
#   tof_batch_correct_rescale(
#     channel_cols = starts_with("channel_")
#   )
#
# sample_tibble |>
#   # select(-random) |>
#   # t() |>
#   # as_tibble() |>
#   # tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = channel_1, color = random)) +
#   geom_density() +
#   labs(subtitle = "pre")
#
# temp_3 |>
#   #select(-random) |>
#   #t() |>
#   #as_tibble() |>
#   #tidyr::pivot_longer(cols = everything(), names_to = "sample", values_to = "values") |>
#   ggplot2::ggplot(ggplot2::aes(x = channel_1, color = random)) +
#   geom_density() +
#   scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
#   labs(subtitle = "post")
#
# temp_4 <-
#   sample_tibble |>
#   tof_batch_correct_rescale(
#     channel_cols = starts_with("channel_"),
#     group_cols = random
#   )
#
# sample_tibble |>
#   ggplot2::ggplot(ggplot2::aes(x = channel_1, color = random)) +
#   geom_density() +
#   labs(subtitle = "pre")
#
# temp_4 |>
#   ggplot2::ggplot(ggplot2::aes(x = channel_1, color = random)) +
#   geom_density() +
#   scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
#   labs(subtitle = "post")
#
# plot_fun <- function(tof_tibble, subtitle = "") {
#   result <-
#     tof_tibble |>
#     ggplot(aes(x = channel_1, color = random)) +
#     geom_density() +
#     scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
#     labs(subtitle = subtitle)
#
#   print(result)
# }
#
# temp_5 <-
#   sample_tibble %T>%
#   plot_fun(subtitle = "raw") |>
#   tof_batch_correct_rescale(channel_cols = channel_1, group_cols = random) %T>%
#   plot_fun(subtitle = "separate") |>
#   tof_batch_correct_rescale(channel_col = channel_1) |>
#   plot_fun(subtitle = "integrated")
#
#
