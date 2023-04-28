# libraries
library(ggplot2)
library(dplyr)
library(forcats)
library(testthat)


# set up data
# levine <-
#   HDCytoData::Levine_32dim_flowSet() |>
#   as_tof_tbl() |>
#   mutate(population_id = as.character(population_id)) |>
#   dplyr::rename_with(
#     .fn = \(x) stringr::str_to_lower(stringr::str_remove(x, "\\|.+"))
#   ) |>
#   dplyr::mutate(dplyr::across(c(file_number, population_id), as.character)) |>
#   tof_preprocess(
#     channel_cols = c(-time, -cell_length, -event_number, -file_number, -population_id)
#   )


# tof_assess_channels() --------------------------------------------------------

set.seed(22)
sim_data <-
  data.frame(
    cd4 = rnorm(n = 100, mean = 5, sd = 0.5),
    cd8 = rnorm(n = 100, mean = 0, sd = 0.1),
    cd33 = rnorm(n = 100, mean = 10, sd = 0.1)
  )

channel_result <-
  sim_data |>
  tof_assess_channels()

channel_result_0 <-
  sim_data |>
  # not possible
  tof_assess_channels(negative_proportion_flag = 2)

channel_result_all <-
  sim_data |>
  tof_assess_channels(negative_proportion_flag = -1)

test_that("channel evaluation has the right columns", {

  expect_true(
    all(
      colnames(channel_result) %in% c("channel", "negative_proportion", "flagged_channel")
    )
  )

  expect_equal(ncol(channel_result), 3L)
})

test_that("channel evaluation has the correct number of rows", {
  expect_equal(ncol(sim_data), nrow(channel_result))
})

test_that("correct number of channels is flagged", {
  expect_equal(sum(channel_result$flagged_channel), 1L)
  expect_equal(sum(channel_result_0$flagged_channel), 0L)
  expect_equal(sum(channel_result_all$flagged_channel), 3L)
})




# tof_calculate_flow_rate() ----------------------------------------------------

rate_data <-
  dplyr::tibble(
    time = 1:10000 %% 500,
    cd45 = rnorm(n = 10000),
    file_number = sample(c("a", "b"), size = 10000, replace = TRUE)
  )

flow_rate <-
  rate_data |>
  tof_calculate_flow_rate(time_col = time, num_timesteps = 100)

test_that("result has the right number of columns and rows", {

  expect_equal(nrow(flow_rate), 100L)
  expect_equal(ncol(flow_rate), 3L)

})

test_that("result has the right column names", {
  expect_true(
    all(colnames(flow_rate) %in% c("timestep", "time_window", "num_cells"))
  )
})

test_that("result has the right number of time windows", {
  expect_equal(100L, length(unique(flow_rate$timestep)))
  expect_equal(100L, length(unique(flow_rate$time_window)))
})

sim_rate <- data.frame(time = 1:100)

flow_rate <-
  sim_rate |>
  tof_calculate_flow_rate(time_col = time, num_timesteps = 100)

test_that("When num_timesteps is the same as the number of input rows, every count is 1", {
  expect_true(all(flow_rate$num_cells == 1))
})

test_that("When num_timesteps > the number of input rows, an error is thrown.", {
  expect_error(tof_calculate_flow_rate(sim_rate, time_col = time, num_timesteps = 200L))
})


# tof_assess_flow_rate_tibble() ------------------------------------------------

flow_rate_assessment <-
  rate_data |>
  tof_assess_flow_rate_tibble(time_col = time, augment = FALSE)

flow_rate_assessment_augment <-
  rate_data |>
  tof_assess_flow_rate_tibble(time_col = time, augment = TRUE)

test_that("result has the right number of rows and columns.", {
  expect_equal(nrow(rate_data), nrow(flow_rate_assessment))
  expect_equal(
    ncol(rate_data) + 2L,
    ncol(flow_rate_assessment_augment)
  )
  expect_equal(3L, ncol(flow_rate_assessment))
})

# tof_assess_flow_rate() -------------------------------------------------------

flow_rate_assessment <-
  rate_data |>
  tof_assess_flow_rate(
    time_col = time,
    group_cols = file_number
  )

flow_rate_assessment_2 <-
  rate_data |>
  tof_assess_flow_rate(
    time_col = time,
    alpha_threshold = 0.01
  )

flow_rate_assessment_augmented <-
  rate_data |>
  tof_assess_flow_rate(
    time_col = time,
    augment = TRUE
  )

flow_rate_assessment_augmented_2 <-
  rate_data |>
  tof_assess_flow_rate(
    time_col = time,
    group_cols = file_number,
    augment = TRUE,
  )

test_that("result has the right number of rows and columns.", {
  expect_equal(nrow(rate_data), nrow(flow_rate_assessment))
  expect_equal(4L, ncol(flow_rate_assessment))
  expect_equal(3L, ncol(flow_rate_assessment_2))
  expect_equal(ncol(rate_data) + 2L, ncol(flow_rate_assessment_augmented))
  expect_equal(ncol(rate_data) + 2L, ncol(flow_rate_assessment_augmented_2))
})



# tof_assess_clusters() --------------------------------------------------------

# tof_assess_clusters_distance()
#
# situation where a clustering algorithm splits the dataset into 3 clusters,
# but in each cluster there ends up being a somewhat more distant subcluster

sim_data_inner <-
    dplyr::tibble(
        cd45 = c(rnorm(n = 600), rnorm(n = 500, mean = -4)),
        cd38 = c(rnorm(n = 100, sd = 0.5), rnorm(n = 500, mean = -3), rnorm(n = 500, mean = 8)),
        cd34 = c(rnorm(n = 100, sd = 0.2, mean = -10), rnorm(n = 500, mean = 4), rnorm(n = 500, mean = 60)),
        cd19 = c(rnorm(n = 100, sd = 0.3, mean = 10), rnorm(n = 1000)),
        cluster_id = c(rep("a", 100), rep("b", 500), rep("c", 500)),
        dataset = "inner"
    )

sim_data_outer <-
  dplyr::tibble(
    cd45 = c(rnorm(n = 10), rnorm(50, mean = 3), rnorm(n = 50, mean = -12)),
    cd38 = c(rnorm(n = 10, sd = 0.5), rnorm(n = 50, mean = -10), rnorm(n = 50, mean = 10)),
    cd34 = c(rnorm(n = 10, sd = 0.2, mean = -15), rnorm(n = 50, mean = 15), rnorm(n = 50, mean = 70)),
    cd19 = c(rnorm(n = 10, sd = 0.3, mean = 19), rnorm(n = 100)),
    cluster_id = c(rep("a", 10), rep("b", 50), rep("c", 50)),
    dataset = "outer"
  )


sim_data <- bind_rows(sim_data_inner, sim_data_outer)

z_result <-
  sim_data |>
  tof_assess_clusters_distance(cluster_col = cluster_id, z_threshold = 2.5)

z_result_augmented <-
  sim_data |>
  tof_assess_clusters_distance(
    cluster_col = cluster_id,
    z_threshold = 2,
    augment = TRUE
  )

test_that("Z-assessment result has the right rows and columns", {
  # rows
  expect_equal(nrow(sim_data), nrow(z_result))
  expect_equal(nrow(sim_data), nrow(z_result_augmented))

  # columns
  expect_equal(ncol(z_result), 3L)
  expect_true(
    all(
      colnames(z_result) %in% c(".mahalanobis_distance", "z_score", "flagged_cell")
    )
  )
  expect_true(
    all(
      colnames(z_result_augmented) %in% c(colnames(sim_data), ".mahalanobis_distance", "z_score", "flagged_cell")
    )
  )
})

flagged_cluster_proportions <-
  z_result_augmented |>
  dplyr::count(cluster_id, flagged_cell, dataset) |>
  dplyr::group_by(cluster_id, dataset) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::filter(flagged_cell)

test_that(
  "Most cells from the Z-assessment outer dataset aka the nonbelonging dataset, are flagged;
          and most cells from the inner dataset aka the belonging dataset, are not",
  {
    expect_true(
      mean(dplyr::filter(flagged_cluster_proportions, dataset == "outer")$prop) > 0.5
    )
    expect_true(
      mean(dplyr::filter(flagged_cluster_proportions, dataset == "inner")$prop) < 0.5
    )
  })


sim_data |>
  tof_plot_cells_embedding(embedding_method = "pca", color_col = cluster_id)

sim_data |>
  tof_reduce_dimensions(method = "pca", augment = TRUE) |>
  ggplot(
    aes(x = .pc1, y = .pc2, fill = cluster_id, color = cluster_id, shape = dataset)
  ) +
  geom_point(size = 2, stroke = 0.2) +
  scale_shape_manual(values = c(21, 24)) +
  theme_bw()

z_result_augmented |>
  tof_plot_cells_embedding(
    compute_embedding_cols = starts_with("cd"),
    embedding_method = "pca",
    color_col = flagged_cell
  ) +
  ggplot2::scale_fill_viridis_c()



# tof_assess_clusters_entropy()
#
# situation where there are ambiguously identified clusters

sim_data <-
  dplyr::tibble(
    cd45 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd38 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd34 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd19 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cluster_id = c(rep("a", 1000), rep("b", 1000), rep("c", 1000))
  )

# sim_data |>
#   tof_reduce_dimensions(method = "pca") |>
#   tof_plot_cells_embedding(
#     embedding_cols = c(.pc1, .pc2),
#     color_col = cluster_id
#   )


sim_data_healthy <-
  sim_data |>
  dplyr::filter(cluster_id %in% c("b", "c"))

sim_data <-
  sim_data |>
  tof_cluster(
    healthy_tibble = sim_data_healthy,
    healthy_label_col = cluster_id,
    method = "ddpr"
  )

entropy_result <-
  sim_data |>
  tof_assess_clusters_entropy(
    cluster_col = .mahalanobis_cluster,
    marker_cols = starts_with("cd"),
    entropy_quantile = 0.8,
    augment = FALSE
  )

entropy_result_augmented <-
  sim_data |>
  tof_assess_clusters_entropy(
    cluster_col = .mahalanobis_cluster,
    marker_cols = starts_with("cd"),
    entropy_quantile = 0.8,
    augment = TRUE
  )

# entropy_assessment |>
#   ggplot(aes(x = entropy, fill = cluster_id)) +
#   geom_density(alpha = 0.4) +
#   theme_bw()

test_that("Entropy result has the right rows and columns", {
  # rows
  expect_equal(nrow(sim_data), nrow(entropy_result))
  expect_equal(nrow(sim_data), nrow(entropy_result_augmented))

  # columns
  expect_equal(ncol(entropy_result), 4L)
  expect_true(
    all(
      colnames(entropy_result) %in% c(".mahalanobis_b", ".mahalanobis_c", "entropy", "flagged_cell")
    )
  )
  expect_true(
    all(
      colnames(
        entropy_result_augmented) %in%
        c(
          colnames(sim_data),
          ".mahalanobis_b",
          ".mahalanobis_c",
          "entropy",
          "flagged_cell"
        )
    )
  )
})


flagged_cluster_proportions <-
  entropy_result_augmented |>
  group_by(cluster_id) |>
  summarize(
    prop_flagged = mean(flagged_cell)
  )

test_that(
  "Most cells from the ambiguous cluster are flagged and most cells from the other clusters are not",
  {
    expect_true(
      mean(dplyr::filter(flagged_cluster_proportions, cluster_id == "a")$prop_flagged) > 0.5
    )
    expect_true(
      mean(dplyr::filter(flagged_cluster_proportions, cluster_id != "a")$prop_flagged) < 0.5
    )
  })




# sim_data |>
#   tof_assess_clusters_distance(cluster_col = .mahalanobis_cluster, augment = TRUE, z_threshold = 1) |>
#   count(cluster_id, flagged_cell) |>
#   group_by(cluster_id) |>
#   mutate(prop = n / sum(n)) |>
#   ungroup() |> filter(flagged_cell)
#
# sim_data |>
#   tof_assess_clusters_distance(cluster_col = .mahalanobis_cluster, augment = TRUE, z_threshold = 2) |>
#   count(cluster_id, flagged_cell) |>
#   group_by(cluster_id) |>
#   mutate(prop = n / sum(n)) |>
#   ungroup() |>
#   filter(flagged_cell)

#
#
# levine <-
#   HDCytoData::Levine_32dim_flowSet() |>
#   as_tof_tbl() |>
#   mutate(population_id = as.character(population_id)) |>
#   dplyr::rename_with(
#     .fn = \(x) stringr::str_to_lower(stringr::str_remove(x, "\\|.+"))
#   ) |>
#   tof_preprocess()
#
# clusters_to_keep <-
#   levine |>
#   dplyr::count(population_id) |>
#   dplyr::slice_max(order_by = n, n = 5L) |>
#   dplyr::arrange(n) |>
#   pull(population_id)
#
# smallest_cluster <- clusters_to_keep[1]
# largest_cluster <- clusters_to_keep[[length(clusters_to_keep)]]
#
# small_levine <-
#   levine |>
#   dplyr::filter(population_id %in% clusters_to_keep)
#
# small_levine <-
#   small_levine |>
#   dplyr::mutate(
#     new_population_id =
#       if_else(
#         population_id %in% smallest_cluster,
#         sample(
#           clusters_to_keep[-which(clusters_to_keep %in% smallest_cluster)],
#           size = nrow(small_levine),
#           replace = TRUE
#         ),
#         population_id
#       )
#   )
#
# entropy_assessment <-
#   small_levine |>
#   tof_assess_clusters_entropy(
#     cluster_col = new_population_id,
#     marker_cols = starts_with("cd"),
#     num_top_clusters = 3L
#   ) |>
#   bind_cols(small_levine) |>
#   mutate(perturbed = if_else(population_id %in% smallest_cluster, "Yes", "No"))
#
# entropy_assessment_augmented <-
#   small_levine |>
#   tof_assess_clusters_entropy(
#     cluster_col = new_population_id,
#     marker_cols = starts_with("cd"),
#     num_top_clusters = 3L
#   ) |>
#   bind_cols(small_levine) |>
#   mutate(perturbed = if_else(population_id %in% smallest_cluster, "Yes", "No"))
#
#
# entropy_assessment |>
#   ggplot(aes(x = entropy, fill = population_id)) +
#   geom_density(alpha = 0.4) +
#   theme_bw()
#
#
# entropy_assessment |>
#   group_by(new_population_id) |>
#   ungroup() |>
#   mutate(population_id = fct_reorder(population_id, entropy)) |>
#   tof_plot_cells_density(
#     marker_col = entropy,
#     group_col = perturbed,
#     use_ggridges = TRUE,
#     scale = 0.25
#   )
#
# entropy_assessment |>
#   mutate(population_id = fct_reorder(population_id, entropy)) |>
#   tof_plot_cells_density(
#     marker_col = entropy,
#     group_col = population_id,
#     use_ggridges = TRUE,
#     scale = 0.25
#   )
#
# entropy_assessment |>
#   ggplot(aes(x = entropy, y = population_id, fill = perturbed)) +
#   ggridges::geom_density_ridges() +
#   theme_bw()
#
# entropy_assessment |>
#   group_by(new_population_id) |>
#   mutate(z_score = (entropy - mean(entropy)) / sd(entropy)) |>
#   ungroup() |>
#   mutate(population_id = fct_reorder(population_id, z_score)) |>
#   tof_plot_cells_density(
#     marker_col = z_score,
#     group_col = population_id,
#     use_ggridges = TRUE
#   ) +
#   ggthemes::scale_fill_tableau() +
#   labs(x = "Z score")
#
# entropy_assessment |>
#   mutate(flagged_cell = entropy > quantile(entropy, prob = 0.9)) |>
#   count(population_id, flagged_cell) |>
#   group_by(population_id) |>
#   mutate(prop = n / sum(n)) |>
#   ungroup() |>
#   filter(flagged_cell) |>
#   mutate(prop_of_flagged_cells = n / sum(n))
#
# entropy_assessment |>
#   mutate(flagged_cell = entropy > quantile(entropy, prob = 0.9))

# tof_assess_clusters_knn()
#
# situation where there are ambiguously identified clusters

# tof_assess_clusters_knn() ----------------------------------------------------

set.seed(2020L)
sim_data <-
  dplyr::tibble(
    cd45 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd38 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd34 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cd19 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
    cluster_id = c(rep("a", 1000), rep("b", 1000), rep("c", 1000))
  )

sim_data <-
  sim_data |>
  dplyr::mutate(
    new_cluster_id =
      dplyr::if_else(
        cluster_id == "a",
        sample(c("b", "c", "d"), size = 3000, replace = TRUE),
        cluster_id
      )
  )

knn_result <-
  sim_data |>
  tof_assess_clusters_knn(
    cluster_col = new_cluster_id,
    num_neighbors = 10
  )

knn_result_augmented <-
  sim_data |>
  tof_assess_clusters_knn(
    cluster_col = new_cluster_id,
    num_neighbors = 10,
    augment = TRUE
  )

flagged_cell_props <-
  knn_result_augmented |>
  dplyr::group_by(cluster_id) |>
  dplyr::summarize(prop_flagged = mean(flagged_cell)) |>
  dplyr::ungroup()

test_that("KNN result has the right rows and columns", {
  # rows
  expect_equal(nrow(sim_data), nrow(knn_result))
  expect_equal(nrow(sim_data), nrow(knn_result_augmented))

  # columns
  expect_equal(ncol(knn_result), 2L)
  expect_true(
    all(
      colnames(knn_result) %in% c(".knn_cluster", "flagged_cell")
    )
  )
  expect_true(
    all(
      colnames(
        knn_result_augmented) %in%
        c(
          colnames(sim_data),
          ".knn_cluster",
          "flagged_cell"
        )
    )
  )
})


flagged_cluster_proportions <-
  entropy_result_augmented |>
  group_by(cluster_id) |>
  summarize(
    prop_flagged = mean(flagged_cell)
  )

test_that(
  "Most cells from the nonsense cluster are flagged and most cells from the other clusters are not",
  {
    expect_true(
      mean(dplyr::filter(flagged_cell_props, cluster_id == "a")$prop_flagged) > 0.5
    )
    expect_true(
      mean(dplyr::filter(flagged_cell_props, cluster_id != "a")$prop_flagged) < 0.5
    )
  })

