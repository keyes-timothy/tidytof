# libraries
library(tidytof)
library(dplyr)
library(testthat)

# set up
reference_tibble <-
  ddpr_data %>%
  slice_sample(prop = 0.1) %>%
  mutate(my_cluster = sample(x = c('a', 'b', 'c'), size = n(), replace = TRUE))

# tof_upsample_distance --------------------------------------------------------

test_that("upsample result is a tibble with a single character vector column of correct length", {

  result <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      return_distances = FALSE
    )

  expect_equal(nrow(ddpr_data), nrow(result))
  expect_true(is.character(result$.upsample_cluster))
  expect_equal(ncol(result), 1L)

})

test_that("return_distances argument works", {
  result <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      return_distances = TRUE
    )

  expect_equal(ncol(result), 4L)

})

test_that("output columns are named correctly", {
  result_mah <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      return_distances = TRUE
    )

  result_cosine <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      distance_function = "cosine",
      return_distances = TRUE,
      parallel_cols = sample_name
    )

  result_pearson <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      distance_function = "pearson",
      return_distances = TRUE
    )

  # all column names should start with "."
  expect_true(all(grepl("^\\.", colnames(result_mah))))
  expect_true(all(grepl("^\\.", colnames(result_cosine))))
  expect_true(all(grepl("^\\.", colnames(result_pearson))))

  # column names for distances should start with ".{distance_function}"
  expect_equal(sum(grepl("\\.mahalanobis", colnames(result_mah))), 3L)
  expect_equal(sum(grepl("\\.cosine", colnames(result_cosine))), 3L)
  expect_equal(sum(grepl("\\.pearson", colnames(result_pearson))), 3L)
})

# tof_upsample_neighbor --------------------------------------------------------

test_that("upsample result is a tibble with a single character vector column of correct length", {

  result <-
    ddpr_data %>%
    tof_upsample_neighbor(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19)
    )

  expect_equal(nrow(ddpr_data), nrow(result))
  expect_true(is.character(result$.upsample_cluster))
  expect_equal(ncol(result), 1L)

})

test_that("output columns are named correctly", {
  result <-
    ddpr_data %>%
    tof_upsample_neighbor(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19)
    )

  expect_identical(colnames(result), ".upsample_cluster")
})

# tof_upsample -----------------------------------------------------------------

test_that("tof_upsample and tof_upsample_distance results are the same", {
  result_1 <-
    ddpr_data %>%
    tof_upsample_distance(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      return_distances = TRUE
    )

  result_2 <-
    ddpr_data %>%
    tof_upsample(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      return_distances = TRUE,
      augment = FALSE,
      method = "distance"
    )

  expect_equal(result_1, result_2)

})

test_that("tof_upsample and tof_upsample_neighbor results are the same", {
  result_1 <-
    ddpr_data %>%
    tof_upsample_neighbor(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
    )

  result_2 <-
    ddpr_data %>%
    tof_upsample(
      reference_tibble = reference_tibble,
      reference_cluster_col = my_cluster,
      upsample_cols = c(cd34, cd45, cd19),
      augment = FALSE,
      method = "neighbor"
    )

  expect_equal(result_1, result_2)

})
