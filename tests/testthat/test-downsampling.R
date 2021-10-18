library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidytof)

# tof_downsample_constant ------------------------------------------------------

### downsampling without grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_constant(num_cells = 100L)

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))

  # number of rows is correct
  expect_identical(nrow(result), 100L)
  })

### downsampling with grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_constant(group_cols = sample_name, num_cells = 100L)

  result_2 <-
    phenograph_data %>%
    tof_downsample_constant(
      group_cols = c(sample_name, phenograph_cluster),
      num_cells = 100L
    )

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))
  expect_identical(ncol(result_2), ncol(phenograph_data))


  # number of rows is correct
  expect_identical(nrow(result), 300L)
  expect_identical(nrow(result), 300L)
})


# tof_downsample_prop ----------------------------------------------------------

### downsampling without grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_prop(prop_cells = 0.1)

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))

  # number of rows is correct
  expect_equal(nrow(result), nrow(phenograph_data) * 0.1)
})

### downsampling with grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_prop(group_cols = sample_name, prop_cells = 0.1)

  result_2 <-
    phenograph_data %>%
    tof_downsample_prop(
      group_cols = c(sample_name, phenograph_cluster),
      prop_cells = 0.1
    )

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))
  expect_identical(ncol(result_2), ncol(phenograph_data))

  # number of rows is correct
  expect_equal(nrow(result), 0.1 * nrow(phenograph_data))
  expect_equal(nrow(result), 0.1 * nrow(phenograph_data))
})


# tof_downsample_density ----------------------------------------------------------



### downsampling without grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_density(target_percentile = 0.05, density_cols = c(cd19, cd11b))

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))

  # number of rows is smaller after density-dependent downsampling
  expect_true(nrow(result) < nrow(phenograph_data))
})

### downsampling with grouping

test_that("Output shape is correct.", {
  result <-
    phenograph_data %>%
    tof_downsample_density(
      group_cols = sample_name,
      target_percentile = 0.05,
      density_cols = c(cd19, cd11b)
    )

  result_2 <-
    phenograph_data %>%
    tof_downsample_density(
      group_cols = c(sample_name, phenograph_cluster),
      target_percentile = 0.05,
      density_cols = c(cd19, cd11b)
    )

  # number of columns is unchanged
  expect_identical(ncol(result), ncol(phenograph_data))
  expect_identical(ncol(result_2), ncol(phenograph_data))

  # number of rows is correct
  expect_true(nrow(result) < nrow(phenograph_data))
  expect_true(nrow(result_2) < nrow(phenograph_data))
})
