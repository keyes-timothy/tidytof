library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidytof)

data(ddpr_data)

# tof_preprocess ---------------------------------------------------------------

test_that("Shape of transformed data is the same as the input data", {
  result <- tof_preprocess(tof_tibble = ddpr_data)
  expect_true(nrow(result) == nrow(ddpr_data))
  expect_true(ncol(result) == ncol(ddpr_data))
})

test_that("tof_preprocess should give different values when undo_noise is TRUE and FALSE", {
  result <- tof_preprocess(tof_tibble = ddpr_data, undo_noise = TRUE)
  result_2 <- tof_preprocess(tof_tibble = ddpr_data, undo_noise = FALSE)

  expect_false(isTRUE(all.equal(result, result_2)))
})

test_that("tof_preprocess should give different values when transform_fun is different", {
  result <- tof_preprocess(tof_tibble = ddpr_data)
  result_2 <- tof_preprocess(tof_tibble = ddpr_data, transform_fun = scale)

  expect_false(isTRUE(all.equal(result, result_2)))
})

test_that("tof_preprocess should give different values when channel_cols is different", {
  result <- tof_preprocess(tof_tibble = ddpr_data)
  result_2 <- tof_preprocess(tof_tibble = ddpr_data, channel_cols = c(cd45, cd20, cd34))

  expect_false(isTRUE(all.equal(result, result_2)))
})


# tof_postprocess --------------------------------------------------------------

test_that("Shape of transformed data is the same as the input data", {
  result <- tof_postprocess(tof_tibble = ddpr_data)
  expect_true(nrow(result) == nrow(ddpr_data))
  expect_true(ncol(result) == ncol(ddpr_data))
})

test_that("tof_postprocess should give different values when redo_noise is TRUE and FALSE", {
  result <- tof_postprocess(tof_tibble = ddpr_data, redo_noise = TRUE)
  result_2 <- tof_postprocess(tof_tibble = ddpr_data, redo_noise = FALSE)

  expect_false(isTRUE(all.equal(result, result_2)))
})

test_that("tof_postprocess should give different values when transform_fun is different", {
  result <- tof_postprocess(tof_tibble = ddpr_data)
  result_2 <- tof_postprocess(tof_tibble = ddpr_data, transform_fun = scale)

  expect_false(isTRUE(all.equal(result, result_2)))
})

test_that("tof_postprocess should give different values when transform_fun is different", {
  result <- tof_postprocess(tof_tibble = ddpr_data)
  result_2 <- tof_postprocess(tof_tibble = ddpr_data, transform_fun = scale)

  expect_false(isTRUE(all.equal(result, result_2)))
})

test_that("tof_postprocess should give different values when channel_cols is different", {
  result <- tof_postprocess(tof_tibble = ddpr_data)
  result_2 <- tof_postprocess(tof_tibble = ddpr_data, channel_cols = c(cd45, cd20, cd34))

  expect_false(isTRUE(all.equal(result, result_2)))
})

