library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidytof)
library(testthat)

# tof_read_fcs -----------------------------------------------------------------

# setup
fcs_path_1 <-
  tidytof_example_data("phenograph") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_2 <-
  tidytof_example_data("ddpr") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_3 <-
  tidytof_example_data("aml") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_4 <-
  tidytof_example_data("surgery") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_5 <-
  tidytof_example_data("surgery") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_6 <-
  tidytof_example_data("surgery") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

fcs_path_list <- unlist(mget(x = paste0("fcs_path_", 1:6)))

csv_path <-
  tidytof_example_data("phenograph_csv") %>%
  dir(full.names = TRUE) %>%
  pluck(1)

# tof_read_fcs -----------------------------------------------------------------

invisible(capture.output(
  ff_tibble <-
    tibble(
      ffs = map(fcs_path_list, flowCore::read.FCS, truncate_max_range = FALSE),
      num_cells = map_int(ffs, nrow),
      num_cols = map_int(ffs, ncol)
    )
))

tof_tibble_list_1 <-
  map(
    .x = fcs_path_list,
    .f = tidytof:::tof_read_fcs
  )

tof_tibble_list_2 <-
  map(
    .x = fcs_path_list,
    .f = tidytof:::tof_read_fcs,
    sep = "___"
  )

test_that("tof_read_fcs() reads an .fcs file into a tof_tbl", {
  class_list <-
    map(.x = tof_tibble_list_1, class)
  lgl_list <- map_lgl(.x = class_list, ~ ("tof_tbl" %in% .x))

  expect_true(all(lgl_list), TRUE)
})

test_that("tof_read_fcs() tof_tbl's have a (correct) panel attribute", {
  panel_tibble <-
    tibble(
      tof_tbl = tof_tibble_list_1,
      panel = map(.x = tof_tbl, .f = tof_get_panel),
      is_panel_tibble = map_lgl(.x = panel, .f = is.tbl),
      panel_has_all_metals =
        map_lgl(.x = panel, .f = ~ all(!is.na(.x$metals)) & all(!(.x$metals == ""))),
      panel_has_all_antigens =
        map_lgl(.x = panel, .f = ~ all(!is.na(.x$antigens)) & all(!(.x$metals == ""))),
      num_channels = map_int(.x = panel, .f = nrow)
    )

  # all panels are a tibble
  expect_true(all(panel_tibble$is_panel_tibble), TRUE)

  # no antigens or metals are NA or "" in any panel
  expect_true(all(panel_tibble$panel_has_all_metals))
  expect_true(all(panel_tibble$panel_has_all_antigens))

  # panel has the expected number of channels
  expect_true(all(panel_tibble$num_channels == ff_tibble$num_cols))

})

test_that("tof_read_fcs() tof_tbl's have the correct number of rows and columns", {
  num_rows <- map_int(tof_tibble_list_1, nrow)
  num_cols <- map_int(tof_tibble_list_1, ncol)

  expect_true(all(num_rows == ff_tibble$num_cells))
  expect_true(all(num_cols == ff_tibble$num_cols))
})

test_that("tof_read_fcs() tof_tbl's have correctly-named columns", {
  colnames_list_1 <- map(tof_tibble_list_1, colnames)
  colnames_list_2 <- map(tof_tibble_list_2, colnames)

  # no colnames are blank
  expect_true(all(map_lgl(colnames_list_1, ~ all(str_length(.x) > 0))))
  expect_true(all(map_lgl(colnames_list_2, ~ all(str_length(.x) > 0))))

  # no colnames are NA
  expect_true(all(map_lgl(colnames_list_1, ~ all(!is.na(.x)))))
  expect_true(all(map_lgl(colnames_list_2, ~ all(!is.na(.x)))))

  # the sep argument works
  expect_true(all(map_lgl(colnames_list_2, ~ all(str_detect(.x, pattern = "___")))))
})

# tof_read_csv -----------------------------------------------------------------

# setup
csv_tibble <-
  csv_path %>%
  tidytof:::tof_read_csv()

csv_tibble_2 <-
  csv_path %>%
  tidytof:::tof_read_csv(
    panel_info =
      tibble(
        antigens = c("CD33", "CD3", "CD333"),
        metals = c("Pd104", "Pd106", "Pd108")
      )
  )

# tests

test_that("tof_read_csv() reads a .csv file into a tof_tbl", {
  expect_true("tof_tbl" %in% class(csv_tibble))
})

test_that("tof_read_csv() tof_tbl's have a (correct) panel attribute", {
  csv_panel <-
    csv_tibble %>%
    tof_get_panel()

  csv_panel_2 <-
    csv_tibble_2 %>%
    tof_get_panel()

  # all panels are a tibble
  expect_true("tbl" %in% class(csv_panel))
  expect_true("tbl" %in% class(csv_panel_2))

  # if the tibble is not empty, it requires the `antigens` and `metals` columns
  expect_true(
    identical(csv_panel, tibble()) | all(c("antigens", "metals") %in% colnames(csv_panel))
  )
  expect_true(
    identical(csv_panel_2, tibble()) | all(c("antigens", "metals") %in% colnames(csv_panel_2))
  )
})

# tof_read_file ----------------------------------------------------------------

# setup

fcs_fcs <-
  fcs_path_1 %>%
  tidytof:::tof_read_fcs()

file_fcs <-
  fcs_path_1 %>%
  tidytof:::tof_read_file()

csv_csv <-
  csv_path %>%
  tidytof:::tof_read_csv()

file_csv <-
  csv_path %>%
  tidytof:::tof_read_file()

test_that("tof_read_fcs and tof_read_file produce identical results for .fcs files", {
  expect_equal(fcs_fcs, file_fcs, ignore_attr = TRUE)
})

test_that("tof_read_csv and tof_read_file produce identical results for .csv files", {
  expect_equal(csv_csv, file_csv, ignore_attr = TRUE)
})



# tof_read_data ----------------------------------------------------------------

test_that("tof_read_data can read in a single .fcs file", {
  expect_equal(tof_read_data(fcs_path_1), tidytof:::tof_read_fcs(fcs_path_1), ignore_attr = TRUE)
})

test_that("tof_read_data can read in a single .csv file", {
  expect_equal(tof_read_data(csv_path), tidytof:::tof_read_csv(csv_path), ignore_attr = TRUE)
})

test_that("tof_read_data can read in multiple .fcs files", {
  tof_tibble <- tof_read_data(tidytof_example_data("phenograph"))
  tof_tibble_mixed <- tof_read_data(tidytof_example_data("mix"))

  # unmixed result is a tof_tibble
  expect_s3_class(tof_tibble, "tof_tbl")

  # unmixed result has a file_name column
  expect_true("file_name" %in% colnames(tof_tibble))

  # mixed result is a normal tibble
  expect_s3_class(tof_tibble_mixed, "tbl_df")

  # mixed result has a panel and data column
  expect_true(all(c("panel", "data") %in% colnames(tof_tibble_mixed)))

  # each entry of the data column of the mixed result is a tof_tibble
  expect_true(all(map_lgl(tof_tibble_mixed$data, ~ "tof_tbl" %in% class(.x))))

  # each entry of the panel column of the mixed result is a normal tibble
  expect_true(all(map_lgl(tof_tibble_mixed$panel, ~ "tbl" %in% class(.x))))

  # mixed result has the correct number of rows
  expect_equal(nrow(tof_tibble_mixed), 2L)
})

test_that("tof_read_data can read in multiple .csv files", {
  tof_tibble <- tof_read_data(tidytof_example_data("phenograph"))

  # result is a tof_tibble
  expect_s3_class(tof_tibble, "tof_tbl")

  # result has a file_name column
  expect_true("file_name" %in% colnames(tof_tibble))

})

test_that("tof_read_data can read in multiple .csv and .fcs files simultaneously", {
  my_panel <-
    tidytof_example_data("phenograph") %>%
    dir(full.names = TRUE) %>%
    pluck(1) %>%
    tidytof:::tof_read_fcs() %>%
    tof_get_panel()

  tof_tibble <- tof_read_data(tidytof_example_data("mix2"), panel = my_panel)

  # result has the right number of rows
  expect_equal(nrow(tof_tibble), 2L)

  # all panels are tibbles
  expect_true(all(map_lgl(tof_tibble$panel, is.tbl)))

  # all data entries are tof_tbls
  expect_true(all(map_lgl(tof_tibble$data, ~"tof_tbl" %in% class(.x))))
})


