library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidytof)
library(testthat)

# setup

dat <-
  ddpr_data %>%
  group_by(sample_name) %>%
  slice_sample(n = 20) %>%
  ungroup() %>%
  mutate(
    condition = str_remove(sample_name, "_Basal"),
    replicate = rep(x = c("1", "2"), times = (nrow(.) / 2)),
    cluster = rep(x = c("a", "a", "b", "b"), times = (nrow(.) / 4))
  )

# testing

# tof_extract_proportion -------------------------------------------------------

prop_1 <-
  dat %>%
  tof_extract_proportion(cluster_col = cluster, group_cols = condition)

prop_2 <-
  dat %>%
  tof_extract_proportion(cluster_col = cluster, group_cols = c(condition, replicate))

prop_no_groups <-
  dat %>%
  tof_extract_proportion(cluster_col = cluster)

prop_long <-
  dat %>%
  tof_extract_proportion(cluster_col = cluster, group_cols = c(condition, replicate), format = "long")

prop_long_no_groups <-
  dat %>%
  tof_extract_proportion(cluster_col = cluster, format = "long")

my_list <-
  ls()[str_detect(ls(), "^prop_")]

global_env <- environment()

test_that("result is a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})

test_that("results have the right shape", {
  expect_equal(dim(prop_1), c(2, 3))
  expect_equal(dim(prop_2), c(4, 4))
  expect_equal(dim(prop_no_groups), c(1, 2))
  expect_equal(dim(prop_long), c(8, 4))
  expect_equal(dim(prop_long_no_groups), c(2, 2))
})

test_that("All proportions are 0.5", {
  expect_true(
    all(
      map_lgl(
        .x = mget(my_list, envir = global_env),
        .f = ~
          .x %>%
          select(where(is.numeric)) %>%
          as.matrix() %>%
          `==`(0.5) %>%
          all()
      )
    )
  )
})


# tof_extract_central_tendency -------------------------------------------------
ct_1 <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34)
    )

ct_2 <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    group_cols = c(condition, replicate),
    marker_cols = c(cd45, cd34)
    )

ct_no_groups <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34)
    )

ct_no_groups_stim <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate
  )

ct_stim <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34),
    group_cols = condition,
    stimulation_col = replicate
  )

ct_long <-
  dat %>%
  tof_extract_central_tendency(
    cluster_col = cluster,
    group_cols = condition,
    format = "long",
    marker_cols = c(cd45, cd34)
    )

my_list <-
  ls()[str_detect(ls(), "^ct")]

global_env <- environment()

test_that("result is a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})

test_that("results have the right shape", {
  expect_equal(dim(ct_1), c(2, 5))
  expect_equal(dim(ct_2), c(4, 6))
  expect_equal(dim(ct_no_groups), c(1, 4))
  expect_equal(dim(ct_long), c(8, 4))
})

# tof_extract_threshold --------------------------------------------------------

thresh_1 <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    thresh = 10
  )

thresh_1b <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    threshold = 2
  )

thresh_2 <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    group_cols = c(condition, replicate),
    marker_cols = c(cd45, cd34)
  )

thresh_3 <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    group_cols = condition,
    stimulation_col = replicate,
    marker_cols = c(cd45, cd34)
  )

thresh_no_groups <-
  dat %>%
  tof_extract_threshold(cluster_col = cluster, marker_cols = c(cd45, cd34))

thresh_no_groups_stim <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate
  )

thresh_long <-
  dat %>%
  tof_extract_threshold(
    cluster_col = cluster,
    group_cols = condition,
    format = "long",
    marker_cols = c(cd45, cd34)
  )

my_list <-
  ls()[str_detect(ls(), "^thresh")]

test_that("result is a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})

test_that("results have the right shape", {
  expect_equal(dim(thresh_1), c(2, 5))
  expect_equal(dim(thresh_1b), c(2, 5))
  expect_equal(dim(thresh_2), c(4, 6))
  expect_equal(dim(thresh_3), c(2, 9))
  expect_equal(dim(thresh_long), c(8, 4))
  expect_equal(dim(thresh_no_groups), c(1, 4))
  expect_equal(dim(thresh_no_groups_stim), c(1, 8))
})

test_that("changing the threshold changes the numeric values", {
  values_1 <-
    thresh_1 %>%
    select(where(is.numeric)) %>%
    as.matrix()
  values_2 <-
    thresh_1b %>%
    select(where(is.numeric)) %>%
    as.matrix()
  expect_false(all(values_1 == values_2))
})

test_that("All proportions are between 0 and 1", {
  expect_true(
    all(
      map_lgl(
        .x = mget(my_list, envir = global_env),
        .f = ~
          .x %>%
          select(where(is.numeric)) %>%
          as.matrix() %>%
          (function(x) x >= 0 & x <= 1)() %>%
          all()
      )
    )
  )
})




# tof_extract_emd --------------------------------------------------------------

emd_1 <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_emd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

emd_1b <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_emd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 30
  )

emd_2 <-
  dat %>%
  mutate(condition_1 = "a", condition_2 = "b") %>%
  tof_extract_emd(
    cluster_col = cluster,
    group_cols = c(condition_1, condition_2),
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

emd_3 <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_emd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "2",
    num_bins = 3
  )

emd_no_groups <-
  dat %>%
  tof_extract_emd(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

emd_long <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_emd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    format = "long",
    num_bins = 3
  )

my_list <-
  ls()[str_detect(ls(), "^emd")]

test_that("result is a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})

test_that("results have the right shape", {
  expect_equal(dim(emd_1), c(1, 5))
  expect_equal(dim(emd_2), c(1, 6))
  expect_equal(dim(emd_3), c(1, 5))
  expect_equal(dim(emd_long), c(4, 5))
  expect_equal(dim(emd_no_groups), c(1, 4))
})

test_that("changing the num_dims changes the numeric values", {
  values_1 <-
    emd_1 %>%
    select(where(is.numeric)) %>%
    as.matrix()
  values_2 <-
    emd_1b %>%
    select(where(is.numeric)) %>%
    as.matrix()
  expect_false(all(values_1 == values_2))
})

test_that("column names change when you change basal_level", {
  expect_false(identical(colnames(emd_1), colnames))
})

test_that("EMDs are symmetric", {
  values_1 <-
    emd_1 %>%
    select(where(is.numeric)) %>%
    as.matrix() %>%
    round(digits = 3) %>%
    as.numeric()
  values_2 <-
    emd_3 %>%
    select(where(is.numeric)) %>%
    as.matrix() %>%
    round(digits = 3) %>%
    as.numeric()
  expect_equal(values_1, values_2)
})

test_that("errors are thrown when required arguments are omitted", {
  # stimulation_col argument is missing
  expect_error(
    tof_extract_emd(dat, cluster_col = cluster)
  )

  # basal_level argument is missing
  expect_error(
    tof_extract_emd(dat, cluster_col = cluster, stimulation_col = replicate)
  )
})


# tof_extract_jsd --------------------------------------------------------------


jsd_1 <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_jsd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

jsd_1b <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_jsd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 30
  )

jsd_2 <-
  dat %>%
  mutate(condition_1 = "a", condition_2 = "b") %>%
  tof_extract_jsd(
    cluster_col = cluster,
    group_cols = c(condition_1, condition_2),
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

jsd_3 <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_jsd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "2",
    num_bins = 3
  )

jsd_no_groups <-
  dat %>%
  tof_extract_jsd(
    cluster_col = cluster,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    num_bins = 3
  )

jsd_long <-
  dat %>%
  mutate(condition = "a") %>%
  tof_extract_jsd(
    cluster_col = cluster,
    group_cols = condition,
    marker_cols = c(cd45, cd34),
    stimulation_col = replicate,
    basal_level = "1",
    format = "long",
    num_bins = 3
  )

my_list <-
  ls()[str_detect(ls(), "^jsd")]

test_that("result is a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})

test_that("results have the right shape", {
  expect_equal(dim(jsd_1), c(1, 5))
  expect_equal(dim(jsd_2), c(1, 6))
  expect_equal(dim(jsd_3), c(1, 5))
  expect_equal(dim(jsd_long), c(4, 5))
  expect_equal(dim(jsd_no_groups), c(1, 4))
})

test_that("changing the num_dims changes the numeric values", {
  values_1 <-
    jsd_1 %>%
    select(where(is.numeric)) %>%
    as.matrix()
  values_2 <-
    jsd_1b %>%
    select(where(is.numeric)) %>%
    as.matrix()
  expect_false(all(values_1 == values_2))
})

test_that("column names change when you change basal_level", {
  expect_false(identical(colnames(jsd_1), colnames))
})

test_that("JSDs are symmetric", {
  values_1 <-
    jsd_1 %>%
    select(where(is.numeric)) %>%
    as.matrix() %>%
    round(digits = 3) %>%
    as.numeric()
  values_2 <-
    jsd_3 %>%
    select(where(is.numeric)) %>%
    as.matrix() %>%
    round(digits = 3) %>%
    as.numeric()
  expect_equal(values_1, values_2)
})

test_that("errors are thrown when required arguments are omitted", {
  # stimulation_col argument is missing
  expect_error(
    tof_extract_jsd(dat, cluster_col = cluster)
  )

  # basal_level argument is missing
  expect_error(
    tof_extract_jsd(dat, cluster_col = cluster, stimulation_col = replicate)
  )
})



# tof_extract_features ---------------------------------------------------------

# feature extraction with threshold signaling method
feat_1 <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    group_cols = group,
    stimulation_col = replicate,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt),
    basal_level = "1"
  )

# same with emd method
feat_2 <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    group_cols = group,
    stimulation_col = replicate,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt),
    signaling_method = "emd",
    basal_level = "1"
  )

# same with ct method
feat_3 <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    group_cols = group,
    stimulation_col = replicate,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt),
    signaling_method = "central tendency",
    basal_level = "1"
  )

# same as 1 without any groups
feat_no_groups <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    stimulation_col = replicate,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt),
    basal_level = "1"
  )

# same as 1 without stimulation
feat_no_stim <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    group_cols = group,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt)
  )

# same as 1 without stimulation or groups
feat_no_stim_no_groups <-
  dat %>%
  mutate(group = "x") %>%
  tof_extract_features(
    cluster_col = cluster,
    lineage_cols = c(cd45, cd19),
    signaling_cols = c(pstat5, pakt)
  )

# tests

my_list <-
  ls()[str_detect(ls(), "^feat")]

test_that("all results are a tibble", {
  expect_true(
    all(
      map_lgl(.x = mget(my_list, envir = global_env), .f = ~ "tbl_df" %in% class(.x))
    )
  )
})


test_that("results have the right shape", {
  expect_equal(dim(feat_1), c(1, 15))
  expect_equal(dim(feat_2), c(1, 15))
  expect_equal(dim(feat_3), c(1, 11))
  expect_equal(dim(feat_no_groups), c(1, 14))
  expect_equal(dim(feat_no_stim), c(1, 11))
  expect_equal(dim(feat_no_stim_no_groups), c(1, 10))
})

test_that("Trying emd or jsd without stimulation_col or basal_level fails", {
  # no stimulation_col
  expect_error(
    tof_extract_features(dat, cluster, signaling_method = "emd")
  )
  expect_error(
    tof_extract_features(dat, cluster, signaling_method = "jsd")
  )

  # no basal_level
  expect_error(
    tof_extract_features(dat, cluster, signaling_col = replicate, signaling_method = "emd")
  )
  expect_error(
    tof_extract_features(dat, cluster, signaling_col = replicate, signaling_method = "jsd")
  )
})


