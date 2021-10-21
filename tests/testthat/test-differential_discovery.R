library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidytof)
library(testthat)

# setup
set.seed(2020)
dd_data <-
  ddpr_data %>%
  mutate(
    stim = sample(c("stim", "basal"), size = nrow(ddpr_data), replace = TRUE),
    stim_2 = sample(c("a", "b", "c"), size = nrow(ddpr_data), replace = TRUE),
    condition = if_else(str_detect(sample_name, "Healthy"), "healthy", "cancer"),
    cell_id = 1:nrow(ddpr_data),
    replicate = as.character(round(cell_id, digits = -3)),
    sample_name = str_c(stim, condition, replicate, sep = "_"),
    sample_name_2 = str_c(stim_2, condition, replicate, sep = "_"),
    pairs = str_c(condition, replicate, sep = "_")
  ) %>%
  select(-cell_id) %>%
  tof_preprocess() %>%
  tof_cluster(method = "flowsom")

# tof_daa_diffcyt --------------------------------------------------------------

### test basic results
diffcyt_glmm <-
  dd_data %>%
  tof_daa_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = condition,
    random_effect_cols = stim,
    method = "glmm"
  )

diffcyt_glmm_fixed <-
  dd_data %>%
  tof_daa_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = condition,
    method = "glmm",
    min_cells = 3,
    min_samples = 2
  )

diffcyt_voom <-
  dd_data %>%
  tof_daa_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = condition,
    random_effect_cols = stim,
    method = "voom"
  )

diffcyt_voom_fixed <-
  dd_data %>%
  tof_daa_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = condition,
    method = "voom"
  )

diffcyt_edgeR <-
  dd_data %>%
  tof_daa_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = condition,
    method = "edgeR"
  )

test_that("diffcyt results tibbles have the right dimensions", {
  # are tibbles
  expect_s3_class(diffcyt_glmm, "tbl_df")
  expect_s3_class(diffcyt_voom, "tbl_df")
  expect_s3_class(diffcyt_edgeR, "tbl_df")

  # have correct number of rows
  expect_equal(nrow(diffcyt_glmm), 5L)
  expect_equal(nrow(diffcyt_voom), 5L)
  expect_equal(nrow(diffcyt_edgeR), 5L)

  # have correct number of columns
  expect_equal(ncol(diffcyt_glmm), 5L)
  expect_equal(ncol(diffcyt_voom), 9L)
  expect_equal(ncol(diffcyt_edgeR), 8L)
})

test_that("diffcyt daa_results have the right content", {
  # glmm
  expect_true(all(
      c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(diffcyt_glmm)
    ))

  # voom
  expect_true(all(
      c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(diffcyt_voom)
    ))
  expect_true(all(
    c("logFC", "AveExpr", "t", "B") %in% colnames(diffcyt_voom)
  ))

  # edgeR
  expect_true(all(
    c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(diffcyt_edgeR)
  ))
  expect_true(all(
    c("logFC", "logCPM", "LR") %in% colnames(diffcyt_edgeR)
  ))
})

### test arguments

test_that("edgeR throws an error when you try to model random effects", {
  expect_error(
    dd_data %>%
      tof_daa_diffcyt(
        sample_col = sample_name,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = condition,
        random_effect_cols = stim,
        method = "edgeR"
      )
  )
})

test_that("OLRE give a warning for any method that isn't GLMMs", {
  expect_warning(
    dd_data %>%
      tof_daa_diffcyt(
        sample_col = sample_name,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = condition,
        method = "edgeR",
        include_observation_level_random_effects = TRUE
      )
  )

  expect_warning(
    dd_data %>%
      tof_daa_diffcyt(
        sample_col = sample_name,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = condition,
        method = "voom",
        include_observation_level_random_effects = TRUE
      )
  )
})

test_that("Using only a fixed-effect doesn't break anything", {
  # are tibbles
  expect_s3_class(diffcyt_glmm_fixed, "tbl_df")
  expect_s3_class(diffcyt_voom_fixed, "tbl_df")

  # have correct number of rows
  expect_equal(nrow(diffcyt_glmm_fixed), nrow(diffcyt_glmm))
  expect_equal(nrow(diffcyt_voom_fixed), nrow(diffcyt_voom))

  # have correct number of columns
  expect_equal(ncol(diffcyt_glmm_fixed), ncol(diffcyt_glmm))
  expect_equal(ncol(diffcyt_voom_fixed), ncol(diffcyt_voom))

  # content
  expect_identical(colnames(diffcyt_glmm), colnames(diffcyt_glmm_fixed))
  expect_identical(colnames(diffcyt_voom), colnames(diffcyt_voom_fixed))
})

test_that("Using a fixed-effect with more than 2 levels results in a nested output tibble", {
  expect_equal(
    dd_data %>%
      tof_daa_diffcyt(
        sample_col = sample_name_2,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = stim_2,
        method = "glmm"
      ) %>%
      nrow(),
    3L
  )
})



# tof_daa_glmm -----------------------------------------------------------------

daa_glmm <-
  dd_data %>%
  tof_daa_glmm(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = stim,
    random_effect_cols = replicate
  )

daa_glmm_fixed <-
  dd_data %>%
  tof_daa_glmm(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    fixed_effect_cols = stim
  )

test_that("glmm results tibbles have the right dimensions", {
  # are tibbles
  expect_s3_class(daa_glmm, "tbl_df")
  expect_s3_class(daa_glmm_fixed, "tbl_df")

  # have correct number of rows
  expect_equal(nrow(daa_glmm), 5L)
  expect_equal(nrow(daa_glmm_fixed), 5L)

  # have correct number of columns
  expect_equal(ncol(daa_glmm), 9L)
  expect_equal(ncol(daa_glmm_fixed), 9L)

})


test_that("glmm daa_results have the right content", {
  expect_true(all(
    c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(daa_glmm)
  ))

  expect_identical(colnames(daa_glmm), colnames(daa_glmm_fixed))
})


test_that("Using a fixed-effect with more than 2 levels results in a larger output tibble", {
  expect_equal(
    dd_data %>%
      tof_daa_glmm(
        sample_col = sample_name_2,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = stim_2
      ) %>%
      nrow(),
    2L
  )
})

# tof_daa_ttest ----------------------------------------------------------------

# setup
daa_t <-
  dd_data %>%
  tof_daa_ttest(
    cluster_col = .flowsom_metacluster,
    effect_col = condition,
    group_cols = replicate,
    min_cells = 2,
    min_samples = 1
  )

daa_t_paired <-
  dd_data %>%
  tof_daa_ttest(
    cluster_col = .flowsom_metacluster,
    effect_col = stim,
    group_cols = c(condition, replicate),
    test_type = "paired"
  )

test_that("results have the right shape", {
  expect_s3_class(daa_t, "tbl_df")
  expect_s3_class(daa_t_paired, "tbl_df")

  expect_equal(ncol(daa_t), ncol(daa_t_paired))
})

test_that("results have the right content", {
  expect_true(
    all(
      c("t", "df", "p_val", "p_adj", "significant", "mean_diff", "mean_fc") %in%
        colnames(daa_t)
    )
  )
  expect_identical(colnames(daa_t), colnames(daa_t_paired))

  # check that the column representing the cluster id of the cluster being
  # tested is a character vector (per tidytof style)
  expect_type(daa_t[,1][[1]], "character")
  expect_type(daa_t_paired[,1][[1]], "character")

  # look for NaN values
  expect_false(any(is.nan(daa_t$mean_diff)))
  expect_false(any(is.nan(daa_t$mean_fc)))

  expect_false(any(is.nan(daa_t_paired$mean_diff)))
  expect_false(any(is.nan(daa_t_paired$mean_fc)))
})


test_that("choosing effect_col with more than 2 levels throws an error", {
  expect_error(
    dd_data %>%
      tof_daa_ttest(
        cluster_col = .flowsom_metacluster,
        effect_col = stim_2,
        group_names = replicate
      )
  )
})

test_that("omitting group_cols should give an error", {
  # unpaired
  expect_error(
      dd_data %>%
      tof_daa_ttest(
        cluster_col = .flowsom_metacluster,
        effect_col = condition,
        test_type = "unpaired",
        min_cells = 1,
        min_samples = 0,
      )
  )

  # paired
  expect_error(
      dd_data %>%
      tof_daa_ttest(
        cluster_col = .flowsom_metacluster,
        effect_col = stim,
        test_type = "paired",
        min_cells = 1,
        min_samples = 0,
      )
  )

})


# tof_dea_diffcyt --------------------------------------------------------------

# setup
diffcyt_lmm <-
  dd_data %>%
  tof_dea_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = condition,
    random_effect_cols = stim,
    method = "lmm"
  )

diffcyt_lmm_fixed <-
  dd_data %>%
  tof_dea_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = condition,
    method = "lmm"
  )

diffcyt_lmm_tidyselection <-
  dd_data %>%
  tof_dea_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = starts_with("cd", ignore.case = FALSE),
    fixed_effect_cols = condition,
    random_effect_cols = stim,
    method = "lmm"
  )

diffcyt_limma <-
  dd_data %>%
  tof_dea_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = condition,
    random_effect_cols = stim,
    method = "limma"
  )

diffcyt_limma_fixed <-
  dd_data %>%
  tof_dea_diffcyt(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = condition,
    method = "limma"
  )


test_that("diffcyt result tibbles have the right dimensions", {
  # are tibbles
  expect_s3_class(diffcyt_lmm, "tbl_df")
  expect_s3_class(diffcyt_lmm_fixed, "tbl_df")
  expect_s3_class(diffcyt_limma, "tbl_df")
  expect_s3_class(diffcyt_limma_fixed, "tbl_df")

  # have correct number of rows
  expect_equal(nrow(diffcyt_lmm), 10L)
  expect_equal(nrow(diffcyt_lmm_fixed), 10L)
  expect_equal(nrow(diffcyt_limma), 10L)
  expect_equal(nrow(diffcyt_limma_fixed), 10L)

  # have correct number of columns
  expect_equal(ncol(diffcyt_lmm), 6L)
  expect_equal(ncol(diffcyt_lmm_fixed), 6L)
  expect_equal(ncol(diffcyt_limma), 10L)
  expect_equal(ncol(diffcyt_limma_fixed), 10L)

})

test_that("diffcyt dea_results have the right content", {
  # lmm
  expect_true(all(
    c(".flowsom_metacluster", "marker", "p_val", "p_adj") %in%
      colnames(diffcyt_lmm)
  ))

  # limma
  expect_true(all(
    c(".flowsom_metacluster", "marker", "p_val", "p_adj") %in%
      colnames(diffcyt_limma)
  ))
  expect_true(all(
    c("logFC", "AveExpr", "t", "B") %in% colnames(diffcyt_limma)
  ))
})

### test arguments

test_that("OLRE give a warning for any method that isn't GLMMs", {
  expect_warning(
    dd_data %>%
      tof_dea_diffcyt(
        sample_col = sample_name,
        cluster_col = .flowsom_metacluster,
        fixed_effect_cols = condition,
        method = "limma",
        include_observation_level_random_effects = TRUE
      )
  )
})

test_that("Using only a fixed-effect doesn't break anything", {
  # are tibbles
  expect_s3_class(diffcyt_lmm_fixed, "tbl_df")
  expect_s3_class(diffcyt_limma_fixed, "tbl_df")

  # content
  expect_true(all(
    c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(diffcyt_lmm_fixed)
  ))
  expect_true(all(
    c(".flowsom_metacluster", "p_val", "p_adj") %in% colnames(diffcyt_limma_fixed)
  ))
  expect_true(all(
    c("logFC", "AveExpr", "t", "B") %in% colnames(diffcyt_limma_fixed)
  ))
})

test_that("Using a fixed-effect with more than 2 levels results in a larger output tibble", {
  expect_equal(
    dd_data %>%
      tof_dea_diffcyt(
        sample_col = sample_name_2,
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        fixed_effect_cols = stim_2,
        method = "lmm"
      ) %>%
      nrow(),
    3L
  )
})


# tof_dea_lmm ------------------------------------------------------------------

# setup

dea_lmm <-
  dd_data %>%
  tof_dea_lmm(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = stim,
    random_effect_cols = replicate
  )

dea_lmm_fixed <-
  dd_data %>%
  tof_dea_lmm(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = stim
  )

dea_lmm_tidyselection <-
  dd_data %>%
  tof_dea_lmm(
    sample_col = sample_name,
    cluster_col = .flowsom_metacluster,
    marker_cols = all_of(c("cd45", "cd19", "cd127")),
    fixed_effect_cols = stim,
    random_effect_cols = replicate
  )


test_that("lmm result tibbles have the right dimensions", {
  # are tibbles
  expect_s3_class(dea_lmm, "tbl_df")
  expect_s3_class(dea_lmm_fixed, "tbl_df")
  expect_s3_class(dea_lmm_tidyselection, "tbl_df")

  # have correct number of rows
  expect_equal(nrow(dea_lmm), 10L)
  expect_equal(nrow(dea_lmm_fixed), 10L)
  expect_equal(nrow(dea_lmm_tidyselection), 15L)

  # have correct number of columns
  expect_equal(ncol(dea_lmm), 11L)
  expect_equal(ncol(dea_lmm_fixed), 10L)
  expect_equal(ncol(dea_lmm_tidyselection), 11L)

})

test_that("lmm dea_results have the right content", {
  # lmm
  expect_true(all(
    c(".flowsom_metacluster", "marker", "p_val", "p_adj", "mean_diff") %in%
      colnames(dea_lmm)
  ))
})


test_that("Using a fixed-effect with more than 2 levels results in a larger output tibble", {
  expect_equal(
    dd_data %>%
      tof_dea_lmm(
        sample_col = sample_name_2,
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        fixed_effect_cols = stim_2
      ) %>%
      nrow(),
    2L
  )
})


# tof_dea_ttest ----------------------------------------------------------------

#setup
dea_t_unpaired <-
  dd_data %>%
  tof_dea_ttest(
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    effect_col = condition,
    group_cols = replicate,
    test_type = "unpaired",
    summary_function = mean,
    quiet = TRUE
  )

dea_t_paired <-
  dd_data %>%
  tof_dea_ttest(
    cluster_col = .flowsom_metacluster,
    marker_cols = c(cd45, cd19),
    effect_col = stim,
    group_cols = replicate,
    test_type = "paired",
    summary_function = mean,
    quiet = TRUE
  )


test_that("t-test dea results are shaped correctly", {
  # rows
  num_clusters <-
    dd_data %>%
    distinct(.flowsom_metacluster) %>%
    nrow()

  expect_equal(nrow(dea_t_unpaired), num_clusters * 2)
  expect_equal(nrow(dea_t_paired), num_clusters * 2)

  # columns
  expect_equal(ncol(dea_t_unpaired), 9L)
  expect_equal(ncol(dea_t_paired), 9L)
})

test_that("t-test dea results have the correct columns", {
  expect_true(all(
    c(".flowsom_metacluster", "marker", "t",
      "df", "p_val", "p_adj", "significant",
      "mean_diff", "mean_fc") %in%
      colnames(dea_t_unpaired)
  ))

  expect_identical(colnames(dea_t_unpaired), colnames(dea_t_paired))
})

test_that("p-values are arranged in ascending order", {
  expect_equal(sort(dea_t_unpaired$p_adj, na.last = TRUE), dea_t_unpaired$p_adj)
  expect_equal(sort(dea_t_paired$p_adj, na.last = TRUE), dea_t_paired$p_adj)
})

test_that("changing the summary function works", {
  expect_s3_class(
    dd_data %>%
      tof_dea_ttest(
        cluster_col = .flowsom_metacluster,
        marker_cols = cd45,
        effect_col = condition,
        group_cols = replicate,
        test_type = "unpaired",
        summary_function = median,
        quiet = TRUE
      ),
    class = "tbl_df"
  )

  expect_s3_class(
    dd_data %>%
      tof_dea_ttest(
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        effect_col = stim,
        group_cols = replicate,
        test_type = "paired",
        summary_function = median,
        quiet = TRUE
      ),
    class = "tbl_df"
  )
})


test_that("omitting group_cols should throw an error", {
  # unpaired
  expect_error(
      dd_data %>%
      tof_dea_ttest(
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        effect_col = condition,
        test_type = "unpaired",
        summary_function = median,
        min_cells = 1,
        min_samples = 0,
      )
  )

  # paired
  expect_error(
      dd_data %>%
      tof_dea_ttest(
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        effect_col = stim,
        test_type = "paired",
        summary_function = median,
        min_cells = 1,
        min_samples = 0,
      )
  )

})

test_that("lack of effect_col throws an error", {
  expect_error(
    dd_data %>%
      tof_dea_ttest(
        cluster_col = .flowsom_metacluster,
        marker_cols = c(cd45, cd19),
        group_cols = replicate,
        test_type = "paired",
        summary_function = median,
        min_cells = 1,
        min_samples = 0,
      )
  )
})
