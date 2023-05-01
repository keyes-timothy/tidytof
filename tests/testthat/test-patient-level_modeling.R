library(dplyr)
library(purrr)
library(tidytof)
library(testthat)

# setup ------------------------------------------------------------------------
feature_tibble <-
  dplyr::tibble(
    sample = as.character(1:100),
    cd45 = runif(n = 100),
    pstat5 = runif(n = 100),
    cd34 = runif(n = 100),
    outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
    class =
      if_else(outcome > median(outcome), "class1", "class2") %>%
      as.factor(),
    multiclass = as.factor(c(rep("class1", 30), rep("class2", 30), rep("class3", 40))),
    event = c(rep(0, times = 50), rep(1, times = 50)),
    time_to_event = rnorm(n = 100, mean = 10, sd = 2)
  )


# tof_split_data ---------------------------------------------------------------

cv_split <-
  feature_tibble %>%
  tof_split_data(
    split_method = "k-fold",
    strata = multiclass
  )

cv_split_2 <-
  feature_tibble %>%
  tof_split_data(split_method = "k-fold", num_cv_repeats = 2L)

boot_split <-
  feature_tibble %>%
  tof_split_data(
    split_method = "bootstrap"
  )

simple_split <-
  feature_tibble %>%
  tof_split_data(
    split_method = "simple"
  )

split_stratum <-
  sample(
    x = c(TRUE, FALSE),
    size = nrow(feature_tibble),
    replace = TRUE,
    prob = c(0.75, 0.25)
  )

simple_split_2 <-
  feature_tibble %>%
  mutate(split_stratum = split_stratum) %>%
  tof_split_data(split_col = split_stratum)

simple_split_3 <-
  feature_tibble %>%
  tof_split_data(split_method = split_stratum)

test_that("S3 classes for resampled objects are correct", {
  expect_s3_class(cv_split, "rset")
  expect_s3_class(boot_split, "rset")
  expect_s3_class(simple_split, "rsplit")
  expect_s3_class(simple_split_2, "rsplit")
})

test_that("Simple split methods give identical results", {
  expect_equal(nrow(simple_split_2), nrow(simple_split_3))
})

test_that("Test that resampled objects have the correct size", {
  expect_equal(nrow(cv_split), 10L)
  expect_equal(nrow(boot_split), 10L)
  expect_equal(nrow(cv_split), 10L)

  expect_identical(colnames(feature_tibble), colnames(rsample::assessment(simple_split)))
  expect_identical(colnames(feature_tibble), colnames(rsample::training(simple_split)))
})

# tof_create_grid --------------------------------------------------------------

# should have  5 * 5 rows
reg_grid_1 <-
  tof_create_grid()

# should have 5 * 4 rows
reg_grid_2 <-
  tof_create_grid(
    penalty_values = c(0, 0.1, 0.2, 0.5)
  )

# should have 3 * 2 rows
reg_grid_3 <-
  tof_create_grid(
    penalty_values = c(0, 0.5, 1),
    mixture_values = c(0.5, 0.8)
  )

test_that("Regular grids have the right size", {
  expect_equal(nrow(reg_grid_1), 25L)
  expect_equal(nrow(reg_grid_2), 20L)
  expect_equal(nrow(reg_grid_3), 6L)

  expect_equal(ncol(reg_grid_1), 2L)
  expect_equal(ncol(reg_grid_2), 2L)
  expect_equal(ncol(reg_grid_3), 2L)
})

test_that("Regular grid columns are correct", {
  expect_equal(colnames(reg_grid_1), c("penalty", "mixture"))
  expect_equal(colnames(reg_grid_2), c("penalty", "mixture"))
  expect_equal(colnames(reg_grid_3), c("penalty", "mixture"))

  expect_true(all(c(is.numeric(reg_grid_1$penalty), is.numeric(reg_grid_1$mixture))))
  expect_true(all(c(is.numeric(reg_grid_2$penalty), is.numeric(reg_grid_2$mixture))))
  expect_true(all(c(is.numeric(reg_grid_3$penalty), is.numeric(reg_grid_3$mixture))))

})


# tof_train_model --------------------------------------------------------------

# linear regression
cv_linear_regression <-
  cv_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear",
    hyperparameter_grid = reg_grid_3
  )

bootstrap_linear_regression <-
  boot_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear",
    hyperparameter_grid = reg_grid_3
  )

simple_linear_regression <-
  simple_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear",
    hyperparameter_grid = reg_grid_3
  )

unsplit_linear_regression <-
  feature_tibble %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear",
    hyperparameter_grid = reg_grid_3
  )

test_that("linear tof_model's have the right components", {
  components <-
    c("model", "recipe", "mixture", "penalty",
      "model_type", "outcome_colnames", "training_data",
      "tuning_metrics")

  # names are correct
  expect_equal(setdiff(components, names(cv_linear_regression)), character())
  expect_equal(setdiff(components, names(bootstrap_linear_regression)), character())
  expect_equal(setdiff(components, names(simple_linear_regression)), character())
  expect_equal(setdiff(components, names(unsplit_linear_regression)), character())

  ### types/classes are correct
  linear_models <- list(cv_linear_regression, bootstrap_linear_regression,
                        simple_linear_regression, unsplit_linear_regression)

  models <- purrr::map_lgl(linear_models, ~ inherits(.x$model, "glmnet"))
  recipes <- purrr::map_lgl(linear_models, ~ inherits(.x$recipe, "recipe"))
  mixtures <- purrr::map_lgl(linear_models, ~ inherits(.x$mixture, "numeric"))
  penalties <- purrr::map_lgl(linear_models, ~ inherits(.x$penalty, "numeric"))
  model_types <-
    purrr::map_lgl(linear_models, ~ inherits(.x$model_type, "character"))
  outcome_colnames <-
    purrr::map_lgl(linear_models, ~ inherits(.x$outcome_colnames, "character"))
  training_datas <-
    purrr::map_lgl(linear_models, ~ inherits(.x$training_data, "tbl_df"))
  tuning_metrics <-
    purrr::map_lgl(linear_models, ~ inherits(.x$training_data, "tbl_df"))

  expect_true(all(models))
  expect_true(all(recipes))
  expect_true(all(mixtures))
  expect_true(all(penalties))
  expect_true(all(model_types))
  expect_true(all(outcome_colnames))
  expect_true(all(training_datas))
  expect_true(all(tuning_metrics))
})

test_that("linear tof_model's training data has the right shape", {

  # number of rows
  expect_equal(
    nrow(feature_tibble),
    nrow(cv_linear_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(bootstrap_linear_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(simple_linear_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(unsplit_linear_regression$training_data)
  )

  # columns
  expect_equal(
    colnames(feature_tibble),
    colnames(cv_linear_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(bootstrap_linear_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(simple_linear_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(unsplit_linear_regression$training_data)
  )
})

test_that("linear tof_model's tuning_metrics has the right shape", {
  linear_tuning_metrics <-
    c("mixture", "penalty", "mse", "mae")

  # columns
  expect_equal(
    c("id", linear_tuning_metrics, ".predictions"),
    colnames(cv_linear_regression$tuning_metrics)
  )
  expect_equal(
    c("id", linear_tuning_metrics, ".predictions"),
    colnames(bootstrap_linear_regression$tuning_metrics)
  )
  expect_equal(
    linear_tuning_metrics,
    colnames(simple_linear_regression$tuning_metrics)
  )
  expect_equal(
    linear_tuning_metrics,
    colnames(unsplit_linear_regression$tuning_metrics)
  )
})


# logistic regression
cv_logistic_regression <-
  cv_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class",
    hyperparameter_grid = reg_grid_3
  )
bootstrap_logistic_regression <-
  boot_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class",
    hyperparameter_grid = reg_grid_3
  )
simple_logistic_regression <-
  simple_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class",
    hyperparameter_grid = reg_grid_3
  )
unsplit_logistic_regression <-
  feature_tibble %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class",
    hyperparameter_grid = reg_grid_3
  )


test_that("logistic tof_model's have the right components", {
  components <-
    c("model", "recipe", "mixture", "penalty",
      "model_type", "outcome_colnames", "training_data",
      "tuning_metrics")

  # names are correct
  expect_equal(setdiff(components, names(cv_logistic_regression)), character())
  expect_equal(setdiff(components, names(bootstrap_logistic_regression)), character())
  expect_equal(setdiff(components, names(simple_logistic_regression)), character())
  expect_equal(setdiff(components, names(unsplit_logistic_regression)), character())

  ### types/classes are correct
  logistic_models <- list(cv_logistic_regression, bootstrap_logistic_regression,
                        simple_logistic_regression, unsplit_logistic_regression)

  models <- purrr::map_lgl(logistic_models, ~ inherits(.x$model, "glmnet"))
  recipes <- purrr::map_lgl(logistic_models, ~ inherits(.x$recipe, "recipe"))
  mixtures <- purrr::map_lgl(logistic_models, ~ inherits(.x$mixture, "numeric"))
  penalties <- purrr::map_lgl(logistic_models, ~ inherits(.x$penalty, "numeric"))
  model_types <-
    purrr::map_lgl(logistic_models, ~ inherits(.x$model_type, "character"))
  outcome_colnames <-
    purrr::map_lgl(logistic_models, ~ inherits(.x$outcome_colnames, "character"))
  training_datas <-
    purrr::map_lgl(logistic_models, ~ inherits(.x$training_data, "tbl_df"))
  tuning_metrics <-
    purrr::map_lgl(logistic_models, ~ inherits(.x$training_data, "tbl_df"))

  expect_true(all(models))
  expect_true(all(recipes))
  expect_true(all(mixtures))
  expect_true(all(penalties))
  expect_true(all(model_types))
  expect_true(all(outcome_colnames))
  expect_true(all(training_datas))
  expect_true(all(tuning_metrics))
})

test_that("logistic tof_model's training data has the right shape", {

  # number of rows
  expect_equal(
    nrow(feature_tibble),
    nrow(cv_logistic_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(bootstrap_logistic_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(simple_logistic_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(unsplit_logistic_regression$training_data)
  )

  # columns
  expect_equal(
    colnames(feature_tibble),
    colnames(cv_logistic_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(bootstrap_logistic_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(simple_logistic_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(unsplit_logistic_regression$training_data)
  )
})

test_that("logistic tof_model's tuning_metrics has the right shape", {
  logistic_tuning_metrics <-
    c("mixture", "penalty", "binomial_deviance", "misclassification_error",
      "roc_auc", "mse", "mae", "accuracy")

  # columns
  expect_equal(
    c("id", logistic_tuning_metrics, ".predictions"),
    colnames(cv_logistic_regression$tuning_metrics)
  )
  expect_equal(
    c("id", logistic_tuning_metrics, ".predictions"),
    colnames(bootstrap_logistic_regression$tuning_metrics)
  )
  expect_equal(
    logistic_tuning_metrics,
    colnames(simple_logistic_regression$tuning_metrics)
  )
  expect_equal(
    logistic_tuning_metrics,
    colnames(unsplit_logistic_regression$tuning_metrics)
  )
})

# multinomial regression
cv_multinomial_regression <-
  cv_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = multiclass,
    model_type = "multiclass",
    hyperparameter_grid = reg_grid_3
  )

bootstrap_multinomial_regression <-
  boot_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = multiclass,
    model_type = "multiclass",
    hyperparameter_grid = reg_grid_3
  )
simple_multinomial_regression <-
  simple_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = multiclass,
    model_type = "multiclass",
    hyperparameter_grid = reg_grid_3
  )
unsplit_multinomial_regression <-
  feature_tibble %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = multiclass,
    model_type = "multiclass",
    hyperparameter_grid = reg_grid_3
  )


test_that("multinomial tof_model's have the right components", {
  components <-
    c("model", "recipe", "mixture", "penalty",
      "model_type", "outcome_colnames", "training_data",
      "tuning_metrics")

  # names are correct
  expect_equal(setdiff(components, names(cv_multinomial_regression)), character())
  expect_equal(setdiff(components, names(bootstrap_multinomial_regression)), character())
  expect_equal(setdiff(components, names(simple_multinomial_regression)), character())
  expect_equal(setdiff(components, names(unsplit_multinomial_regression)), character())

  ### types/classes are correct
  multinomial_models <- list(cv_multinomial_regression, bootstrap_multinomial_regression,
                             simple_multinomial_regression, unsplit_multinomial_regression)

  models <- purrr::map_lgl(multinomial_models, ~ inherits(.x$model, "glmnet"))
  recipes <- purrr::map_lgl(multinomial_models, ~ inherits(.x$recipe, "recipe"))
  mixtures <- purrr::map_lgl(multinomial_models, ~ inherits(.x$mixture, "numeric"))
  penalties <- purrr::map_lgl(multinomial_models, ~ inherits(.x$penalty, "numeric"))
  model_types <-
    purrr::map_lgl(multinomial_models, ~ inherits(.x$model_type, "character"))
  outcome_colnames <-
    purrr::map_lgl(multinomial_models, ~ inherits(.x$outcome_colnames, "character"))
  training_datas <-
    purrr::map_lgl(multinomial_models, ~ inherits(.x$training_data, "tbl_df"))
  tuning_metrics <-
    purrr::map_lgl(multinomial_models, ~ inherits(.x$training_data, "tbl_df"))

  expect_true(all(models))
  expect_true(all(recipes))
  expect_true(all(mixtures))
  expect_true(all(penalties))
  expect_true(all(model_types))
  expect_true(all(outcome_colnames))
  expect_true(all(training_datas))
  expect_true(all(tuning_metrics))
})

test_that("multinomial tof_model's training data has the right shape", {

  # number of rows
  expect_equal(
    nrow(feature_tibble),
    nrow(cv_multinomial_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(bootstrap_multinomial_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(simple_multinomial_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(unsplit_multinomial_regression$training_data)
  )

  # columns
  expect_equal(
    colnames(feature_tibble),
    colnames(cv_multinomial_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(bootstrap_multinomial_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(simple_multinomial_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(unsplit_multinomial_regression$training_data)
  )
})


test_that("multinomial tof_model's tuning_metrics has the right shape", {
  multinomial_tuning_metrics <-
    c("mixture", "penalty", "multinomial_deviance", "misclassification_error",
      "mse", "mae", "roc_auc", "accuracy")

  # columns
  expect_equal(
    c("id", multinomial_tuning_metrics, ".predictions"),
    colnames(cv_multinomial_regression$tuning_metrics)
  )
  expect_equal(
    c("id", multinomial_tuning_metrics, ".predictions"),
    colnames(bootstrap_multinomial_regression$tuning_metrics)
  )
  expect_equal(
    multinomial_tuning_metrics,
    colnames(simple_multinomial_regression$tuning_metrics)
  )
  expect_equal(
    multinomial_tuning_metrics,
    colnames(unsplit_multinomial_regression$tuning_metrics)
  )
})

# survival regression
cv_survival_regression <-
  cv_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    time_col = time_to_event,
    event_col = event,
    model_type = "survival",
    hyperparameter_grid = reg_grid_3
  )

bootstrap_survival_regression <-
  boot_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    time_col = time_to_event,
    event_col = event,
    model_type = "survival",
    hyperparameter_grid = reg_grid_3
  )

simple_survival_regression <-
  simple_split %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    time_col = time_to_event,
    event_col = event,
    model_type = "survival",
    hyperparameter_grid = reg_grid_3
  )

unsplit_survival_regression <-
  feature_tibble %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    time_col = time_to_event,
    event_col = event,
    model_type = "survival",
    hyperparameter_grid = reg_grid_3
  )


test_that("survival tof_model's have the right components", {
  components <-
    c("model", "recipe", "mixture", "penalty",
      "model_type", "outcome_colnames", "training_data",
      "tuning_metrics", "log_rank_thresholds", "best_log_rank_threshold")

  # names are correct
  expect_equal(setdiff(components, names(cv_survival_regression)), character())
  expect_equal(setdiff(components, names(bootstrap_survival_regression)), character())
  expect_equal(setdiff(components, names(simple_survival_regression)), character())
  expect_equal(setdiff(components, names(unsplit_survival_regression)), character())

  ### types/classes are correct
  survival_models <- list(cv_survival_regression, bootstrap_survival_regression,
                          simple_survival_regression, unsplit_survival_regression)

  models <- purrr::map_lgl(survival_models, ~ inherits(.x$model, "glmnet"))
  recipes <- purrr::map_lgl(survival_models, ~ inherits(.x$recipe, "recipe"))
  mixtures <- purrr::map_lgl(survival_models, ~ inherits(.x$mixture, "numeric"))
  penalties <- purrr::map_lgl(survival_models, ~ inherits(.x$penalty, "numeric"))
  model_types <-
    purrr::map_lgl(survival_models, ~ inherits(.x$model_type, "character"))
  outcome_colnames <-
    purrr::map_lgl(survival_models, ~ inherits(.x$outcome_colnames, "character"))
  training_datas <-
    purrr::map_lgl(survival_models, ~ inherits(.x$training_data, "tbl_df"))
  tuning_metrics <-
    purrr::map_lgl(survival_models, ~ inherits(.x$training_data, "tbl_df"))
  log_rank_thresholds <-
    purrr::map_lgl(survival_models, ~ inherits(.x$log_rank_thresholds, "tbl_df"))
  best_log_rank_thresholds <-
    purrr::map_lgl(survival_models, ~ inherits(.x$best_log_rank_threshold, "numeric"))

  expect_true(all(models))
  expect_true(all(recipes))
  expect_true(all(mixtures))
  expect_true(all(penalties))
  expect_true(all(model_types))
  expect_true(all(outcome_colnames))
  expect_true(all(training_datas))
  expect_true(all(tuning_metrics))
  expect_true(all(log_rank_thresholds))
  expect_true(all(best_log_rank_thresholds))
})

test_that("survival tof_model's training data has the right shape", {

  # number of rows
  expect_equal(
    nrow(feature_tibble),
    nrow(cv_survival_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(bootstrap_survival_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(simple_survival_regression$training_data)
  )
  expect_equal(
    nrow(feature_tibble),
    nrow(unsplit_survival_regression$training_data)
  )

  # columns
  expect_equal(
    colnames(feature_tibble),
    colnames(cv_survival_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(bootstrap_survival_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(simple_survival_regression$training_data)
  )
  expect_equal(
    colnames(feature_tibble),
    colnames(unsplit_survival_regression$training_data)
  )
})

test_that("survival tof_model's tuning_metrics has the right shape", {
  survival_tuning_metrics <-
    c("mixture", "penalty", "neg_log_partial_likelihood", "concordance_index")

  # columns
  expect_equal(
    c("id", survival_tuning_metrics, ".predictions"),
    colnames(cv_survival_regression$tuning_metrics)
  )
  expect_equal(
    c("id", survival_tuning_metrics, ".predictions"),
    colnames(bootstrap_survival_regression$tuning_metrics)
  )
  expect_equal(
    survival_tuning_metrics,
    colnames(simple_survival_regression$tuning_metrics)
  )
  expect_equal(
    survival_tuning_metrics,
    colnames(unsplit_survival_regression$tuning_metrics)
  )
})

# tof_predict ------------------------------------------------------------------

# linear regression
linear_predictions_0 <-
  bootstrap_linear_regression %>%
  tof_predict(
    new_data = dplyr::arrange(bootstrap_linear_regression$training_data, sample)
  )

linear_predictions_1 <-
  bootstrap_linear_regression %>%
  tof_predict(new_data = dplyr::arrange(feature_tibble, sample))

linear_predictions_2 <-
  bootstrap_linear_regression %>%
  tof_predict(new_data = feature_tibble[1:10,])

linear_predictions_3 <-
  bootstrap_linear_regression %>%
  tof_predict(new_data = feature_tibble[1:10,], prediction_type = "link")

test_that("linear regression predictions have the correct shape", {
  # number of rows
  expect_equal(nrow(linear_predictions_0), 100L)
  expect_equal(nrow(linear_predictions_1), 100L)
  expect_equal(nrow(linear_predictions_2), 10L)
  expect_equal(nrow(linear_predictions_3), 10L)

  # number of columns
  expect_equal(ncol(linear_predictions_0), 1L)
  expect_equal(ncol(linear_predictions_1), 1L)
  expect_equal(ncol(linear_predictions_2), 1L)
  expect_equal(ncol(linear_predictions_3), 1L)
})

test_that("linear regression predictions are identical where expected", {
  # identical data
  expect_identical(linear_predictions_0, linear_predictions_1)
  expect_identical(linear_predictions_2, linear_predictions_3)
})

test_that("linear regression types are correct", {
  expect_type(linear_predictions_1$.pred, "double")
  expect_type(linear_predictions_2$.pred, "double")
  expect_type(linear_predictions_3$.pred, "double")
})

# logistic regression
logistic_predictions_0 <-
  bootstrap_logistic_regression %>%
  tof_predict(
    new_data = dplyr::arrange(bootstrap_logistic_regression$training_data, sample)
  )

logistic_predictions_1 <-
  bootstrap_logistic_regression %>%
  tof_predict(new_data = dplyr::arrange(feature_tibble, sample))

logistic_predictions_2 <-
  bootstrap_logistic_regression %>%
  tof_predict(new_data = feature_tibble[1:10,])

logistic_predictions_3 <-
  bootstrap_logistic_regression %>%
  tof_predict(new_data = feature_tibble[1:10,], prediction_type = "class")

test_that("logistic regression predictions have the correct shape", {
  # number of rows
  expect_equal(nrow(logistic_predictions_0), 100L)
  expect_equal(nrow(logistic_predictions_1), 100L)
  expect_equal(nrow(logistic_predictions_2), 10L)
  expect_equal(nrow(logistic_predictions_3), 10L)

  # number of columns
  expect_equal(ncol(logistic_predictions_0), 1L)
  expect_equal(ncol(logistic_predictions_1), 1L)
  expect_equal(ncol(logistic_predictions_2), 1L)
  expect_equal(ncol(logistic_predictions_3), 1L)
})

test_that("logistic regression predictions are identical where expected", {
  # identical data
  expect_identical(logistic_predictions_0, logistic_predictions_1)
})

test_that("logistic regression types are correct", {
  expect_type(logistic_predictions_1$.pred, "double")
  expect_type(logistic_predictions_2$.pred, "double")
  expect_type(logistic_predictions_3$.pred, "character")
})

# multinomial regression
multinomial_predictions_0 <-
  bootstrap_multinomial_regression %>%
  tof_predict(dplyr::arrange(bootstrap_multinomial_regression$training_data, sample))

multinomial_predictions_1 <-
  bootstrap_multinomial_regression %>%
  tof_predict(arrange(feature_tibble, sample))

multinomial_predictions_2 <-
  bootstrap_multinomial_regression %>%
  tof_predict(new_data = feature_tibble[1:10,])

multinomial_predictions_3 <-
  bootstrap_multinomial_regression %>%
  tof_predict(new_data = feature_tibble[1:10,], prediction_type = "class")

multinomial_predictions_4 <-
  bootstrap_multinomial_regression %>%
  tof_predict()

test_that("multinomial regression predictions have the correct shape", {
  # number of rows
  expect_equal(nrow(multinomial_predictions_0), 100L)
  expect_equal(nrow(multinomial_predictions_1), 100L)
  expect_equal(nrow(multinomial_predictions_2), 10L)
  expect_equal(nrow(multinomial_predictions_3), 10L)
  expect_equal(nrow(multinomial_predictions_4), 100L)

  # number of columns
  expect_equal(ncol(multinomial_predictions_0), 3L)
  expect_equal(ncol(multinomial_predictions_1), 3L)
  expect_equal(ncol(multinomial_predictions_2), 3L)
  expect_equal(ncol(multinomial_predictions_3), 1L)
  expect_equal(ncol(multinomial_predictions_4), 3L)
})

test_that("multinomial regression predictions are identical where expected", {
  # identical data
  expect_identical(multinomial_predictions_0, multinomial_predictions_1)
})

test_that("multinomial regression types are correct", {
  expect_type(multinomial_predictions_1[[1]], "double")
  expect_type(multinomial_predictions_1[[2]], "double")
  expect_type(multinomial_predictions_1[[3]], "double")

  expect_type(multinomial_predictions_2[[1]], "double")
  expect_type(multinomial_predictions_2[[2]], "double")
  expect_type(multinomial_predictions_2[[3]], "double")

  expect_type(multinomial_predictions_3$.pred, "character")

  expect_type(multinomial_predictions_4[[1]], "double")
  expect_type(multinomial_predictions_4[[2]], "double")
  expect_type(multinomial_predictions_4[[3]], "double")
})

# survival regression
survival_predictions_0 <-
  bootstrap_survival_regression %>%
  tof_predict(dplyr::arrange(bootstrap_survival_regression$training_data, sample))

survival_predictions_1 <-
  bootstrap_survival_regression %>%
  tof_predict(arrange(feature_tibble, sample))

survival_predictions_2 <-
  bootstrap_survival_regression %>%
  tof_predict(new_data = feature_tibble[1:10,])

survival_predictions_3 <-
  bootstrap_survival_regression %>%
  tof_predict(new_data = feature_tibble[1:10,], prediction_type = "link")

survival_predictions_4 <-
  bootstrap_survival_regression %>%
  tof_predict(
    new_data = feature_tibble[1:10,],
    prediction_type = "survival curve"
  )

survival_predictions_5 <-
  bootstrap_survival_regression %>%
  tof_predict()

test_that("survival regression predictions have the correct shape", {
  # number of rows
  expect_equal(nrow(survival_predictions_0), 100L)
  expect_equal(nrow(survival_predictions_1), 100L)
  expect_equal(nrow(survival_predictions_2), 10L)
  expect_equal(nrow(survival_predictions_3), 10L)
  expect_equal(nrow(survival_predictions_4), 10L)
  expect_equal(nrow(survival_predictions_5), 100)

  # number of columns
  expect_equal(ncol(survival_predictions_0), 1L)
  expect_equal(ncol(survival_predictions_1), 1L)
  expect_equal(ncol(survival_predictions_2), 1L)
  expect_equal(ncol(survival_predictions_3), 1L)
})

test_that("survival regression predictions are identical where expected", {
  # identical data
  expect_identical(survival_predictions_0, survival_predictions_1)
})

test_that("survival regression types are correct", {
  expect_type(survival_predictions_1$.pred, "double")
  expect_type(survival_predictions_2$.pred, "double")
  expect_type(survival_predictions_3$.pred, "double")
  expect_type(survival_predictions_4$.survival_curve, "list")
  expect_type(survival_predictions_5$.pred, "double")
  expect_type(survival_predictions_0$.pred, "double")
})


# tof_assess_model -------------------------------------------------------------

# linear regression
linear_a_0 <-
  bootstrap_linear_regression %>%
  tof_assess_model()

linear_a_1 <-
  bootstrap_linear_regression %>%
  tof_assess_model(new_data = feature_tibble)

linear_a_2 <-
  bootstrap_linear_regression %>%
  tof_assess_model(new_data = feature_tibble[1:10,])

# add a random variable that wasn't in the original dataset
linear_a_3 <-
  bootstrap_linear_regression %>%
  tof_assess_model(new_data = mutate(feature_tibble, new_var = "a"))

linear_a_tuning <-
  bootstrap_linear_regression %>%
  tof_assess_model(new_data = "tuning")

linear_assessments <-
  list(linear_a_0, linear_a_1, linear_a_2, linear_a_3, linear_a_tuning)

test_that("linear assessment results give a single model_metrics table", {

  # list component is model_metrics and nothing else
  has_model_metrics <-
    purrr::map(linear_assessments, ~ setdiff(names(.x), "model_metrics")) %>%
    purrr::map_lgl(.f = ~ identical(.x, character()))
  expect_true(all(has_model_metrics))

  # all model_metrics inherit the tbl_df class
  is_tbl_df <-
    purrr::map_lgl(linear_assessments, ~inherits(.x$model_metrics, what = "tbl_df"))
  expect_true(all(is_tbl_df))

  # contents of the model_metrics table are correct
  metric_cols_correct <-
    linear_assessments %>%
    map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("metric", "value"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(metric_cols_correct))

  metric_row_num <-
    linear_assessments %>%
    purrr::map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map_int(nrow)

  expect_true(all(metric_row_num == 2L))

})

# logistic regression
logistic_a_1 <-
  bootstrap_logistic_regression %>%
  tof_assess_model(new_data = feature_tibble)

# add a random variable that wasn't in the original dataset
logistic_a_2 <-
  bootstrap_logistic_regression %>%
  tof_assess_model(new_data = mutate(feature_tibble, new_var = "a"))

logistic_a_tuning <-
  bootstrap_logistic_regression %>%
  tof_assess_model(new_data = "tuning")

logistic_assessments <-
  list(logistic_a_1, logistic_a_2, logistic_a_tuning)

test_that("logistic assessment results give a model_metrics table, an roc_curve, and a confusion matrix", {

  # list component is model_metrics and nothing else
  has_right_components <-
    purrr::map(
      logistic_assessments,
      ~ setdiff(names(.x), c("model_metrics", "roc_curve", "confusion_matrix"))
    ) %>%
    purrr::map_lgl(.f = ~ identical(.x, character()))
  expect_true(all(has_right_components))

  # all model_metrics inherit the tbl_df class
  is_tbl_df <-
    purrr::map_lgl(
      logistic_assessments,
      ~inherits(.x$model_metrics, what = "tbl_df")
    )
  expect_true(all(is_tbl_df))

  # contents of the model_metrics table are correct
  metric_cols_correct <-
    logistic_assessments %>%
    map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("metric", "value"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(metric_cols_correct))

  metric_row_num <-
    logistic_assessments %>%
    purrr::map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map_int(nrow)

  expect_true(all(metric_row_num == 6L))

  # contents of the roc_curve table are correct
  roc_cols_correct <-
    logistic_assessments %>%
    map(~purrr::pluck(.x, "roc_curve")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c(".threshold", "specificity", "sensitivity", "tpr", "fpr"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(roc_cols_correct))

  # contents of the confusion_matrix table are correct
  confusion_cols_correct <-
    logistic_assessments %>%
    purrr::map(~purrr::pluck(.x, "confusion_matrix")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("true_outcome", "predicted_outcome", "num_observations"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(confusion_cols_correct))

})

test_that("Warning if only 1 class is present in new_data", {
  # expect a warning
    expect_warning(
      bootstrap_logistic_regression %>%
        tof_assess_model(new_data = filter(feature_tibble, class == "class1"))
    )
})

# multinomial regression
multinomial_a_1 <-
  bootstrap_multinomial_regression %>%
  tof_assess_model(new_data = feature_tibble)


# add a random variable that wasn't in the original dataset
multinomial_a_2 <-
  bootstrap_multinomial_regression %>%
  tof_assess_model(new_data = mutate(feature_tibble, new_var = "a"))

multinomial_a_tuning <-
  bootstrap_multinomial_regression %>%
  tof_assess_model(new_data = "tuning")

multinomial_assessments <-
  list(multinomial_a_1, multinomial_a_2, multinomial_a_tuning)


test_that("multinomial assessment results give a model_metrics table, an roc_curve, and a confusion matrix", {

  # list component is model_metrics and nothing else
  has_right_components <-
    purrr::map(
      multinomial_assessments,
      ~ setdiff(names(.x), c("model_metrics", "roc_curve", "confusion_matrix"))
    ) %>%
    purrr::map_lgl(.f = ~ identical(.x, character()))
  expect_true(all(has_right_components))

  # all model_metrics inherit the tbl_df class
  is_tbl_df <-
    purrr::map_lgl(
      multinomial_assessments,
      ~inherits(.x$model_metrics, what = "tbl_df")
    )
  expect_true(all(is_tbl_df))

  # contents of the model_metrics table are correct
  metric_cols_correct <-
    multinomial_assessments %>%
    map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("metric", "value"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(metric_cols_correct))

  metric_row_num <-
    multinomial_assessments %>%
    purrr::map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map_int(nrow)

  expect_true(all(metric_row_num == 6L))

  # contents of the roc_curve table are correct
  roc_cols_correct <-
    multinomial_assessments %>%
    map(~purrr::pluck(.x, "roc_curve")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c(".level", ".threshold", "specificity", "sensitivity", "tpr", "fpr"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(roc_cols_correct))

  # contents of the confusion_matrix table are correct
  confusion_cols_correct <-
    multinomial_assessments %>%
    purrr::map(~purrr::pluck(.x, "confusion_matrix")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("true_outcome", "predicted_outcome", "num_observations"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(confusion_cols_correct))

})

test_that("Warning if only 1 class is present in new_data", {
  expect_warning(
    expect_warning(
      bootstrap_multinomial_regression %>%
        tof_assess_model(new_data = feature_tibble[1:10,])
    )
  )
})

# survival regression
survival_a_1 <-
  bootstrap_survival_regression %>%
  tof_assess_model(new_data = feature_tibble)

survival_a_2 <-
  bootstrap_survival_regression %>%
  tof_assess_model(new_data = feature_tibble[1:60,])

# add a random variable that wasn't in the original dataset
survival_a_3 <-
  bootstrap_survival_regression %>%
  tof_assess_model(new_data = mutate(feature_tibble, new_var = "a"))

survival_a_tuning <-
  bootstrap_survival_regression %>%
  tof_assess_model(new_data = "tuning")

survival_assessments <-
  list(survival_a_1, survival_a_2, survival_a_3, survival_a_tuning)

test_that("survival assessment results give a model_metrics table and survival curves", {

  # list component is model_metrics and nothing else
  has_right_components <-
    purrr::map(
      survival_assessments,
      ~ setdiff(names(.x), c("model_metrics", "survival_curves"))
    ) %>%
    purrr::map_lgl(.f = ~ identical(.x, character()))
  expect_true(all(has_right_components))

  # all model_metrics inherit the tbl_df class
  is_tbl_df <-
    purrr::map_lgl(
      survival_assessments,
      ~inherits(.x$model_metrics, what = "tbl_df")
    )
  expect_true(all(is_tbl_df))

  # contents of the model_metrics table are correct
  metric_cols_correct <-
    survival_assessments %>%
    map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map(colnames) %>%
    purrr::map(~setdiff(.x, c("metric", "value"))) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(metric_cols_correct))

  metric_row_num <-
    survival_assessments %>%
    purrr::map(~purrr::pluck(.x, "model_metrics")) %>%
    purrr::map_int(nrow)

  expect_true(all(metric_row_num == 3L))

  # contents of the survival_curves table are correct
  survival_cols_correct <-
    survival_assessments[-length(survival_assessments)] %>%
    map(~purrr::pluck(.x, "survival_curves")) %>%
    purrr::map(colnames) %>%
    purrr::map(
      ~setdiff(
        .x,
        c("row_index", "survival_curve", "risk_group", "relative_risk", "time_to_event", "event")
      )
    ) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(survival_cols_correct))

  tuning_survival_cols_correct <-
    survival_assessments[length(survival_assessments)] %>%
    map(~purrr::pluck(.x, "survival_curves")) %>%
    purrr::map(colnames) %>%
    purrr::map(
      ~setdiff(
        .x,
        c("survival_curve", "risk_group", "relative_risk", "time_to_event", "event")
      )
    ) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(tuning_survival_cols_correct))

  # contents of the survival_curves for each sample are correct
  survival_curves_correct <-
    survival_assessments %>%
    map(~purrr::pluck(.x, "survival_curves")) %>%
    purrr::map(~tidyr::unnest(.x, cols = "survival_curve")) %>%
    purrr::map(colnames) %>%
    purrr::map(
      ~setdiff(
        .x,
        c("row_index", "probability", "time", "risk_group", "relative_risk", "time_to_event", "event")
      )
    ) %>%
    purrr::map_lgl(~identical(.x, character()))
  expect_true(all(survival_curves_correct))

})
