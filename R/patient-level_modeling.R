# patient-level_modeling.R
# This file contains functions relevant to performing patient- or sample-level
# predictive modeling on aggregated single-cell CyTOF data.

# Model Setup  -----------------------------------------------------------------

#' Split CyTOF data into a training and test set
#'
#' @param feature_tibble A tibble in which each row represents a sample- or patient-
#' level observation, such as those produced by \code{tof_extract_features}.
#'
#' @param split_method Either a string or a logical vector specifying how to perform
#' the split. If a string, valid options include k-fold cross validation
#' ("k-fold"; the default), bootstrapping ("bootstrap"), or
#' a single binary split ("simple"). If a logical vector, it should contain one entry
#' for each row in `feature_tibble` indicating if that row should be included in the
#' training set (TRUE) or excluded for the validation/test set (FALSE).
#' Ignored if `split_col` is specified.
#'
#' @param split_col The unquoted column name of the logical column in `feature_tibble`
#' indicating if each row should be included in the
#' training set (TRUE) or excluded for the validation/test set (FALSE).
#'
#' @param simple_prop A numeric value between 0 and 1 indicating what proportion of the data
#' should be used for training. Defaults to 3/4. Ignored if split_method is not "simple".
#'
#' @param num_cv_folds An integer indicating how many cross-validation folds should be used.
#' Defaults to 10. Ignored if split_method is not "k-fold".
#'
#' @param num_cv_repeats An integer indicating how many independent
#' cross-validation replicates should be used (i.e. how many num_cv_fold splits
#' should be performed). Defaults to 1. Ignored if split_method is not "k-fold".
#'
#' @param num_bootstraps An integer indicating how many independent bootstrap
#' replicates should be used. Defaults to 25. Ignored if split_method is not
#' "bootstrap".
#'
#' @param strata An unquoted column name representing the column in \code{feature_tibble}
#' that should be used to stratify the data splitting. Defaults to NULL (no stratification).
#'
#' @param ... Optional additional arguments to pass to \code{\link[rsample]{vfold_cv}}
#' for k-fold cross validation, \code{\link[rsample]{bootstraps}} for bootstrapping,
#' or \code{\link[rsample]{initial_split}} for simple splitting.
#'
#' @return If for k-fold cross validation and bootstrapping, an "rset" object;
#' for simple splitting, an "rsplit" object. For details, see
#' \code{\link[rsample]{rsample}}.
#'
#' @export
#'
tof_split_data <-
  function(
    feature_tibble,
    split_method = c("k-fold", "bootstrap", "simple"),
    split_col,
    simple_prop = 3/4,
    num_cv_folds = 10,
    num_cv_repeats = 1,
    num_bootstraps = 10,
    strata = NULL,
    ...
  ) {

    # check split_col
    if (!missing(split_col)) {
      split_method <-
        feature_tibble %>%
        dplyr::pull({{split_col}})

      if (!is.logical(split_method)) {
        stop("If `split_col` points to a column in `feature_tibble`, all values
               in the column must be logical.")
      }
    }

    # check split method
    if (is.character(split_method)) {
      split_method <- rlang::arg_match(split_method)
    } else if (is.logical(split_method)) {
      rows_match <- nrow(feature_tibble) == length(split_method)
      if (!rows_match) {
        stop("If `split_method` is a logical vector, it must have one entry for each
             row of `feature_tibble.`")
      }
    } else {
      stop("`split_method` must be either a valid string or a logical vector.")
    }

    # perform split
    if (is.logical(split_method)) {
      training_ids <- which(split_method)
      validation_ids <- which(!split_method)
      split_result <-
        list(data = feature_tibble, in_id = training_ids, out_id = NA, id = "Resample1")
      class(split_result) <- "rsplit"

    } else if (split_method == "simple") {
      split_result <-
        feature_tibble %>%
        rsample::initial_split(prop = simple_prop, strata = {{strata}}, ...)

    } else if (split_method == "k-fold") {
      split_result <-
        feature_tibble %>%
        rsample::vfold_cv(
          v = num_cv_folds,
          repeats = num_cv_repeats,
          strata = {{strata}},
          ...
        )

    } else {
      split_result <-
        feature_tibble %>%
        rsample::bootstraps(times = num_bootstraps, strata = {{strata}}, ...)
    }

    # return result
    return(split_result)
  }


#' Create a {glmnet} hyperparameter search grid of a specified size
#'
#' This function creates a hyperparameter search grid (in the form of a
#' \code{\link[tibble]{tibble}}) specifying the search space for the two
#' hyperparameters of a generalized linear model using the {glmnet} package:
#' the regularization penalty term
#' and the lasso/ridge regression mixture term. Either a regular or maximum-entropy
#' search grid (of a specified size) or a grid with a user-specified range of
#' hyperparameter values can be created.
#'
#'
#' @param penalty_values A numeric vector of all penalty values to include in the
#' hyperparamter grid. If unspecified, {tidytof}'s default values will be used:
#' 10^(-10), 10^(-7.5), 10^(-5), 10^(-2.5), 10^(0).
#'
#' @param mixture_values  A numeric vector of all penalty values to include in the
#' hyperparamter grid. If unspecified, {tidytof}'s default values will be used:
#' `seq(0, 1, length.out = 10)`.
#'
#' @param grid_type A string indicating if a regular grid ("regular"; the default) or maximum-entropy
#' ("entropy") search grid should be returned. Ignored if both `penalty_values` and
#' `mixture_values` are provided. If "entropy" is chosen, the {dials} package
#' is required.
#'
#' @param grid_size An integer indicating how many combinations of penalty and mixture
#' values should be searched during model tuning. Defaults to 25 (5 unique values
#' of `penalty` and 5 unique values of `mixture`).
#'
#' @return A tibble with two numeric columns: `penalty` and `mixture`.
#'
#' @export
#'
#' @importFrom rlang arg_match
#' @importFrom tidyr expand_grid
#'
#' @examples
#' tof_create_grid()
#'
tof_create_grid <-
  function(
    penalty_values,
    mixture_values,
    grid_type = c("regular", "entropy"),
    grid_size = 25
  ) {

    # check grid_type argument
    grid_type <- rlang::arg_match(grid_type)

    # create grid
    if (!missing(penalty_values) & !missing(mixture_values)) {
      hyperparam_grid <-
        tidyr::expand_grid(penalty_values, mixture_values) %>%
        dplyr::rename(
          penalty = penalty_values,
          mixture = mixture_values
        )

    } else if (grid_type == "regular") {
      hyperparam_grid <-
        tidyr::expand_grid(
          penalty = 10^(seq(-10, 0, length.out = round(sqrt(grid_size)))),
          mixture = seq(0, 1, length.out = round(sqrt(grid_size)))
        )

    } else if (grid_type == "entropy") {
      # conditionally load the dials package for an entropy grid
      has_dials <- requireNamespace(package = "dials")
      if (!has_dials) {
        stop(
          "Creating a maximum-entropy search grid requires the {dials} package. Install it with this code:\n
           install.packages(\"dials\")"
        )
      }

      hyperparams <-
        dials::parameters(list(dials::penalty(), dials::mixture()))

      hyperparam_grid <-
        dials::grid_max_entropy(hyperparams, size = grid_size)
    } else {
      stop("If both `penalty_values` and `mixture_values` are unspecified, the grid
           type must be either regular or entropy.")
    }

    return(hyperparam_grid)
  }



# Building Models --------------------------------------------------------------

tof_train_model <-
  function(
    split_data,
    predictor_cols,
    response_col, # classification and regression
    time_col, # survival
    event_col, # survival
    model_type = c("linear", "two-class", "multiclass", "survival"),
    hyperparameter_grid = tof_create_grid(),
    standardize_predictors = TRUE,
    remove_zv_predictors = FALSE,
    impute_missing_predictors = FALSE,
    optimization_metric = "tidytof_default",
    best_model_type = c("best", "best with sparsity"),
    num_cores = 1
  ) {

    # check arguments ----------------------------------------------------------

    # split_data
    if (inherits(split_data, "rsplit")) {
      feature_tibble <-
        split_data %>%
        rsample::training()
    } else if (inherits(split_data, "rset")) {
      feature_tibble <-
        split_data$splits[[1]] %>%
        rsample::training()
    } else {
      stop("split_data must be either an rsplit or an rset object.")
    }

    # string arguments
    model_type <- rlang::arg_match(model_type)
    best_model_type <- rlang::arg_match(best_model_type)

    # outcome variables
    if (model_type %in% c("linear", "two-class", "multiclass") & !missing(response_col)) {
      response <-
        feature_tibble %>%
        dplyr::pull({{response_col}})

      if (model_type == "linear") {
        if (!is.numeric(response)) {
          stop("`response_col` must specify a numeric column for linear regression.")
        }
      } else if (model_type %in% c("two-class", "multiclass")) {
        if (!is.factor(response)) {
          stop("`response_col` must specify a factor column for logistic and multinomial regression.")
        }
      }
    } else if (model_type == "survival" & !missing(time_col) & !missing(event_col)) {
      time <-
        feature_tibble %>%
        dplyr::pull({{time_col}})

      event <-
        feature_tibble %>%
        dplyr::pull({{event_col}})

      if (!is.numeric(time)) {
        stop("All values in time_col must be numeric - they represent the time
           to event outcome for each patient (i.e. each row of feature_tibble.")
      }

      if (!all(event %in% c(0, 1))) {
        stop("All values in event_col must be either 0 or 1 (with 1 indicating
           the adverse event) or FALSE and TRUE (with TRUE indicating the
           adverse event.")
      }

    } else {
      stop(
        "Either `response_col` (for linear, two-class, and multiclass models)
           or both `time_col` and `event_col` (for survival models) must be specified."
      )
    }

    # create and prep recipe ---------------------------------------------------

    # create recipe
    if (model_type %in% c("linear", "two-class", "multiclass")) {
      # if any type of model other than survival
      unprepped_recipe <-
        feature_tibble %>%
        tof_create_recipe(
          predictor_cols = {{predictor_cols}},
          outcome_cols = {{response_col}},
          standardize_predictors = standardize_predictors,
          remove_zv_predictors = remove_zv_predictors,
          impute_missing_predictors = impute_missing_predictors
        )
    } else {
      # if survival model
      unprepped_recipe <-
        feature_tibble %>%
        tof_create_recipe(
          predictor_cols = {{predictor_cols}},
          outcome_cols = c({{time_col}}, {{event_col}}), # see if there is a way to do this using tidyeval
          standardize_predictors = standardize_predictors,
          remove_zv_predictors = remove_zv_predictors,
          impute_missing_predictors = impute_missing_predictors
        )
    }

    # prep recipe(s) - will be a single recipe if split_data is a single rsplit,
    # but will be a list of recipes if split_data is an rset object.
    prepped_recipe <-
      split_data %>%
      tof_prep_recipe(unprepped_recipe = unprepped_recipe)

    # tune the model using the resampling, the recipes, and the
    # hyperparameter grid ------------------------------------------------------
    if (model_type == "survival") {
      tuning_metrics <-
        tof_tune_glmnet(
          split_data = split_data,
          prepped_recipe = prepped_recipe,
          hyperparameter_grid = hyperparameter_grid,
          model_type = model_type,
          outcome_cols = c({{time_col}}, {{event_col}}),
          optimization_metric = optimization_metric,
          num_cores = num_cores
        )

    } else {
      tuning_metrics <-
        tof_tune_glmnet(
          split_data = split_data,
          prepped_recipe = prepped_recipe,
          hyperparameter_grid = hyperparameter_grid,
          model_type = model_type,
          outcome_cols = {{response_col}},
          optimization_metric = optimization_metric,
          num_cores = num_cores
        )
    }

    # find the best model using the tuning metrics -----------------------------

    return(tuning_metrics)





    # pull out some values in the format glmnet wants them ---------------------

    # the ids of which fold each row of feature_tibble is "in", i.e. the
    # fold in which each row is in the assessment set
    # currently this will only work for cross-validation
    fold_ids <-
      setup_list %>%
      purrr::pluck("split_data") %>%
      rsample::tidy() %>%
      dplyr::filter(Data == "Assessment") %>%
      dplyr::arrange(Row) %>%
      dplyr::pull(Fold) %>%
      stringr::str_remove("Fold") %>%
      as.numeric()

    # the preprocessed feature_tibble (without standardization, which glmnet
    # will do for each fold internally)
    preprocessed_data <-
      setup_list %>%
      # there is some information leak here because all folds are being
      # used to impute and remove nzv predictors
      purrr::pluck("model_recipe") %>%
      recipes::prep() %>%
      recipes::juice()

    # vector containing all alpha values to try
    alphas <-
      setup_list %>%
      purrr::pluck("hyperparam_grid") %>%
      dplyr::pull(mixture)

    # vector containing all lambda values to try
    # lambdas <-
    #   setup_list %>%
    #   purrr::pluck("hyperparam_grid") %>%
    #   dplyr::pull(penalty) %>%
    #   sort(decreasing = TRUE)

    # predictors and survival outcome objects in the form glmnet needs them
    x <-
      preprocessed_data %>%
      dplyr::select(-{{time_col}}, -{{event_col}}) %>%
      as.matrix()

    y <-
      survival::Surv(
        time = time,
        event = status
      )

    # fit full lambda path for each alpha value
    if (run_parallel) {
      num_cores <- min(parallel::detectCores() - 1, nrow(setup_list$split_data))
      doParallel::registerDoParallel(num_cores)
    }

    surv_models <-
      purrr::map(
        .x = alphas,
        .f = ~ glmnet::cv.glmnet(
          x = x,
          y = y,
          alpha = .x,
          family = "cox",
          standardize = TRUE
        )
      )

    if (run_parallel) {
      foreach::registerDoSEQ()
    }

    tuning_results <-
      tibble(
        mixture = alphas,
        models = surv_models
      )

    # extract some performance metrics from each model

    tuning_perf_metrics <-
      tuning_results %>%
      dplyr::mutate(
        .metrics = purrr::map(.x = models, .f = broom::tidy)
      ) %>%
      dplyr::select(mixture, .metrics) %>%
      tidyr::unnest(cols = .metrics) %>%
      dplyr::rename_with(.fn = ~stringr::str_replace(.x, "\\.", "_")) %>%
      dplyr::mutate(
        .metric = "partial likelihood deviance",
        n = nrow(setup_list$split_data)
      ) %>%
      dplyr::rename(
        penalty = lambda,
        mean = estimate,
        num_nonzero_coefficients = nzero,
        std_err = std_error
      ) %>%
      dplyr::relocate(
        penalty,
        mixture,
        .metric,
        mean,
        n,
        std_err
      )

    # find best model overall
    best_model_info <-
      tuning_perf_metrics %>%
      dplyr::slice_min(order_by = mean, with_ties = FALSE, n = 1) %>%
      dplyr::select(penalty, mixture)

    best_alpha <- best_model_info$mixture
    best_lambda <- best_model_info$penalty

    if (best_model_type == "best with sparsity") {
      sparse_lambda <-
        tuning_results %>%
        dplyr::filter(mixture == best_alpha) %>%
        dplyr::mutate(
          best_penalty = purrr::map_dbl(.x = models, .f = ~.x$lambda.1se),
        ) %>%
        dplyr::pull(best_penalty)

      best_model_info$penalty <- sparse_lambda

    }

    # fit model with best hyperparameters on all data
    final_model <-
      glmnet::glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = best_model_info$mixture,
        lambda = best_model_info$penalty
      )

    # make predictions with final model
    final_mod_predictions <-
      final_model %>%
      stats::predict(newx = x, s = best_model_info$penalty)

    # find final model performance metrics
    final_mod_deviance <- final_model$dev.ratio
    final_mod_perf_metrics <-
      tibble::tibble(
        .metric = "deviance ratio",
        .estimator = "glmnet",
        .estimate = final_mod_deviance
      )

    # extract coefficients from final model
    final_coefficients <-
      final_model %>%
      stats::coef() %>%
      as.matrix() %>%
      tibble::as_tibble(rownames = "feature_name") %>%
      dplyr::rename(beta = s0) %>%
      dplyr::filter(abs(beta) > 0)

    results <-
      list(
        tuning_results = tuning_results,
        tuning_perf_metrics = tuning_perf_metrics,
        final_mod_predictions = final_mod_predictions,
        final_mod_hyperparams = best_model_info,
        final_mod_perf_metrics = final_mod_perf_metrics,
        final_model = final_model,
        final_coefficients = final_coefficients
      )

    return(results)
  }


#' Build a regression model using aggregated CyTOF data.
#'
#' @param feature_tibble A tibble in which each row represents a sample- or patient-
#' level observation, such as those produced by \code{tof_extract_features}.
#'
#' @param response_col An unquoted column name representing the column in \code{feature_tibble}
#' that that contains the (continuous) outcome variable being predicted by the
#' regression model.
#'
#' @param predictor_cols An unquoted column name representing the column in \code{feature_tibble}
#' that that contains the (continuous) outcome variable being predicted by the
#' regression model.
#'
#' @param hyperparameter_grid
#'
#' @param split_method TO DO
#'
#' @param ... TO DO
#'
#' @param standardize_predictors TO DO
#'
#' @param remove_zv_predictors TO DO
#'
#' @param impute_missing_predictors TO DO
#'
#' @param run_parallel TO DO
#'
#'
#' @return TO DO
#'
#' @export
#'
#'
tof_train_regression <-
  function(
    feature_tibble,
    response_col,
    predictor_cols,
    hyperparameter_grid = tof_create_grid(),
    split_method = c("k-fold", "simple", "bootstrap"),
    ..., # options for the data splitting function
    standardize_predictors = TRUE,
    remove_zv_predictors = TRUE,
    impute_missing_predictors = FALSE,
    run_parallel = FALSE
  ) {

    # perform all the setup for a glmnet model from any family
    setup_list <-
      tof_setup_glmnet_mod(
        feature_tibble = feature_tibble,
        response_col = {{response_col}},
        predictor_cols = {{predictor_cols}},
        split_method = split_method,
        ...,
        standardize_predictors = standardize_predictors,
        remove_zv_predictors = remove_zv_predictors,
        impute_missing_predictors = impute_missing_predictors,
        grid_type = grid_type,
        grid_size = grid_size,
        model_type = "linear"
      )

    # perform model tuning, feature selection, and final model fitting
    result <-
      tof_fit_glmnet_mod(
        feature_tibble = feature_tibble,
        response_col = {{response_col}},
        setup_list = setup_list,
        model_type = "linear",
        run_parallel = run_parallel
      )

    return(result)

  }

#' Title
#'
#' @param feature_tibble A tibble in which each row represents a sample- or patient-
#' level observation, such as those produced by \code{tof_extract_features}.
#'
#' @param response_col TO DO
#'
#' @param predictor_cols TO DO
#'
#' @param split_method TO DO
#'
#' @param ... TO DO
#'
#' @param model_type TO DO
#'
#' @param standardize_predictors TO DO
#'
#' @param remove_zv_predictors TO DO
#'
#' @param impute_missing_predictors TO DO
#'
#' @param grid_type TO DO
#'
#' @param grid_size TO DO
#'
#' @param run_parallel TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_train_classifier <-
  function(
    feature_tibble,
    response_col,
    predictor_cols = -{{response_col}}, #check that this doesn't cause problems
    split_method = c("k-fold", "simple", "bootstrap"), # simple doesn't work
    ..., # options for the data splitting function
    model_type = c("two-class", "multiclass"),
    standardize_predictors = TRUE,
    remove_zv_predictors = TRUE,
    impute_missing_predictors = FALSE,
    grid_type = c("entropy", "regular"),
    grid_size = 100,
    run_parallel = FALSE
  ) {

    # check mode argument
    model_type <- rlang::arg_match(model_type)

    # perform all the setup for a glmnet model from any family
    setup_list <-
      tof_setup_glmnet_mod(
        feature_tibble = feature_tibble,
        response_col = {{response_col}},
        predictor_cols = {{predictor_cols}},
        split_method = split_method,
        ...,
        standardize_predictors = standardize_predictors,
        remove_zv_predictors = remove_zv_predictors,
        impute_missing_predictors = impute_missing_predictors,
        grid_type = grid_type,
        grid_size = grid_size,
        model_type = model_type
      )

    return(setup_list)

    # perform model tuning, feature selection, and final model fitting
    result <-
      tof_fit_glmnet_mod(
        feature_tibble = feature_tibble,
        response_col = {{response_col}},
        setup_list = setup_list,
        model_type = model_type,
        run_parallel = run_parallel
      )

    return(result)
  }




#' Title
#'
#' @param feature_tibble A tibble in which each row represents a sample- or patient-
#' level observation, such as those produced by \code{tof_extract_features}.
#'
#' @param time_col TO DO
#'
#' @param event_col TO DO
#'
#' @param predictor_cols TO DO
#'
#' @param split_method TO DO
#'
#' @param ... TO DO
#'
#' @param remove_zv_predictors TO DO
#'
#' @param impute_missing_predictors TO DO
#'
#' @param grid_type TO DO
#'
#' @param grid_size TO DO
#'
#' @param best_model_type TO DO
#'
#' @param run_parallel TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_train_survival <-
  function(
    feature_tibble,
    time_col,
    event_col,
    predictor_cols = -c({{time_col}}, {{event_col}}), #check that this doesn't cause problems
    split_method = c("k-fold", "simple", "bootstrap"),
    ..., # options for the data splitting function
    remove_zv_predictors = TRUE,
    impute_missing_predictors = FALSE,
    grid_type = c("entropy", "regular"),
    grid_size = 100,
    best_model_type = c("best", "best with sparsity"),
    run_parallel = FALSE
  ) {

    # check best model type
    best_model_type <- rlang::arg_match(best_model_type)

    # extract (and check) the two outcome columns
    time <-
      feature_tibble %>%
      dplyr::pull({{time_col}})

    if (!is.numeric(time)) {
      stop("All values in time_col must be numeric - they represent the time
           to event outcome for each patient (i.e. each row of feature_tibble.")
    }

    status <-
      feature_tibble %>%
      dplyr::pull({{event_col}})

    if (!all(status %in% c(0, 1))) {
      stop("All values in event_col must be either 0 and 1 (with 1 indicating
           the adverse event) or FALSE and TRUE (with TRUE indicating the
           adverse event.")
    }

    # perform all the setup for a glmnet model from any family -----------------
    setup_list <-
      tof_setup_glmnet_mod(
        feature_tibble = feature_tibble,
        response_col = c({{time_col}}, {{event_col}}),
        predictor_cols = {{predictor_cols}},
        split_method = split_method,
        ...,
        # glmnet will standardize all the predictors in each fold for us
        standardize_predictors = FALSE,
        remove_zv_predictors = remove_zv_predictors,
        impute_missing_predictors = impute_missing_predictors,
        grid_type = grid_type,
        grid_size = grid_size, # number of alpha values to try, as glmnet provides the full lambda path
        model_type = "survival"
      )

    # pull out some values in the format glmnet wants them ---------------------

    # the ids of which fold each row of feature_tibble is "in", i.e. the
    # fold in which each row is in the assessment set
    # currently this will only work for cross-validation
    fold_ids <-
      setup_list %>%
      purrr::pluck("split_data") %>%
      rsample::tidy() %>%
      dplyr::filter(Data == "Assessment") %>%
      dplyr::arrange(Row) %>%
      dplyr::pull(Fold) %>%
      stringr::str_remove("Fold") %>%
      as.numeric()

    # the preprocessed feature_tibble (without standardization, which glmnet
    # will do for each fold internally)
    preprocessed_data <-
      setup_list %>%
      # there is some information leak here because all folds are being
      # used to impute and remove nzv predictors
      purrr::pluck("model_recipe") %>%
      recipes::prep() %>%
      recipes::juice()

    # vector containing all alpha values to try
    alphas <-
      setup_list %>%
      purrr::pluck("hyperparam_grid") %>%
      dplyr::pull(mixture)

    # vector containing all lambda values to try
    # lambdas <-
    #   setup_list %>%
    #   purrr::pluck("hyperparam_grid") %>%
    #   dplyr::pull(penalty) %>%
    #   sort(decreasing = TRUE)

    # predictors and survival outcome objects in the form glmnet needs them
    x <-
      preprocessed_data %>%
      dplyr::select(-{{time_col}}, -{{event_col}}) %>%
      as.matrix()

    y <-
      survival::Surv(
        time = time,
        event = status
      )

    # fit full lambda path for each alpha value
    if (run_parallel) {
      num_cores <- min(parallel::detectCores() - 1, nrow(setup_list$split_data))
      doParallel::registerDoParallel(num_cores)
    }

    surv_models <-
      purrr::map(
        .x = alphas,
        .f = ~ glmnet::cv.glmnet(
          x = x,
          y = y,
          alpha = .x,
          family = "cox",
          standardize = TRUE
        )
      )

    if (run_parallel) {
      foreach::registerDoSEQ()
    }

    tuning_results <-
      tibble(
        mixture = alphas,
        models = surv_models
      )

    # extract some performance metrics from each model

    tuning_perf_metrics <-
      tuning_results %>%
      dplyr::mutate(
        .metrics = purrr::map(.x = models, .f = broom::tidy)
      ) %>%
      dplyr::select(mixture, .metrics) %>%
      tidyr::unnest(cols = .metrics) %>%
      dplyr::rename_with(.fn = ~stringr::str_replace(.x, "\\.", "_")) %>%
      dplyr::mutate(
        .metric = "partial likelihood deviance",
        n = nrow(setup_list$split_data)
      ) %>%
      dplyr::rename(
        penalty = lambda,
        mean = estimate,
        num_nonzero_coefficients = nzero,
        std_err = std_error
      ) %>%
      dplyr::relocate(
        penalty,
        mixture,
        .metric,
        mean,
        n,
        std_err
      )

    # find best model overall
    best_model_info <-
      tuning_perf_metrics %>%
      dplyr::slice_min(order_by = mean, with_ties = FALSE, n = 1) %>%
      dplyr::select(penalty, mixture)

    best_alpha <- best_model_info$mixture
    best_lambda <- best_model_info$penalty

    if (best_model_type == "best with sparsity") {
      sparse_lambda <-
        tuning_results %>%
        dplyr::filter(mixture == best_alpha) %>%
        dplyr::mutate(
          best_penalty = purrr::map_dbl(.x = models, .f = ~.x$lambda.1se),
        ) %>%
        dplyr::pull(best_penalty)

      best_model_info$penalty <- sparse_lambda

    }

    # fit model with best hyperparameters on all data
    final_model <-
      glmnet::glmnet(
        x = x,
        y = y,
        family = "cox",
        alpha = best_model_info$mixture,
        lambda = best_model_info$penalty
      )

    # make predictions with final model
    final_mod_predictions <-
      final_model %>%
      stats::predict(newx = x, s = best_model_info$penalty)

    # find final model performance metrics
    final_mod_deviance <- final_model$dev.ratio
    final_mod_perf_metrics <-
      tibble::tibble(
        .metric = "deviance ratio",
        .estimator = "glmnet",
        .estimate = final_mod_deviance
      )

    # extract coefficients from final model
    final_coefficients <-
      final_model %>%
      stats::coef() %>%
      as.matrix() %>%
      tibble::as_tibble(rownames = "feature_name") %>%
      dplyr::rename(beta = s0) %>%
      dplyr::filter(abs(beta) > 0)

    results <-
      list(
        tuning_results = tuning_results,
        tuning_perf_metrics = tuning_perf_metrics,
        final_mod_predictions = final_mod_predictions,
        final_mod_hyperparams = best_model_info,
        final_mod_perf_metrics = final_mod_perf_metrics,
        final_model = final_model,
        final_coefficients = final_coefficients
      )

    return(results)

  }


# applying models to new data --------------------------------------------------

tof_predict <-
  function(
    model,
    new_data,
    ...
  ) {
    # extract the model from the results list
    NULL
    stop("Prediction is not yet supported.")
    # check if the model is a parsnip model or a coxnet model
    NULL

    # make predictions

  }

# assessing model performance --------------------------------------------------

tof_assess <-
  function(
    model,
    predictions,
    ...
  ) {
    # extract the model from the results list
    NULL
    stop("Assessment is not yet supported.")
    # check if the model is a parsnip model or a coxnet model
    NULL

    # depending on the type of model, provide some performance metrics

  }
