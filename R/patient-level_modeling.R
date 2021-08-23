# patient-level_modeling.R
# This file contains functions relevant to performing patient- or sample-level
# predictive modeling on aggregated single-cell CyTOF data.

# data splitting  ---------------

tof_split_data <-
  function(
    feature_tibble,
    split_method = c("k-fold", "bootstrap", "simple"),
    simple_prop = 3/4,
    num_cv_folds = 10,
    num_cv_repeats = 1,
    num_bootstraps = 25,
    strata = NULL,
    ... # optional additional arguments to the rsample function
  ) {
    # check split method
    split_method = rlang::arg_match(split_method)

    # perform split
    if (split_method == "simple") {
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


# building models -------------------------------------

#' Title
#'
#' @param feature_tibble TO DO
#'
#' @param response_col TO DO
#'
#' @param predictor_cols TO DO
#'
#' @param split_method TO DO
#'
#' @param ... TO DO
#'
#' @param normalize_predictors TO DO
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
tof_train_regression <-
  function(
    feature_tibble,
    response_col,
    predictor_cols = -{{response_col}}, #check that this doesn't cause problems
    split_method = c("k-fold", "simple", "bootstrap"),
    ..., # options for the data splitting function
    normalize_predictors = TRUE,
    remove_zv_predictors = TRUE,
    impute_missing_predictors = FALSE,
    grid_type = c("entropy", "regular"),
    grid_size = 100,
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
        normalize_predictors = normalize_predictors,
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
#' @param feature_tibble TO DO
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
#' @param normalize_predictors TO DO
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
    split_method = c("k-fold", "simple", "bootstrap"),
    ..., # options for the data splitting function
    model_type = c("two-class", "multiclass"),
    normalize_predictors = TRUE,
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
        normalize_predictors = normalize_predictors,
        remove_zv_predictors = remove_zv_predictors,
        impute_missing_predictors = impute_missing_predictors,
        grid_type = grid_type,
        grid_size = grid_size,
        model_type = model_type
      )

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
#' @param feature_tibble TO DO
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
        normalize_predictors = FALSE,
        remove_zv_predictors = remove_zv_predictors,
        impute_missing_predictors = impute_missing_predictors,
        grid_type = grid_type,
        grid_size = grid_size, # number of alpha values to try, as glmnet provides the full lambda path
        model_type = "survival"
      )

    # pull out some values in the format glmnet wants them ---------------------

    # the ids of which fold each row of feature_tibble is "in", i.e. the
    # fold in which each row is in the assessment set
    fold_ids <-
      setup_list %>%
      purrr::pluck("split_data") %>%
      generics::tidy() %>%
      dplyr::filter(Data == "Assessment") %>%
      dplyr::arrange(Row) %>%
      dplyr::pull(Fold) %>%
      stringr::str_remove("Fold") %>%
      as.numeric()

    # the untrained recipe
    surv_recipe <- setup_list$model_recipe

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


# applying models to new data -------------------------------------


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

    # depending on the type of model, provide some performance metrics

  }


