# modeling_helpers.R
#
# This file contains helper functions for performing sample- and patient-level
# modeling on tibble containing aggregated single-cell data.



#' Title
#'
#' TO DO
#'
#' @param feature_tibble TO DO
#' @param predictor_cols TO DO
#' @param outcome_cols TO DO
#' @param standardize_predictors TO DO
#' @param remove_zv_predictors TO DO
#' @param impute_missing_predictors TO DO
#'
#' @importFrom rlang arg_match
#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#' @importFrom recipes recipe
#' @importFrom recipes step_dummy
#' @importFrom recipes step_nzv
#' @importFrom recipes step_impute_knn
#' @importFrom recipes all_numeric_predictors
#' @importFrom recipes all_predictors
#'
#' @export
#'
tof_create_recipe <-
  function(
    feature_tibble,
    predictor_cols,
    outcome_cols,
    standardize_predictors = TRUE,
    remove_zv_predictors = FALSE,
    impute_missing_predictors = FALSE
  ) {

    # extract predictor column names
    predictor_colnames <-
      rlang::enquo(predictor_cols) %>%
      tidyselect::eval_select(data = feature_tibble) %>%
      names()

    outcome_colnames <-
      rlang::enquo(outcome_cols) %>%
      tidyselect::eval_select(data = feature_tibble) %>%
      names()

    # build recipe
    roles <-
      c(
        rep("outcome", times = length(outcome_colnames)),
        rep("predictor", times = length(predictor_colnames))
      )

    model_recipe <-
      recipes::recipe(
        x = feature_tibble,
        vars = c(outcome_colnames, predictor_colnames),
        roles = roles
      ) %>%
      recipes::step_select(
        dplyr::any_of(outcome_colnames),
        dplyr::any_of(predictor_colnames)
      ) %>%
      recipes::step_dummy(recipes::all_nominal_predictors())

    if (remove_zv_predictors) {
      model_recipe <-
        model_recipe %>%
        recipes::step_nzv(recipes::all_predictors())
    }

    if (impute_missing_predictors) {
      model_recipe <-
        model_recipe %>%
        recipes::step_impute_knn(recipes::all_numeric_predictors())
    }

    if (standardize_predictors) {
      model_recipe <-
        model_recipe %>%
        recipes::step_normalize(recipes::all_numeric_predictors(), na_rm = TRUE)
    }

    return(model_recipe)

  }





#' Title
#'
#' @param split_data TO DO
#'
#' @param unprepped_recipe TO DO
#'
#' @return TO DO
#'
#' @importFrom purrr map
#' @importFrom recipes prepper
#' @importFrom recipes prep
#'
#'
#' @export
#'
tof_prep_recipe <-
  function(
    split_data,
    unprepped_recipe
  ) {

    if (inherits(split_data, "rset")) {
      result <-
        purrr::map(
          .x = split_data$splits,
          .f = recipes::prepper,
          recipe = unprepped_recipe,
          retain = FALSE
        )

    } else if (inherits(split_data, "rsplit")) {
      result <-
        unprepped_recipe %>%
        recipes::prep(retain = FALSE)

    } else {
      stop("split_data must be either an rset or an rsplit object")
    }

    return(result)

  }

#' Title
#'
#' TO DO
#'
#' @param split_data TO DO
#' @param prepped_recipe TO DO
#' @param hyperparameter_grid TO DO
#' @param model_type TO DO
#' @param outcome_cols TO DO
#' @param optimization_metric TO DO
#' @param num_cores TO DO
#'
#' @return TO DO
#'
#'
#'
tof_tune_glmnet <-
  function(
    split_data,
    prepped_recipe,
    hyperparameter_grid,
    model_type,
    outcome_cols,
    optimization_metric = "tidytof_default",
    num_cores = 1
  ) {

    if (inherits(split_data, "rsplit")) {
      # extract output column names
      outcome_colnames <-
        rlang::enquo(outcome_cols) %>%
        tidyselect::eval_select(data = rsample::training(split_data)) %>%
        names()

      performance_metrics <-
        tof_fit_split(split_data, prepped_recipe, hyperparameter_grid, model_type, outcome_colnames)

    } else if (inherits(split_data, "rset")) {

      # extract output column names
      outcome_colnames <-
        rlang::enquo(outcome_cols) %>%
        tidyselect::eval_select(data = rsample::training(split_data$splits[[1]])) %>%
        names()

      # initialize cluster
      my_cluster <- NULL

      if (num_cores == 1) {
        performance_metrics <-
          purrr::map2(
            .x = split_data$splits,
            .y = prepped_recipe,
            .f = ~tof_fit_split(
              rsplit_data = .x,
              prepped_recipe = .y,
              hyperparameter_grid = hyperparameter_grid,
              model_type = model_type,
              outcome_colnames = outcome_colnames
            )
          )

      } else {
        my_cluster <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(my_cluster)
        `%my_do%` <- foreach::`%dopar%`

        performance_metrics <-
          foreach::foreach(
            split_data = split_data$splits,
            prepped_recipe = prepped_recipe,
            .combine = list,
            .packages =
              c("dplyr", "purrr", "survival", "glmnet", "yardstick", "recipes"),
            .export =
              c("tof_fit_split", "tof_create_grid", "tof_clean_metric_names"),
            .multicombine = TRUE
          ) %my_do%
          tof_fit_split(
            rsplit_data = split_data,
            prepped_recipe = prepped_recipe,
            hyperparameter_grid = hyperparameter_grid,
            model_type = model_type,
            outcome_colnames = outcome_colnames
          )

        # stop cluster if it was set up
        if (!is.null(my_cluster)) {
          parallel::stopCluster(my_cluster)
        }
      }

      performance_metrics <-
        tibble::tibble(
          performance_metrics = performance_metrics
        ) %>%
        tidyr::unnest(cols = performance_metrics) %>%
        dplyr::group_by(mixture, penalty) %>%
        dplyr::summarize(
          dplyr::across(
            dplyr::everything(),
            .fns = mean,
            #.fns = list(mean = mean, sd = stats::sd),
            #.names = "{.col}__tidytof__{.fn}"
            )
        ) %>%
        dplyr::ungroup() #%>%
        # tidyr::pivot_longer(
        #   cols = c(-mixture, -penalty),
        #   names_to = "metric",
        #   values_to = "value"
        # ) %>%
        # tidyr::separate(metric, into = c("metric", "statistic"), sep = "__tidytof__")

    } else {
      stop("split_data must be an rsplit or rset object.")
    }

    # extract best hyperparameters
    best_model_parameters <-
      tof_find_best(performance_metrics, model_type, optimization_metric)

    # finalize model

    tof_model <- tof_finalize_model(best_model_parameters, split_data)

    return(tof_model)
  }

#' Fit a glmnet model and calculate performance metrics using a single rsplit object
#'
#' This function trains a glmnet model on the training set of an rsplit object, then
#' calculates performance metrics of that model on the validation/holdout set
#' at all combinations of the mixture and
#' penalty hyperparameters provided in a hyperparameter grid.
#'
#' @param rsplit_data An rsplit object produced by the \code{\link[rsample]{rsample}} package
#' @param prepped_recipe A trained \code{\link[recipes]{recipe}}
#' @param hyperparameter_grid A tibble containing the hyperparameter values to tune.
#' Can be created using \code{\link{tof_create_grid}}
#' @param model_type A string representing the type of glmnet model being fit.
#' @param outcome_colnames Quoted column names indicating which columns in the data
#' being fit represent the outcome variables (with all others assumed to be predictors).
#'
#' @return A tibble with the same number of rows as the input hyperparameter grid.
#' Each row represents a combination of mixture and penalty, and each column contains
#' a performance metric for the fitted glmnet model on `rsplit_data`'s holdout set.
#' The specific performance metrics depend on the type of model being fit:
#' \describe{
#' \item{"linear"}{mean-squared error (`mse`) and mean absolute error (`mae`)}
#' \item{"two-class"}{binomial deviance (`binomial_deviance`); misclassification error rate
#' `misclassification_error`; the area under the receiver-operating curve (`roc_auc`);
#' and `mse` and `mse` as above}
#' \item{"multiclass"}{multinomial deviance (`multinomial_deviance`); misclassification error rate
#' `misclassification_error`; the area under the receiver-operating curve (`roc_auc`)
#' computed using the Hand-Till method in \code{\link[yardstick]{roc_auc}};
#' and `mse` and `mse` as above}
#' \item{"survival"}{the negative log2-transformed partial likelihood (`neg_log_partial_likelihood`)
#' and Harrel's concordance index (often simply called "C"; `concordance_index`)}
#' }
#'
#' @references Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) Tutorial in biostatistics: multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing error, Statistics in Medicine, 15, pages 361–387.
#'
#'
tof_fit_split <-
  function(
    rsplit_data,
    prepped_recipe,
    hyperparameter_grid,
    model_type,
    outcome_colnames
  ) {

    # extract training set
    training_data <-
      recipes::bake(
        object = prepped_recipe,
        new_data = rsample::training(rsplit_data)
      )

    # # extract outcome column names
    # outcome_colnames <-
    #   rlang::enquo(outcome_cols) %>%
    #   tidyselect::eval_select(data = training_data) %>%
    #   names()

    # separate training and testing x and y (as needed by glmnet)
    training_x <-
      training_data %>%
      dplyr::select(-dplyr::any_of(outcome_colnames)) %>%
      as.matrix()

    training_y <-
      training_data %>%
      dplyr::select(dplyr::any_of(outcome_colnames))

    # extract test set
    validation_data <-
      recipes::bake(
        object = prepped_recipe,
        new_data = rsample::testing(rsplit_data)
      )

    validation_x <-
      validation_data %>%
      dplyr::select(-dplyr::any_of(outcome_colnames)) %>%
      as.matrix()

    validation_y <-
      validation_data %>%
      dplyr::select(dplyr::any_of(outcome_colnames))

    # convert outcome into the form glmnet needs it
    if (model_type == "survival") {
      training_y <-
        survival::Surv(
          time = training_y[[1]],
          event = training_y[[2]]
        )

      validation_y <-
        survival::Surv(
          time = validation_y[[1]],
          event = validation_y[[2]]
        )

    } else {
      training_y <-
        training_y[[outcome_colnames]]

      validation_y <-
        validation_y[[outcome_colnames]]
    }

    if (!missing(hyperparameter_grid)) {
      alphas <-
        hyperparameter_grid %>%
        dplyr::pull(mixture) %>%
        unique()

      lambdas <-
        hyperparameter_grid %>%
        dplyr::pull(penalty) %>%
        unique()

    } else {
      alphas <-
        tof_create_grid() %>%
        dplyr::pull(mixture) %>%
        unique()

      lambdas <- NULL
    }

    # "linear", "two-class", "multiclass", "survival"
    glmnet_family <-
      switch(
        model_type,
        "linear" = "gaussian",
        "two-class" = "binomial",
        "multiclass" = "multinomial",
        "survival" = "cox"
      )

    glmnet_models <-
      purrr::map(
        .x = alphas,
        .f = ~ glmnet::glmnet(
          x = training_x,
          y = training_y,
          alpha = .x,
          lambda = lambdas, # decide if this should be specified or automatic via glmnet
          family = glmnet_family,
          standardize = FALSE
        )
      )

    if (is.null(lambdas)) {
      lambdas <- NULL
    }

    # compute error metrics for each model using the holdout set
    model_metrics <-
      purrr::map(
        .x = glmnet_models,
        .f = ~ glmnet::assess.glmnet(
          object = .x,
          newx = validation_x,
          newy = validation_y,
          s = lambdas
        )  %>%
          tibble::as_tibble() %>%
          dplyr::mutate(penalty = lambdas)
      )

    # tidy model metrics for each alpha and lambda value used to fit a model
    result <-
      tibble::tibble(
        mixture = alphas,
        model_metrics = model_metrics
      ) %>%
      tidyr::unnest(cols = model_metrics) %>%
      dplyr::relocate(mixture, penalty, dplyr::everything())


    # use yardstick to compute the hand-till roc_auc for multiclass problems
    if (model_type == "multiclass") {
      prediction_matrices <-
        purrr::map(
          .x = glmnet_models,
          .f = ~ predict(.x, newx = validation_x, type = "response", s = lambdas)
        )

      roc_auc_tibble <-
        hyperparameter_grid %>%
        dplyr::arrange(mixture, penalty) %>%
        dplyr::mutate(
          predictions =
            purrr::map2(
              .x = sort(rep(1:length(alphas), times = length(lambdas))),
              .y = rep(1:length(lambdas), times = length(alphas)),
              .f = ~ tibble::as_tibble(prediction_matrices[[.x]][,,.y])
            ),
          roc_auc =
            purrr::map_dbl(
              .x = predictions,
              .f = ~
                yardstick::roc_auc(
                  data = .x,
                  truth = validation_y,
                  dplyr::all_of(levels(validation_y))
                ) %>%
                dplyr::pull(.estimate)
            )
        ) %>%
        dplyr::select(penalty, mixture, roc_auc)

      # combine multiclass roc_auc with other performance metrics
      result <-
        result %>%
        dplyr::left_join(roc_auc_tibble, by = c("mixture", "penalty"))

    }

    # clean metric names to make them more descriptive for all model types
    result <-
      result %>%
      tof_clean_metric_names(model_type = model_type)

    return(result)

  }

#' Title
#'
#' TO DO
#'
#' @param metric_tibble TO DO
#' @param model_type TO DO
#'
#' @return TO DO
#'
#'
tof_clean_metric_names <- function(metric_tibble, model_type) {

  if (model_type == "linear") {
    # mse = mse (mean squared error)
    # mae = mae (mean absolute error)
    new_metric_tibble <-
      metric_tibble

  } else if (model_type == "two-class") {
    # binomial_deviance = deviance
    # misclassification_error = class
    # roc_auc = auc
    # mse = mse
    # mae = mae
    new_metric_tibble <-
      metric_tibble %>%
      dplyr::rename(
        binomial_deviance = deviance,
        misclassification_error = class,
        roc_auc = auc
      )

  } else if (model_type == "multiclass") {
    # multinomial_deviance = deviance
    # misclassification_error = class
    # mse = mse
    # mae = mae
    new_metric_tibble <-
      metric_tibble %>%
      dplyr::rename(
        multinomial_deviance = deviance,
        misclassification_error = class
      )

  } else if (model_type == "survival") {
    # partial_likelihood deviance = deviance
    # concordance_index = C
    new_metric_tibble <-
      metric_tibble %>%
      dplyr::rename(
        neg_log_partial_likelihood = deviance,
        concordance_index = C
      )
  }

  return(new_metric_tibble)

}





tof_find_best <- function(performance_metrics, model_type, optimization_metric) {

  # of optimization metric is default, give it tidytof's automatic preference
  # based on the model type
  if (optimization_metric == "tidytof_default") {
    optimization_metric <-
      switch(
        model_type,
        "linear" = "mse",
        "two-class" = "roc_auc",
        "multiclass" = "roc_auc",
        "survival" = "concordance_index"
      )
  }

  # check optimization metric
  if (model_type == "linear") {
    optimization_metric <-
      rlang::arg_match0(optimization_metric, values = c("mse", "mae"))

  } else if (model_type == "two-class") {
    # binomial_deviance = deviance
    # misclassification_error = class
    # roc_auc = auc
    # mse = mse
    # mae = mae
    optimization_metric <-
      rlang::arg_match0(
        optimization_metric,
        values =
          c("roc_auc", "misclassification_error", "binomial_deviance", "mse", "mae")
      )

  } else if (model_type == "multiclass") {
    # multinomial_deviance = deviance
    # misclassification_error = class
    # mse = mse
    # mae = mae
    optimization_metric <-
      rlang::arg_match0(
        optimization_metric,
        values =
          c("roc_auc", "misclassification_error", "multinomial_deviance", "mse", "mae")
      )

  } else if (model_type == "survival") {
    # partial_likelihood deviance = deviance
    # concordance_index = C
    optimization_metric <-
      rlang::arg_match0(
        optimization_metric,
        values =
          c("partial_likelihood_deviance", "concordance_index")
      )
  }

  # find best model parameters
  performance_metrics$optimizer <- performance_metrics[[optimization_metric]]

  if (optimization_metric %in% c("roc_auc", "concordance_index")) {
    best_metrics <-
      performance_metrics %>%
      dplyr::slice_max(order_by = optimizer, n = 1)
  } else {
      best_metrics <-
        performance_metrics %>%
        dplyr::slice_min(order_by = optimizer, n = 1)
  }

  best_hyperparameters <-
    best_metrics %>%
    dplyr::select(penalty, mixture)

  return(best_hyperparameters)

}

tof_finalize_model <-
  function(
    best_model_parameters,
    split_data
  ) {

  }

tof_fit_glmnet <- function(training_data) {
  NULL
}

# new_tof_model ----------------------------------------------------------------

#' Constructor for a tof_model.
#'
#' @param model A glmnet model.
#'
#' @param recipe A prepped recipe object.
#'
#' @return A `tof_model`, an S3 class that includes a trained glmnet model and
#' the recipe used to perform its associated preprocessing.
#'
new_tof_model <- function(model, recipe) {

  stopifnot(inherits(recipe, "recipe"))
  stopifnot(inherits(recipe, "glmnet"))

  tof_model <- list(model = model, recipe = recipe)

  return(tof_model)
}

