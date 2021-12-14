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
#' @param num_penalty_values TO DO
#'
#' @param num_mixture_values TO DO
#'
#' @param entropy_grid_size TO DO
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
    num_penalty_values = 5,
    num_mixture_values = 5,
    grid_type = c("regular", "entropy"),
    entropy_grid_size = 25
  ) {

    # check grid_type argument
    grid_type <- rlang::arg_match(grid_type)

    if (grid_type == "entropy") {
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
        dials::grid_max_entropy(hyperparams, size = entropy_grid_size)

    } else if (grid_type == "regular") {

      if (missing(penalty_values)) {
        penalty_values <- 10^(seq(-10, 0, length.out = num_penalty_values))
      }

      if (missing(mixture_values)) {
        mixture_values <- seq(0, 1, length.out = num_mixture_values)
      }

      hyperparam_grid <-
        tidyr::expand_grid(penalty_values, mixture_values) %>%
        dplyr::rename(
          penalty = penalty_values,
          mixture = mixture_values
        )

    } else {
      stop("The grid type must be either regular or entropy.")
    }

    return(hyperparam_grid)
  }



# Building Models --------------------------------------------------------------

#' TO DO
#'
#' TO DO
#'
#' @param split_data TO DO
#' @param predictor_cols TO DO
#' @param response_col TO DO
#' @param time_col TO DO
#' @param event_col TO DO
#' @param model_type TO DO
#' @param hyperparameter_grid TO DO
#' @param standardize_predictors TO DO
#' @param remove_zv_predictors TO DO
#' @param impute_missing_predictors TO DO
#' @param optimization_metric TO DO
#' @param best_model_type TO DO
#' @param num_cores TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_train_model <-
  function(
    split_data,
    predictor_cols,
    response_col = NULL, # classification and regression
    time_col = NULL, # survival
    event_col = NULL, # survival
    model_type = c("linear", "two-class", "multiclass", "survival"),
    hyperparameter_grid = tof_create_grid(),
    standardize_predictors = TRUE,
    remove_zv_predictors = FALSE,
    impute_missing_predictors = FALSE,
    optimization_metric = "tidytof_default",
    best_model_type = c("best", "best with sparsity"), # not currently used
    num_cores = 1
  ) {

    # check arguments ----------------------------------------------------------
    feature_tibble <-
      tof_check_model_args(
        split_data = split_data,
        model_type = model_type,
        best_model_type = best_model_type,
        response_col = {{response_col}},
        time_col = {{time_col}},
        event_col = {{event_col}}
      )

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
    tuning_metrics <-
      tof_tune_glmnet(
        split_data = split_data,
        prepped_recipe = prepped_recipe,
        hyperparameter_grid = hyperparameter_grid,
        model_type = model_type,
        outcome_cols = c({{response_col}}, {{time_col}}, {{event_col}}),
        optimization_metric = optimization_metric,
        num_cores = num_cores
      )

    # find the best hyperparameters using the tuning metrics -------------------


    # find hyperparameters with best performance
    best_model_parameters <-
      tuning_metrics %>%
      tof_find_best(model_type = model_type, optimization_metric = optimization_metric)

    #return(list(tuning_metrics = tuning_metrics, best_model_parameters = best_model_parameters))

    # combine data into a single preprocessing pipeline and glmnet model--------
    all_data <- tof_all_data(split_data = split_data)

    final_recipe <-
      all_data %>%
      tof_prep_recipe(unprepped_recipe = unprepped_recipe)

    all_data_glmnet <-
      recipes::bake(
        object = final_recipe,
        new_data = all_data
      ) %>%
      tof_setup_glmnet_xy(
        outcome_cols = c({{response_col}}, {{time_col}}, {{event_col}}),
        model_type = model_type
      )

    # extract response column names
    response_colname <- time_colname <- event_colname <- NULL

    if (model_type %in% c("linear", "two-class", "multiclass")) {
      response_colname <- rlang::as_name(rlang::ensym(response_col))
    } else {
      time_colname <- rlang::as_name(rlang::ensym(time_col))
      event_colname <- rlang::as_name(rlang::ensym(event_col))
    }

    tof_model <-
      tof_finalize_model(
        best_model_parameters = best_model_parameters,
        recipe = final_recipe,
        glmnet_x = all_data_glmnet$x,
        glmnet_y = all_data_glmnet$y,
        model_type = model_type,
        outcome_colnames = c(response_colname, time_colname, event_colname)
      )

    return(tof_model)

  }


# applying models to new data --------------------------------------------------

#' TO DO
#'
#' TO DO
#'
#' @param tof_model TO DO
#' @param new_data TO DO
#' @param prediction_type TO DO
#' @param ... TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_predict <-
  function(
    tof_model,
    new_data,
    prediction_type = c("response", "class", "link"),
    ...
  ) {

    # set up -------------------------------------------------------------------
    # check prediction_type
    prediction_type <- rlang::arg_match(prediction_type)

    # extract the recipe and the glmnet model
    model <- tof_model$model
    recipe <- tof_model$recipe
    lambda <- tof_get_model_penalty(tof_model)

    # extract the model type and outcome colnames
    model_type <- tof_get_model_type(tof_model)
    outcome_colnames <- tof_get_model_outcomes(tof_model)

    # preprocess the new_data
    preprocessed_data <-
      recipes::bake(
        recipe,
        new_data = new_data
      ) %>%
      tof_setup_glmnet_xy(
        outcome_cols = dplyr::any_of(outcome_colnames),
        model_type = model_type
      )

    # make predictions ---------------------------------------------------------
    predictions <-
      predict(
        object = model,
        newx = preprocessed_data$x,
        s = lambda,
        type = prediction_type,
        exact = TRUE,
        x = tof_get_model_x(tof_model),
        y = tof_get_model_y(tof_model)
      ) %>%
      as.vector()

    return(predictions)

  }

# assessing model performance --------------------------------------------------
#' TO DO
#'
#' TO DO
#'
#' @param tof_model TO DO
#' @param new_data TO DO
#' @param ... TO DO
#'
#' @return TO DO
#'
#' @export
#'
tof_assess_model <-
  function(
    tof_model,
    new_data = NULL,
    ...
  ) {

    # set up -------------------------------------------------------------------
    # extract the recipe and the glmnet model
    model <- tof_model$model
    recipe <- tof_model$recipe
    lambda <- tof_get_model_penalty(tof_model)

    # extract the model type and outcome colnames
    model_type <- tof_get_model_type(tof_model)
    outcome_colnames <- tof_get_model_outcomes(tof_model)

    # preprocess the new_data
    if (!missing(new_data)) {
      preprocessed_data <-
        recipes::bake(
          recipe,
          new_data = new_data
        ) %>%
        tof_setup_glmnet_xy(
          outcome_cols = dplyr::any_of(outcome_colnames),
          model_type = model_type
        )
    } else {
      preprocessed_data <-
        list(
          x = tof_get_model_x(tof_model),
          y = tof_get_model_y(tof_model)
        )
    }

    # calculate metrics --------------------------------------------------------

    model_metrics <-
      glmnet::assess.glmnet(
        object = model,
        newx = preprocessed_data$x,
        newy = preprocessed_data$y,
        family = tof_find_glmnet_family(model_type),
        s = lambda
      ) %>%
      tibble::as_tibble() %>%
      tof_clean_metric_names(model_type = model_type)

    return(model_metrics)


  }
