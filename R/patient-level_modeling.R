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
#' @importFrom dplyr pull
#' @importFrom rlang arg_match
#' @importFrom rsample initial_split
#' @importFrom rsample vfold_cv
#' @importFrom rsample bootstraps
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

#' Create an elastic net hyperparameter search grid of a specified size
#'
#' This function creates a hyperparameter search grid (in the form of a
#' \code{\link[tibble]{tibble}}) specifying the search space for the two
#' hyperparameters of a generalized linear model using the {glmnet} package:
#' the regularization penalty term
#' and the lasso/ridge regression mixture term. Either a regular or maximum-entropy
#' search grid (of a specified size) or a grid with a user-specified range of
#' hyperparameter values can be created.
#'
#' @param penalty_values A numeric vector of the unique elastic net penalty values ("lambda")
#' to include in the
#' hyperparameter grid. If unspecified, a regular grid with `num_penalty_values` between
#' 10^(-10) and 10^(0) will be used.
#'
#' @param mixture_values  A numeric vector of all elastic net mixture values ("alpha") to include in the
#' hyperparameter grid. If unspecified, a regular grid with `num_mixture_values` between
#' 0 and 1 will be used.
#'
#' @param grid_type A string indicating if a regular grid ("regular"; the default) or maximum-entropy
#' ("entropy") search grid should be returned. Ignored if both `penalty_values` and
#' `mixture_values` are provided. If "entropy" is chosen, the {dials} package
#' is required (and a message will prompt you to install it).
#'
#' @param num_penalty_values Optional. If `penalty_values` is not supplied, `num_penalty_values`
#' (an integer) can be given to specify how many equally-spaced penalty values between
#' 10^(-10) and 1 should be included in the hyperparameter grid. If this method is used,
#' the regular grid will always be returned. Defaults to 5.
#'
#' @param num_mixture_values Optional. If `mixture_values` is not supplied, `num_mixture_values`
#' (an integer) can be given to specify how many equally-spaced penalty values between
#' 0 (ridge regression) and 1 (lasso) should be included in the hyperparameter grid. If this method is used,
#' the regular grid will always be returned. Defaults to 5.
#'
#' @param entropy_grid_size An integer indicating how many combinations of
#' penalty and mixture should be generated by the max entropy hyperparameter grid. Ignored
#' if `grid_type` is not set to "entropy".
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

#' Train an elastic net model to predict sample-level phenomena using CyTOF data.
#'
#' This function uses a training set/test set paradigm to tune and fit an
#' elastic net model using a variety of user-specified details. Tuning can be
#' performed using either a simple training vs. test set split, k-fold cross-validation,
#' or bootstrapping, and multiple preprocessing options are available.
#'
#' @param split_data An `rsplit` or `rset` object from the \code{\link[rsample]{rsample}}
#' package. The easiest way to generate this is to use \code{\link{tof_split_data}}
#'
#' @param predictor_cols Unquoted column names indicating which columns in the
#' data contained in `split_data` should be used as predictors in the elastic net model.
#' Supports tidyselect helpers.
#'
#' @param response_col Unquoted column name indicating which column in the data
#' contained in `split_data` should be used as the outcome in a "two-class", "multiclass",
#' or "linear" elastic net model. Must be a factor for "two-class" and "multiclass"
#' models and must be a numeric for "linear" models. Ignored if `model_type` is "survival".
#'
#' @param time_col Unquoted column name indicating which column in the data
#' contained in `split_data` represents the time-to-event outcome in a "survival"
#' elastic net model. Must be numeric. Ignored if `model_type` is "two-class", "multiclass",
#' or "linear".
#'
#' @param event_col Unquoted column name indicating which column in the data
#' contained in `split_data` represents the time-to-event outcome in a "survival"
#' elastic net model. Must be a binary column - all values should be either 0 or 1
#' (with 1 indicating the adverse event) or FALSE and TRUE (with TRUE indicating the
#' adverse event). Ignored if `model_type` is "two-class", "multiclass",
#' or "linear".
#'
#' @param model_type A string indicating which kind of elastic net model to build.
#' If a continuous response is being predicted, use "linear" for linear regression;
#' if a categorical response with only 2 classes is being predicted, use
#' "two-class" for logistic regression; if a categorical response with more than 2
#' levels is being predicted, use "multiclass" for multinomial regression; and if
#' a time-to-event outcome is being predicted, use "survival" for Cox regression.
#'
#' @param hyperparameter_grid A hyperparameter grid indicating which values of
#' the elastic net penalty (lambda) and the elastic net mixture (alpha) hyperparamters
#' should be used during model tuning. Generate this grid using \code{\link{tof_create_grid}}.
#'
#' @param standardize_predictors A logical value indicating if numeric predictor columns
#' should be standardized (centered and scaled) before model fitting. Defaults to TRUE.
#'
#' @param remove_zv_predictors A logical value indicating if predictor columns
#' with near-zero variance should be removed before model fitting using
#' \code{\link[recipes]{step_nzv}}. Defaults to FALSE.
#'
#' @param impute_missing_predictors A logical value indicating if predictor columns
#' should have missing values imputed using k-nearest neighbors before model fitting (see
#' \code{\link[recipes]{step_impute_knn}}). Imputation is performed using an observation's
#' 5 nearest-neighbors. Defaults to FALSE.
#'
#' @param optimization_metric A string indicating which optimization metric
#' should be used for hyperparameter selection during model tuning. Valid values
#' depend on the model_type.
#'
#' \itemize{
#' \item For "linear" models, choices are "mse" (the mean squared error
#' of the predictions; the default) and "mae" (the mean absolute error of the predcitions).
#' \item For "two-class" models, choices are "roc_auc" (the area under the Receiver-Operating
#' Curve for the classification; the default), "misclassification error" (the proportion of
#' misclassified observations), "binomial_deviance" (see \code{\link[glmnet]{deviance.glmnet}}),
#' "mse" (the mean squared error of the logit function), and "mae" (the mean absolute error of the
#' logit function).
#' \item For "multiclass" models, choices are "roc_auc" (the area under the Receiver-Operating
#' Curve for the classification using the Hand-Till generalization of the ROC AUC
#' for multiclass models in \code{\link[yardstick]{roc_auc}}; the default), "misclassification error"
#' (the proportion of
#' misclassified observations), "multinomial_deviance" (see \code{\link[glmnet]{deviance.glmnet}}),
#' and "mse" and "mae" as above.
#' \item For "survival" models, choices are "concordance_index" (Harrel's C index;
#' see \code{\link[glmnet]{deviance.glmnet}}) and "partial_likelihood_deviance"
#' (see \code{\link[glmnet]{deviance.glmnet}}).
#' }
#'
#' @param best_model_type Currently unused.
#'
#' @param num_cores Integer indicating how many cores should be used for parallel
#' processing when fitting multiple models. Defaults to 1. Overhead to separate
#' models across multiple cores can be high, so significant speedup is unlikely
#' to be observed unless many large models are being fit.
#'
#' @return A `tof_model`, an S3 class that includes the elastic net model with
#' the best performance (assessed via cross-validation, bootstrapping, or simple splitting
#' depending on `split_data`) across all tested hyperparameter value combinations. `tof_models`
#' store the following information:
#'
#' \describe{
#' \item{model}{The final elastic net ("glmnet") model}
#' \item{recipe}{The \code{\link[recipes]{recipe}} used for data preprocessing}
#' \item{mixture}{The optimal mixture hyperparameter (alpha) for the glmnet model}
#' \item{penalty}{The optimal penalty hyperparameter (lambda) for the glmnet model}
#' \item{model_type}{A string indicating which type of glmnet model was fit}
#' \item{outcome_colnames}{A character vector representing the names of the columns in the training data modeled as outcome variables}
#' \item{training_data}{A tibble containing the (not preprocessed) data used to train the model}
#' }
#'
#' @export
#'
#' @importFrom rlang as_name
#' @importFrom rlang ensym
#'
tof_train_model <-
  function(
    split_data,
    predictor_cols,
    response_col = NULL,
    time_col = NULL,
    event_col = NULL,
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

    # extract response column names
    response_colname <- time_colname <- event_colname <- NULL

    if (model_type %in% c("linear", "two-class", "multiclass")) {
      response_colname <- rlang::as_name(rlang::ensym(response_col))
    } else {
      time_colname <- rlang::as_name(rlang::ensym(time_col))
      event_colname <- rlang::as_name(rlang::ensym(event_col))
    }

    tof_model <-
      all_data %>%
      tof_finalize_model(
        best_model_parameters = best_model_parameters,
        recipe = final_recipe,
        model_type = model_type,
        outcome_colnames = c(response_colname, time_colname, event_colname)
      )

    return(tof_model)

  }


# applying models to new data --------------------------------------------------

#' Use a trained elastic net model to predict fitted values from new data
#'
#' This function uses a trained `tof_model` to make predictions on new data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations for which predictions should be made.
#' If new_data isn't provided, predictions will be made for the training data used to
#' fit the model.
#'
#' @param prediction_type A string indicating which type of prediction should be
#' provided by the model:
#'
#' \describe{
#' \item{"response" (the default)}{For "linear" models, the predicted response
#' for each observation. For "two-class" and "multiclass" models, the fitted
#' probabilities of each class for each observation. For "survival" models, the fitted relative-risk
#' for each observation.}
#' \item{"class}{Only applies to "two-class" and "multiclass" models. For both,
#' the class label corresponding to the class with the maximum fitted probability.}
#' \item{"link"}{The linear predictions of the model
#' (the output of the link function for each model family.)}
#' }
#'
#' @return A \code{link[tibble]{tibble}} with a single column (`.pred`) containing
#' the predictions. Each row in the output corresponds to a row in `new_data` (
#' or, if `new_data` is not provided, to a row in the `tof_model`'s training data.
#'
#' @export
#'
#' @importFrom dplyr any_of
#' @importFrom tibble tibble
#' @importFrom rlang arg_match
#'
tof_predict <-
  function(
    tof_model,
    new_data,
    prediction_type = c("response", "class", "link")
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

    # if no new data is provided, simply return the training data predictions
    if (missing(new_data)) {
      new_data <- tof_get_model_training_data(model)
    }

    # preprocess the new_data
    preprocessed_data <-
      new_data %>%
      tof_setup_glmnet_xy(
        recipe = recipe,
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

    return(tibble::tibble(.pred = predictions))

  }

# assessing model performance --------------------------------------------------

#' Assess a trained elastic net model when applied to new data
#'
#' This function assesses a trained `tof_model`'s performance on new data by
#' computing model type-specific performance measurements.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations that should be used to evaluate `tof_model`.
#' If new_data isn't provided, model evaluation will will be performed using the
#' training data used to fit the model.
#'
#' @return A list of performance metrics that depends on the model type.
#'
#' @export
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr everything
#' @importFrom glmnet assess.glmnet
#' @importFrom glmnet confusion.glmnet
#' @importFrom glmnet roc.glmnet
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom survival survfit
#' @importFrom purrr discard
#'
tof_assess_model <-
  function(
    tof_model,
    new_data
  ) {

    # set up -------------------------------------------------------------------
    # extract the recipe and the glmnet model
    model <- tof_model$model
    recipe <- tof_model$recipe
    lambda <- tof_model$penalty

    # extract the model type and outcome colnames
    model_type <- tof_get_model_type(tof_model)
    outcome_colnames <- tof_get_model_outcomes(tof_model)

    # if no new data is supplied, simply return the training data assessment
    if (missing(new_data)) {
      new_data <- tof_get_model_training_data(tof_model)
    }

    # preprocess the new data
    new_xy <-
      new_data %>%
      tof_setup_glmnet_xy(
        recipe = recipe,
        outcome_cols = dplyr::any_of(outcome_colnames),
        model_type = model_type
      )

    # calculate metrics --------------------------------------------------------

    model_metrics <-
      glmnet::assess.glmnet(
        object = model,
        newx = new_xy$x,
        newy = new_xy$y,
        family = tof_find_glmnet_family(model_type),
        s = lambda
      ) %>%
      tibble::as_tibble() %>%
      tof_clean_metric_names(model_type = model_type) %>%
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = "metric",
        values_to = "value"
      )

    if (model_type == "two-class") {
      auc <- model_metrics$roc_auc

      roc_curve <-
        glmnet::roc.glmnet(
          object = model,
          newx = new_xy$x,
          newy = new_xy$y,
          family = "binomial",
          s = lambda
        ) %>%
        tibble::as_tibble()

    } else {
      roc_curve <- NULL
    }

    if (model_type %in% c("two-class", "multiclass")) {
      confusion_matrix <-
        glmnet::confusion.glmnet(
          object = model,
          newx = new_xy$x,
          newy = new_xy$y,
          family = tof_find_glmnet_family(model_type),
          s = lambda
        ) %>%
        tibble::as_tibble() %>%
        dplyr::transmute(
          true_outcome = True,
          predicted_outcome = Predicted,
          num_observations = n
        )

    } else {
      confusion_matrix <- NULL
    }

    if (model_type == "survival") {
      survfit_result <-
        survival::survfit(
          model,
          s = lambda,
          x = tof_get_model_x(tof_model),
          y = tof_get_model_y(tof_model),
          newx = new_xy$x
        )

      times <-
        tibble::tibble(
          time = survfit_result$time,
          .timepoint_index = 1:length(survfit_result$time)
        )

      survival_curves <-
        survfit_result$surv %>%
        as_tibble() %>%
        dplyr::mutate(.timepoint_index = 1:nrow(survfit_result$surv)) %>%
        tidyr::pivot_longer(
          cols = -.timepoint_index,
          names_to = "row_index",
          values_to = "probability"
        ) %>%
        dplyr::left_join(times, by = ".timepoint_index") %>%
        dplyr::select(-.timepoint_index) %>%
        tidyr::nest(survival_curve = -row_index)

      return(survival_curves)

    } else {
      survival_curve <- NULL
    }

    # return result ------------------------------------------------------------

    result <-
      list(
        model_metrics = model_metrics,
        roc_curve = roc_curve,
        confusion_matrix = confusion_matrix,
        survival_curve = survival_curve
      )

    result <- purrr::discard(result, is.null)

    return(result)
  }
