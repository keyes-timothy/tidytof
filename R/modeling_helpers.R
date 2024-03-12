# modeling_helpers.R
#
# This file contains helper functions for performing sample- and patient-level
# modeling on tibble containing aggregated single-cell data.


#' Create a recipe for preprocessing sample-level cytometry data for an elastic net model
#'
#' @param feature_tibble A tibble in which each row represents a sample- or patient-
#' level observation, such as those produced by \code{tof_extract_features}.
#'
#' @param predictor_cols Unquoted column names indicating which columns in the
#' data contained in `feature_tibble` should be used as predictors in the elastic net model.
#' Supports tidyselect helpers.
#'
#' @param outcome_cols Unquoted column names indicating which columns in
#' `feature_tibble` should be used as outcome variables in the elastic net model.
#' Supports tidyselect helpers.
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
#' @return A \code{\link[recipes]{recipe}} object.
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr select
#'
#'
tof_create_recipe <-
    function(
        feature_tibble,
        predictor_cols,
        outcome_cols,
        standardize_predictors = TRUE,
        remove_zv_predictors = FALSE,
        impute_missing_predictors = FALSE) {
        # Check that recipes is installed
        has_recipes <- requireNamespace("recipes", quietly = TRUE)
        if (!has_recipes) {
            stop(
                "This function requires the {recipes} package. Install it with this code:\n",
                "install.packages(\"recipes\")\n"
            )
        }

        # extract predictor column names
        predictor_colnames <-
            feature_tibble |>
            dplyr::select({{ predictor_cols }}) |>
            colnames()

        outcome_colnames <-
            feature_tibble |>
            dplyr::select({{ outcome_cols }}) |>
            colnames()

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
            ) |>
            recipes::step_select(
                dplyr::any_of(outcome_colnames),
                dplyr::any_of(predictor_colnames)
            ) |>
            recipes::step_dummy(recipes::all_nominal_predictors())

        if (remove_zv_predictors) {
            model_recipe <-
                model_recipe |>
                recipes::step_nzv(recipes::all_predictors())
        }

        if (impute_missing_predictors) {
            model_recipe <-
                model_recipe |>
                recipes::step_impute_knn(recipes::all_numeric_predictors())
        }

        if (standardize_predictors) {
            model_recipe <-
                model_recipe |>
                recipes::step_normalize(recipes::all_numeric_predictors(), na_rm = TRUE)
        }

        return(model_recipe)
    }



#' Train a recipe or list of recipes for preprocessing sample-level cytometry data
#'
#' @param split_data An `rsplit` or `rset` object from the \code{\link[rsample]{rsample}}
#' package containing the sample-level data to use for modeling.
#' The easiest way to generate this is to use \code{\link{tof_split_data}}.
#' Alternatively, an unsplit tbl_df, though this is not recommended.
#'
#' @param unprepped_recipe A \code{\link[recipes]{recipe}} object (if `split_data`
#' is an `rsplit` object or a `tbl_df`) or list of recipes
#' (if `split_data` is an `rset` object).
#'
#' @return If split_data is an "rsplit" or "tbl_df" object, will return a single
#' prepped recipe. If split_data is an "rset" object, will return a list of prepped
#' recipes specific for each fold of the resampling procedure.
#'
#' @importFrom purrr map
#'
#'
tof_prep_recipe <-
    function(
        split_data,
        unprepped_recipe) {
        # Check that recipes is installed
        has_recipes <- requireNamespace("recipes", quietly = TRUE)
        if (!has_recipes) {
            stop(
                "This function requires the {recipes} package. Install it with this code:\n",
                "install.packages(\"recipes\")\n"
            )
        }

        # Check that rsample is installed only if split_data is "rset" or "rsplit"
        if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
            has_rsample <- requireNamespace("rsample", quietly = TRUE)
            if (!has_rsample) {
                stop(
                    "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                    "Install it with:\ninstall.packages(\"rsample\")\n"
                )
            }
        }

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
                unprepped_recipe |>
                recipes::prep(training = rsample::training(split_data), retain = FALSE)
        } else if (inherits(split_data, "tbl_df")) {
            result <-
                unprepped_recipe |>
                recipes::prep(training = split_data, retain = FALSE)
        } else {
            stop("split_data must be a tibble, an rset, or an rsplit object")
        }

        return(result)
    }

#' Tune an elastic net model's hyperparameters over multiple resamples
#'
#' @param split_data An `rsplit` or `rset` object from the \code{\link[rsample]{rsample}}
#' package. The easiest way to generate this is to use \code{\link{tof_split_data}}.
#' Alternatively, an unsplit tbl_df can be provided, though this is not recommended.
#'
#' @param prepped_recipe Either a single \code{\link[recipes]{recipe}} object
#' (if `split_data` is an `rsplit` object or a `tbl_df`) or list of recipes
#' (if `split_data` is an `rset` object) such that each entry in the list
#' corresponds to a resample in `split_data`.
#'
#' @param hyperparameter_grid A hyperparameter grid indicating which values of
#' the elastic net penalty (lambda) and the elastic net mixture (alpha) hyperparameters
#' should be used during model tuning. Generate this grid using \code{\link{tof_create_grid}}.
#'
#' @param model_type A string indicating which kind of elastic net model to build.
#' If a continuous response is being predicted, use "linear" for linear regression;
#' if a categorical response with only 2 classes is being predicted, use
#' "two-class" for logistic regression; if a categorical response with more than 2
#' levels is being predicted, use "multiclass" for multinomial regression; and if
#' a time-to-event outcome is being predicted, use "survival" for Cox regression.
#'
#' @param outcome_cols Unquoted column name(s) indicating which column(s) in the data
#' contained in `split_data` should be used as the outcome in the elastic net model.
#' For survival models, two columns should be selected; for all others, only one
#' column should be selected.
#'
#' @param optimization_metric A string indicating which optimization metric
#' should be used for hyperparameter selection during model tuning. Valid values
#' depend on the model_type.
#'
#' @param num_cores Integer indicating how many cores should be used for parallel
#' processing when fitting multiple models. Defaults to 1. Overhead to separate
#' models across multiple cores can be high, so significant speedup is unlikely
#' to be observed unless many large models are being fit.
#'
#' @return A tibble containing a summary of the model's performance in each
#' resampling iteration across all hyperparameter combinations. Will contain
#' 3 columns: "splits" (a list-col containing each resampling iteration's
#' `rsplit` object), "id" (the name of the resampling iteration), and
#' "performance_metrics" (a list-col containing the performance metrics for each
#' resampling iteration. Each row of "performance_metrics" is a tibble with
#' the columns "mixture" and "penalty" and several additional columns containing the
#' performance metrics of the model for each mixture/penalty combination).
#' See \code{\link{tof_fit_split}} for additional details.
#'
#' @importFrom doParallel registerDoParallel
#'
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @importFrom purrr map2
#'
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
        num_cores = 1) {
        # Check that rsample is installed only if split_data is "rset" or "rsplit"
        if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
            has_rsample <- requireNamespace("rsample", quietly = TRUE)
            if (!has_rsample) {
                stop(
                    "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                    "Install it with:\ninstall.packages(\"rsample\")\n"
                )
            }
        }

        # if the input data is a single training/test split
        if (inherits(split_data, "rsplit")) {
            # extract output column names
            outcome_colnames <-
                rlang::enquo(outcome_cols) |>
                tidyselect::eval_select(data = rsample::training(split_data)) |>
                names()

            # compute performance metrics for single split
            performance_metrics <-
                tof_fit_split(
                    split_data,
                    prepped_recipe,
                    hyperparameter_grid,
                    model_type,
                    outcome_colnames
                )

            # if the input data is multiple training/test splits
        } else if (inherits(split_data, "rset")) {
            # extract output column names
            outcome_colnames <-
                split_data$splits[[1]] |>
                rsample::training() |>
                dplyr::select({{ outcome_cols }}) |>
                colnames()

            # compute performance metrics
            if (num_cores == 1) {
                performance_metrics <-
                    purrr::map2(
                        .x = split_data$splits,
                        .y = prepped_recipe,
                        .f = ~ tof_fit_split(
                            split_data = .x,
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
                        split_data = split_data,
                        prepped_recipe = prepped_recipe,
                        hyperparameter_grid = hyperparameter_grid,
                        model_type = model_type,
                        outcome_colnames = outcome_colnames
                    )

                parallel::stopCluster(my_cluster)
            }

            performance_metrics <-
                split_data |>
                dplyr::mutate(
                    performance_metrics = performance_metrics
                )

            # if input data is an unsplit tbl_df
        } else if (inherits(split_data, "tbl_df")) {
            # extract output column names
            outcome_colnames <-
                split_data |>
                dplyr::select({{ outcome_cols }}) |>
                colnames()

            # compute performance metrics
            performance_metrics <-
                tof_fit_split(
                    split_data,
                    prepped_recipe,
                    hyperparameter_grid,
                    model_type,
                    outcome_colnames
                )
        } else {
            stop("split_data must be an rsplit object, an rset object, or a tbl_df.")
        }

        return(performance_metrics)
    }


#' Fit a glmnet model and calculate performance metrics using a single rsplit object
#'
#' This function trains a glmnet model on the training set of an rsplit object, then
#' calculates performance metrics of that model on the validation/holdout set
#' at all combinations of the mixture and
#' penalty hyperparameters provided in a hyperparameter grid.
#'
#' @param split_data An `rsplit` object from the \code{\link[rsample]{rsample}}
#' package.
#' Alternatively, an unsplit tbl_df can be provided, though this is not recommended.
#'
#' @param prepped_recipe A trained \code{\link[recipes]{recipe}}
#'
#' @param hyperparameter_grid A tibble containing the hyperparameter values to tune.
#' Can be created using \code{\link{tof_create_grid}}
#'
#' @param model_type A string representing the type of glmnet model being fit.
#'
#' @param outcome_colnames Quoted column names indicating which columns in the data
#' being fit represent the outcome variables (with all others assumed to be predictors).
#'
#' @return A tibble with the same number of rows as the input hyperparameter grid.
#' Each row represents a combination of mixture and penalty, and each column contains
#' a performance metric for the fitted glmnet model on `split_data`'s holdout set.
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
#' @importFrom dplyr any_of
#' @importFrom dplyr arrange
#' @importFrom dplyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr relocate
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
#' @importFrom glmnet glmnet
#' @importFrom glmnet assess.glmnet
#'
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr map_dbl
#'
#' @importFrom survival Surv
#'
#' @importFrom tidyr unnest
#'
#' @importFrom yardstick roc_auc
#'
#'
tof_fit_split <-
    function(
        split_data,
        prepped_recipe,
        hyperparameter_grid,
        model_type,
        outcome_colnames) {
        # Check that recipes is installed
        has_recipes <- requireNamespace("recipes", quietly = TRUE)
        if (!has_recipes) {
            stop(
                "This function requires the {recipes} package. Install it with this code:\n",
                "install.packages(\"recipes\")\n"
            )
        }

        # Check that rsample is installed only if split_data is "rset" or "rsplit"
        if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
            has_rsample <- requireNamespace("rsample", quietly = TRUE)
            if (!has_rsample) {
                stop(
                    "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                    "Install it with:\ninstall.packages(\"rsample\")\n"
                )
            }
        }


        # extract training and validation set
        if (inherits(split_data, "rsplit")) {
            training_data <- rsample::training(split_data)
            validation_data <- rsample::testing(split_data)
        } else if (inherits(split_data, "tbl_df")) {
            training_data <- split_data
            validation_data <- split_data
        }

        training_data <-
            recipes::bake(object = prepped_recipe, new_data = training_data)

        validation_data <-
            recipes::bake(object = prepped_recipe, new_data = validation_data)

        # separate training and testing x and y (as needed by glmnet)
        training_x <-
            training_data |>
            dplyr::select(-dplyr::any_of(outcome_colnames)) |>
            as.matrix()

        training_y <-
            training_data |>
            dplyr::select(dplyr::any_of(outcome_colnames))

        validation_x <-
            validation_data |>
            dplyr::select(-dplyr::any_of(outcome_colnames)) |>
            as.matrix()

        validation_y <-
            validation_data |>
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
                hyperparameter_grid |>
                dplyr::pull(.data$mixture) |>
                unique()

            lambdas <-
                hyperparameter_grid |>
                dplyr::pull(.data$penalty) |>
                unique()
        } else {
            alphas <-
                tof_create_grid() |>
                dplyr::pull(.data$mixture) |>
                unique()

            lambdas <- NULL
        }

        # convert model type into a glmnet family
        glmnet_family <-
            switch(model_type,
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
                    lambda = lambdas,
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
                ) |>
                    dplyr::as_tibble() |>
                    dplyr::mutate(penalty = lambdas)
            )

        # tidy model metrics for each alpha and lambda value used to fit a model
        result <-
            dplyr::tibble(
                mixture = alphas,
                model_metrics = model_metrics
            ) |>
            tidyr::unnest(cols = model_metrics) |>
            dplyr::relocate("mixture", "penalty", dplyr::everything())


        # use yardstick to compute the hand-till roc_auc for multiclass problems
        if (model_type == "multiclass") {
            prediction_matrices <-
                purrr::map(
                    .x = glmnet_models,
                    .f = ~ predict(.x, newx = validation_x, type = "response", s = lambdas)
                )

            roc_auc_tibble <-
                hyperparameter_grid |>
                dplyr::arrange(.data$mixture, .data$penalty) |>
                dplyr::mutate(
                    predictions =
                        purrr::map2(
                            .x =
                                sort(
                                    rep(seq_len(length(alphas)), times = length(lambdas))
                                ),
                            .y = rep(seq_len(length(lambdas)), times = length(alphas)),
                            .f = ~ dplyr::as_tibble(prediction_matrices[[.x]][, , .y])
                        ),
                    roc_auc =
                        purrr::map_dbl(
                            .x = .data$predictions,
                            .f = ~
                                yardstick::roc_auc_vec(
                                    estimate = as.matrix(.x[levels(validation_y)]),
                                    truth = validation_y
                                )
                        )
                ) |>
                dplyr::select("penalty", "mixture", "roc_auc")

            # combine multiclass roc_auc with other performance metrics
            result <-
                result |>
                dplyr::left_join(roc_auc_tibble, by = c("mixture", "penalty"))
        }

        # clean metric names to make them more descriptive for all model types
        result <-
            result |>
            tof_clean_metric_names(model_type = model_type)

        attr(result, which = "models") <-
            dplyr::tibble(
                mixture = alphas,
                model = glmnet_models
            )

        return(result)
    }


#' Calculate and store the predicted outcomes for each validation set observation during model tuning
#'
#' @param split_data An `rsplit` object from the \code{\link[rsample]{rsample}}
#' package.
#' Alternatively, an unsplit tbl_df can be provided, though this is not recommended.
#'
#' @param prepped_recipe A trained \code{\link[recipes]{recipe}}
#'
#' @param lambda A single numeric value indicating which penalty (lambda) value
#' should be used to make the predictions
#'
#' @param alpha A single numeric value indicating which mixture (alpha) value
#' should be used to make the predictions
#'
#' @param model_type A string indicating which kind of elastic net model to build.
#' If a continuous response is being predicted, use "linear" for linear regression;
#' if a categorical response with only 2 classes is being predicted, use
#' "two-class" for logistic regression; if a categorical response with more than 2
#' levels is being predicted, use "multiclass" for multinomial regression; and if
#' a time-to-event outcome is being predicted, use "survival" for Cox regression.
#'
#' @param outcome_colnames Quoted column names indicating which columns in the data
#' being fit represent the outcome variables (with all others assumed to be predictors).
#'
#' @return A tibble containing the predicted and true values for the outcome
#' for each of the validation observations in `split_data`.
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
#' @importFrom glmnet glmnet
#'
#'
#' @importFrom stats C
#' @importFrom stats deviance
#' @importFrom stats predict
#'
#' @importFrom survival Surv
#' @importFrom survival survfit
#'
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#'
tof_find_cv_predictions <-
    function(
        split_data,
        prepped_recipe,
        lambda,
        alpha,
        model_type,
        outcome_colnames) {
        # Check that recipes is installed
        has_recipes <- requireNamespace("recipes", quietly = TRUE)
        if (!has_recipes) {
            stop(
                "This function requires the {recipes} package. Install it with this code:\n",
                "install.packages(\"recipes\")\n"
            )
        }


        # Check that rsample is installed only if split_data is "rset" or "rsplit"
        if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
            has_rsample <- requireNamespace("rsample", quietly = TRUE)
            if (!has_rsample) {
                stop(
                    "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                    "Install it with:\ninstall.packages(\"rsample\")\n"
                )
            }
        }

        # extract training set
        training_data <-
            recipes::bake(
                object = prepped_recipe,
                new_data = rsample::training(split_data)
            )

        # separate training and testing x and y (as needed by glmnet)
        training_x <-
            training_data |>
            dplyr::select(-dplyr::any_of(outcome_colnames)) |>
            as.matrix()

        training_y <-
            training_data |>
            dplyr::select(dplyr::any_of(outcome_colnames))

        # extract test set
        validation_data <-
            recipes::bake(
                object = prepped_recipe,
                new_data = rsample::testing(split_data)
            )

        validation_x <-
            validation_data |>
            dplyr::select(-dplyr::any_of(outcome_colnames)) |>
            as.matrix()

        validation_y <-
            validation_data |>
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

        # convert model type into a glmnet family
        glmnet_family <-
            switch(model_type,
                "linear" = "gaussian",
                "two-class" = "binomial",
                "multiclass" = "multinomial",
                "survival" = "cox"
            )

        glmnet_model <-
            glmnet::glmnet(
                x = training_x,
                y = training_y,
                alpha = alpha,
                family = glmnet_family,
                standardize = FALSE
            )

        response <-
            stats::predict(
                object = glmnet_model,
                newx = validation_x,
                s = lambda,
                type = "response"
            )

        if (model_type %in% c("linear", "two-class")) {
            predictions <-
                dplyr::tibble(
                    response = as.numeric(response),
                    truth = validation_y
                )

            if (model_type == "two-class") {
                outcome_levels <- levels(predictions$truth)
            }
        } else if (model_type == "multiclass") {
            outcome_levels <- levels(validation_y)

            predictions <-
                (response[, , 1]) |>
                as.matrix() |>
                dplyr::as_tibble() |>
                dplyr::rename_with(.fn = ~ paste0("prob_", .x)) |>
                dplyr::mutate(truth = validation_y)
        } else {
            predictions <-
                dplyr::tibble(
                    response = as.numeric(response),
                    true_time_to_event = as.matrix(validation_y)[, 1],
                    true_event = as.matrix(validation_y)[, 2]
                )
        }

        if (model_type %in% c("multiclass", "two-class")) {
            class <-
                stats::predict(
                    object = glmnet_model,
                    newx = validation_x,
                    s = lambda,
                    type = "class"
                ) |>
                factor(levels = outcome_levels)

            predictions <-
                predictions |>
                dplyr::mutate(class = class)
        }

        if (model_type == "survival") {
            survfit_result <-
                survival::survfit(
                    glmnet_model,
                    s = lambda,
                    x = training_x,
                    y = training_y,
                    newx = validation_x
                )

            times <-
                dplyr::tibble(
                    time = survfit_result$time,
                    .timepoint_index = seq_len(length(survfit_result$time))
                )

            survival_curves <-
                survfit_result$surv |>
                dplyr::as_tibble() |>
                dplyr::mutate(.timepoint_index = seq_len(nrow(survfit_result$surv))) |>
                tidyr::pivot_longer(
                    cols = -".timepoint_index",
                    names_to = "row_index",
                    values_to = "probability"
                ) |>
                dplyr::left_join(times, by = ".timepoint_index") |>
                dplyr::select(-".timepoint_index") |>
                tidyr::nest(survival_curve = -"row_index")

            predictions <-
                predictions |>
                dplyr::mutate(survival_curve = survival_curves$survival_curve) |>
                dplyr::rename(relative_risk = "response")
        }

        return(predictions)
    }

#' Rename glmnet's default model evaluation metrics to make them more interpretable
#'
#' @param metric_tibble A tibble in which each column represents a glmnet
#' model evaluation metric with its default name.
#'
#' @param model_type A string indicating which type of glmnet model was trained.
#'
#' @return A tibble in which each column represents a glmnet
#' model evaluation metric with its "cleaned" name.
#'
#' @importFrom dplyr rename
#'
#'
tof_clean_metric_names <- function(metric_tibble, model_type) {
    if (model_type == "linear") {
        # mse = mse (mean squared error)
        # mae = mae (mean absolute error)
        new_metric_tibble <- metric_tibble
    } else if (model_type == "two-class") {
        # binomial_deviance = deviance
        # misclassification_error = class
        # roc_auc = auc
        # mse = mse
        # mae = mae
        new_metric_tibble <-
            metric_tibble |>
            dplyr::rename(
                binomial_deviance = "deviance",
                misclassification_error = "class",
                roc_auc = "auc"
            ) |>
            dplyr::mutate(accuracy = 1 - .data$misclassification_error)
    } else if (model_type == "multiclass") {
        # multinomial_deviance = deviance
        # misclassification_error = class
        # mse = mse
        # mae = mae
        new_metric_tibble <-
            metric_tibble |>
            dplyr::rename(
                multinomial_deviance = "deviance",
                misclassification_error = "class"
            ) |>
            dplyr::mutate(accuracy = 1 - .data$misclassification_error)
    } else if (model_type == "survival") {
        # partial_likelihood deviance = deviance
        # concordance_index = C
        new_metric_tibble <-
            metric_tibble |>
            dplyr::rename(
                neg_log_partial_likelihood = "deviance",
                concordance_index = "C"
            )
    }

    return(new_metric_tibble)
}





#' Find the optimal hyperparameters for an elastic net model from candidate performance metrics
#'
#' @param performance_metrics A tibble of performance metrics for an elastic
#' net model (in wide format)
#'
#' @param model_type A string indicating which type of glmnet model was trained.
#'
#' @param optimization_metric A string indicating which performance metric should
#' be used to select the optimal model.
#'
#' @return A tibble with 3 columns: "mixture", "penalty", and a column containing
#' the chosen optimization metric. If the returned tibble has more than 1 column,
#' it means that more than 1 mixture/penalty combination yielded the optimal
#' result (i.e. the tuning procedure resulted in a tie).
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr select
#' @importFrom dplyr slice_max
#' @importFrom dplyr slice_min
#'
#' @importFrom rlang arg_match0
#'
tof_find_best <- function(performance_metrics, model_type, optimization_metric) {
    # if optimization metric is default, give it tidytof's automatic preference
    # based on the model type
    if (optimization_metric == "tidytof_default") {
        optimization_metric <-
            switch(model_type,
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
            performance_metrics |>
            dplyr::slice_max(order_by = .data$optimizer, n = 1)
    } else {
        best_metrics <-
            performance_metrics |>
            dplyr::slice_min(order_by = .data$optimizer, n = 1)
    }

    best_hyperparameters <-
        best_metrics |>
        dplyr::select(
            "penalty",
            "mixture",
            dplyr::any_of(optimization_metric)
        )

    return(best_hyperparameters)
}

tof_find_glmnet_family <- function(model_type) {
    glmnet_family <-
        switch(model_type,
            "linear" = "gaussian",
            "two-class" = "binomial",
            "multiclass" = "multinomial",
            "survival" = "cox"
        )
    return(glmnet_family)
}

tof_finalize_model <-
    function(
        feature_tibble,
        best_model_parameters,
        recipe,
        model_type,
        outcome_colnames) {
        # find glmnet family
        glmnet_family <- tof_find_glmnet_family(model_type)

        # convert data into the format glmnet uses
        features_glmnet <-
            feature_tibble |>
            tof_setup_glmnet_xy(
                outcome_cols = dplyr::any_of(outcome_colnames),
                recipe = recipe,
                model_type = model_type
            )

        glmnet_x <- features_glmnet$x
        glmnet_y <- features_glmnet$y

        # if there is more than 1 optimal set of hyperparamters, choose the sparsest
        # model and then the one that is closest to the lasso (alpha = 1).
        best_models <-
            best_model_parameters |>
            dplyr::slice_max(order_by = .data$penalty) |>
            dplyr::slice_max(order_by = .data$mixture)

        best_alpha <- best_models$mixture
        best_lambda <- best_models$penalty

        # fit best glmnet model
        best_model <-
            glmnet::glmnet(
                x = glmnet_x,
                y = glmnet_y,
                alpha = best_alpha,
                family = glmnet_family,
                standardize = FALSE
            )

        # hack glmnet a bit to fix a bug in the predict.glmnet function
        new_call <- best_model$call
        new_call$family <- tof_find_glmnet_family(model_type)
        new_call$alpha <- best_alpha

        best_model$call <- new_call

        # construct final model
        tof_model <-
            new_tof_model(
                model = best_model,
                recipe = recipe,
                penalty = best_lambda,
                mixture = best_alpha,
                model_type = model_type,
                outcome_colnames = outcome_colnames,
                training_data = feature_tibble
            )

        return(tof_model)
    }

# accepts a feature tibble (i.e. produced by an extract function) and a
# recipe for preprocessing and converts the raw data into a preprocessed x
# matrix and y vector (or survival object for model_type == "survival")
# per glmnet's preferences
tof_setup_glmnet_xy <-
    function(feature_tibble, recipe, outcome_cols, model_type) {
        # Check that recipes is installed
        has_recipes <- requireNamespace("recipes", quietly = TRUE)
        if (!has_recipes) {
            stop(
                "This function requires the {recipes} package. Install it with this code:\n",
                "install.packages(\"recipes\")\n"
            )
        }

        # preprocess the input features using the recipe ---------------------------
        preprocessed_features <-
            recipe |>
            recipes::bake(new_data = feature_tibble)

        # separate the predictors (x) and outcomes (y) from one another ------------

        y <-
            preprocessed_features |>
            dplyr::select({{ outcome_cols }})

        predictor_colnames <-
            setdiff(colnames(preprocessed_features), colnames(y))

        x <-
            preprocessed_features |>
            dplyr::select(dplyr::any_of(predictor_colnames)) |>
            as.matrix()

        if (model_type == "survival") {
            # create a survival object for the outcome in a cox model
            y <-
                survival::Surv(
                    time = y[[1]],
                    event = y[[2]]
                )
        } else {
            y <- y[[colnames(y)]]
        }

        return(list(x = x, y = y))
    }

#' Check argument specifications for a glmnet model.
#'
#' @param split_data An `rsplit` or `rset` object from the \code{\link[rsample]{rsample}}
#' package containing the sample-level data to use for modeling. Alternatively,
#' an unsplit tbl_df can be provided, though this is not recommended.
#'
#' @param model_type A string indicating which kind of elastic net model to build.
#' If a continuous response is being predicted, use "linear" for linear regression;
#' if a categorical response with only 2 classes is being predicted, use
#' "two-class" for logistic regression; if a categorical response with more than 2
#' levels is being predicted, use "multiclass" for multinomial regression; and if
#' a time-to-event outcome is being predicted, use "survival" for Cox regression.
#'
#' @param best_model_type Currently unused.
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
#' @return A tibble. If arguments are specified correctly, this tibble can be
#' used to create a recipe for preprocessing.
#'
#'
tof_check_model_args <-
    function(
        split_data,
        model_type = c("linear", "two-class", "multiclass", "survival"),
        best_model_type = c("best", "best with sparsity"),
        response_col,
        time_col,
        event_col) {
        # Check that rsample is installed only if split_data is "rset" or "rsplit"
        if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
            has_rsample <- requireNamespace("rsample", quietly = TRUE)
            if (!has_rsample) {
                stop(
                    "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                    "Install it with:\ninstall.packages(\"rsample\")\n"
                )
            }
        }

        # check split_data
        if (inherits(split_data, "rsplit")) {
            feature_tibble <-
                split_data |>
                rsample::training()
        } else if (inherits(split_data, "rset")) {
            feature_tibble <-
                split_data$splits[[1]] |>
                rsample::training()
        } else if (inherits(split_data, "tbl_df")) {
            feature_tibble <- split_data
        } else {
            stop("split_data must be either an rsplit or an rset object.")
        }

        # check string arguments
        model_type <- rlang::arg_match(model_type)
        best_model_type <- rlang::arg_match(best_model_type)

        # check outcome variables
        if (model_type %in% c("linear", "two-class", "multiclass") & !missing(response_col)) {
            response <-
                feature_tibble |>
                dplyr::pull({{ response_col }})

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
                feature_tibble |>
                dplyr::pull({{ time_col }})

            event <-
                feature_tibble |>
                dplyr::pull({{ event_col }})

            if (!is.numeric(time)) {
                stop("All values in time_col must be numeric - they represent the time
                to event outcome for each patient (i.e. each row of of the input data.")
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

        return(feature_tibble)
    }

# Create a tibble containing only unique observations from an rset or rsplit object.
tof_all_data <- function(split_data) {
    # Check that rsample is installed only if split_data is "rset" or "rsplit"
    if (inherits(split_data, "rset") || inherits(split_data, "rsplit")) {
        has_rsample <- requireNamespace("rsample", quietly = TRUE)
        if (!has_rsample) {
            stop(
                "This function requires the {rsample} package for processing 'rset' or 'rsplit' objects. ",
                "Install it with:\ninstall.packages(\"rsample\")\n"
            )
        }
    }

    if (inherits(split_data, "rsplit")) {
        feature_tibble <-
            rsample::training(split_data) |>
            dplyr::bind_rows(rsample::testing(split_data))

        if (inherits(split_data, "bootstraps")) {
            feature_tibble <- dplyr::distinct(feature_tibble)
        }
    } else if (inherits(split_data, "rset")) {
        feature_tibble <-
            rsample::training(split_data$splits[[1]]) |>
            dplyr::bind_rows(rsample::testing(split_data$splits[[1]]))

        if (inherits(split_data, "bootstraps")) {
            feature_tibble <- dplyr::distinct(feature_tibble)
        }
    } else if (inherits(split_data, "tbl_df")) {
        feature_tibble <- split_data
    } else {
        stop("split_data must be either an rsplit or an rset object (or a tbl_df).")
    }
    return(feature_tibble)
}

# tof_model utilities ----------------------------------------------------------

#' Constructor for a tof_model.
#'
#' @param model A glmnet model.
#'
#' @param recipe A prepped recipe object.
#'
#' @param penalty A double indicating which lambda value should be used within the
#' glmnet path.
#'
#' @param mixture A double indicating which alpha value was used to fit the glmnet model.
#'
#' @param model_type A string indicating which type of glmnet model is being fit.
#'
#' @param outcome_colnames TO DO
#'
#' @param training_data TO DO
#'
#' @return A `tof_model`, an S3 class that includes a trained glmnet model and
#' the recipe used to perform its associated preprocessing.
#'
new_tof_model <-
    function(
        model,
        recipe,
        penalty,
        mixture,
        model_type = c("linear", "two-class", "multiclass", "survival"),
        outcome_colnames,
        training_data) {
        # check arguments
        stopifnot(inherits(recipe, "recipe"))
        stopifnot(inherits(model, "glmnet"))
        stopifnot(is.numeric(penalty))
        stopifnot(is.numeric(mixture))
        stopifnot(is.data.frame(training_data))
        model_type <- rlang::arg_match(model_type)

        # assemble tof_model
        tof_model <-
            list(
                model = model,
                recipe = recipe,
                mixture = mixture,
                penalty = penalty,
                model_type = model_type,
                outcome_colnames = outcome_colnames,
                training_data = training_data
            )

        # add class attribute and return
        class(tof_model) <- "tof_model"
        return(tof_model)
    }

#' @exportS3Method
print.tof_model <- function(x, ...) {
    model_type <- x$model_type
    mixture <- x$mixture
    penalty <- x$penalty

    if (model_type != "multiclass") {
        coefficients <-
            x$model |>
            glmnet::coef.glmnet(s = penalty) |>
            as.matrix() |>
            dplyr::as_tibble(rownames = "feature")

        colnames(coefficients) <- c("feature", "coefficient")

        coefficients <-
            dplyr::filter(coefficients, .data$coefficient != 0) |>
            dplyr::arrange(-abs(.data$coefficient))
    } else {
        coefficients <-
            x$model |>
            glmnet::coef.glmnet(s = penalty) |>
            purrr::map(
                .f = function(x) {
                    result <- x |>
                        as.matrix() |>
                        dplyr::as_tibble(rownames = "feature")

                    colnames(result) <- c("feature", "coefficient")

                    result <-
                        dplyr::filter(result, .data$coefficient != 0) |>
                        dplyr::arrange(-abs(.data$coefficient))

                    return(result)
                }
            )

        names(coefficients) <- levels(x$training_data[[tof_get_model_outcomes(x)]])
    }

    cat(
        "A", model_type, "`tof_model` with a mixture parameter (alpha) of",
        round(mixture, 3),
        "and a penalty parameter (lambda) of",
        format(penalty, digits = 4, scientific = TRUE), "\n"
    )

    print(coefficients)
}

#' Get a `tof_model`'s training data
#'
#'
#' @param tof_model A tof_model
#'
#' @return A tibble of (non-preprocessed) training data used to fit the model
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_training_data(regression_model)
#'
tof_get_model_training_data <-
    function(tof_model) {
        return(tof_model$training_data)
    }

#' Get a `tof_model`'s optimal penalty (lambda) value
#'
#' @param tof_model A tof_model
#'
#' @return A numeric value
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_penalty(regression_model)
#'
tof_get_model_penalty <-
    function(tof_model) {
        return(tof_model$penalty)
    }

#' Get a `tof_model`'s optimal mixture (alpha) value
#'
#' @param tof_model A tof_model
#'
#' @return A numeric value
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_mixture(regression_model)
#'
tof_get_model_mixture <-
    function(tof_model) {
        return(tof_model$mixture)
    }

#' Get a `tof_model`'s model type
#'
#' @param tof_model A tof_model
#'
#' @return A string
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_type(regression_model)
#'
tof_get_model_type <-
    function(tof_model) {
        return(tof_model$model_type)
    }

#' Get a `tof_model`'s outcome variable name(s)
#'
#' @param tof_model A tof_model
#'
#' @return A character vector
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_outcomes(regression_model)
#'
tof_get_model_outcomes <-
    function(tof_model) {
        return(tof_model$outcome_colnames)
    }

#' Get a `tof_model`'s processed predictor matrix (for glmnet)
#'
#' @param tof_model A tof_model
#'
#' @return An x value formatted for glmnet
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_x(regression_model)
#'
tof_get_model_x <-
    function(tof_model) {
        xy <-
            tof_model |>
            tof_get_model_training_data() |>
            tof_setup_glmnet_xy(
                recipe = tof_model$recipe,
                outcome_cols = tof_get_model_outcomes(tof_model),
                model_type = tof_get_model_type(tof_model)
            )

        return(xy$x)
    }

#' Get a `tof_model`'s processed outcome variable matrix (for glmnet)
#'
#' @param tof_model A tof_model
#'
#' @return A y value formatted for glmnet
#'
#' @export
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             ),
#'         multiclass =
#'             as.factor(
#'                 c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
#'             ),
#'         event = c(rep(0, times = 30), rep(1, times = 70)),
#'         time_to_event = rnorm(n = 100, mean = 10, sd = 2)
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a regression model
#' regression_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = outcome,
#'         model_type = "linear"
#'     )
#'
#' tof_get_model_y(regression_model)
#'
tof_get_model_y <-
    function(tof_model) {
        xy <-
            tof_model |>
            tof_get_model_training_data() |>
            tof_setup_glmnet_xy(
                recipe = tof_model$recipe,
                outcome_cols = tof_get_model_outcomes(tof_model),
                model_type = tof_get_model_type(tof_model)
            )

        return(xy$y)
    }

#' Access a trained elastic net model's performance metrics using its tuning data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @return A list of performance metrics whose components depend on the model type.
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_rows
#' @importFrom dplyr count
#' @importFrom dplyr everything
#' @importFrom dplyr if_else
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#' @importFrom dplyr transmute
#'
#' @importFrom glmnet assess.glmnet
#' @importFrom glmnet confusion.glmnet
#' @importFrom glmnet roc.glmnet
#'
#' @importFrom purrr discard
#'
#' @importFrom survival Surv
#' @importFrom survival survfit
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
tof_assess_model_tuning <-
    function(tof_model) {
        model_type <- tof_model$model_type

        tuning_data <-
            tof_model$tuning_metrics |>
            dplyr::select(".predictions") |>
            tidyr::unnest(cols = ".predictions")

        if (model_type %in% c("linear", "two-class", "survival")) {
            if (model_type %in% c("linear", "two-class")) {
                prediction_matrix <- as.matrix(tuning_data$response)
            } else {
                # remove infinite estimations of relative risk, which glmnet
                # can compute sometimes and that always cause bugs
                max_non_infinite <-
                    max(purrr::discard(.x = tuning_data$relative_risk, .p = is.infinite))
                tuning_data$relative_risk <-
                    dplyr::if_else(
                        is.infinite(tuning_data$relative_risk),
                        true = max_non_infinite,
                        false = tuning_data$relative_risk
                    )
                prediction_matrix <- as.matrix(tuning_data$relative_risk)
            }

            if (model_type != "survival") {
                new_y <- tuning_data$truth
            } else {
                new_y <-
                    survival::Surv(
                        time = tuning_data$true_time_to_event,
                        event = tuning_data$true_event
                    )
            }

            # if model_type is "multiclass"
        } else {
            prediction_matrix <-
                tuning_data |>
                dplyr::select(-"truth", -"class") |>
                as.matrix()

            col_names <-
                tuning_data |>
                dplyr::select(-"truth", -"class") |>
                colnames()

            prediction_matrix <-
                array(
                    data = prediction_matrix,
                    dim = c(
                        dim(prediction_matrix)[[1]],
                        dim(prediction_matrix)[[2]],
                        1L
                    )
                )

            predicted_outcome <- tuning_data$class

            new_y <- tuning_data$truth

            outcome_levels <- unique(tuning_data$truth)
            confusion_matrix <-
                tuning_data |>
                dplyr::mutate(
                    truth = factor(.data$truth, levels = outcome_levels),
                    class = factor(.data$class, levels = outcome_levels)
                ) |>
                dplyr::count(.data$truth, .data$class, .drop = FALSE) |>
                dplyr::transmute(
                    true_outcome = as.character(.data$truth),
                    predicted_outcome = as.character(.data$class),
                    num_observations = n
                )
        }

        model_metrics <-
            glmnet::assess.glmnet(
                object = prediction_matrix,
                newy = new_y,
                family = tof_find_glmnet_family(model_type)
            ) |>
            dplyr::as_tibble() |>
            tof_clean_metric_names(model_type = model_type) |>
            tidyr::pivot_longer(
                cols = dplyr::everything(),
                names_to = "metric",
                values_to = "value"
            )

        if (model_type == "two-class") {
            outcome_levels <-
                tuning_data |>
                dplyr::group_by(.data$class) |>
                dplyr::summarize(mean_prob = mean(.data$response)) |>
                dplyr::arrange(.data$mean_prob) |>
                dplyr::pull(.data$class)

            roc_curve <-
                tuning_data |>
                dplyr::mutate(truth = factor(.data$truth, levels = outcome_levels)) |>
                tof_make_roc_curve(
                    truth_col = "truth",
                    prob_cols = "response"
                )

            confusion_matrix <-
                tuning_data |>
                dplyr::mutate(
                    truth = factor(.data$truth, levels = outcome_levels),
                    class = factor(.data$class, levels = outcome_levels)
                ) |>
                dplyr::count(.data$truth, .data$class, .drop = FALSE) |>
                dplyr::transmute(
                    true_outcome = as.character(.data$truth),
                    predicted_outcome = as.character(.data$class),
                    num_observations = .data$n
                )
        } else if (model_type == "multiclass") {
            outcome_levels <- unique(tuning_data$truth)

            prediction_colnames <- paste0("prob_", outcome_levels)

            roc_auc <-
                tuning_data |>
                dplyr::mutate(truth = as.factor(.data$truth)) |>
                yardstick::roc_auc(
                    truth = "truth",
                    dplyr::any_of(prediction_colnames)
                )

            roc_auc_tibble <-
                dplyr::tibble(metric = "roc_auc", value = roc_auc$.estimate)

            model_metrics <-
                model_metrics |>
                dplyr::bind_rows(roc_auc_tibble)

            roc_curve <-
                tuning_data |>
                dplyr::rename_with(
                    cols = dplyr::everything(),
                    .fn = ~ gsub(pattern = "prob_", x = .x, replacement = "")
                ) |>
                tof_make_roc_curve(
                    truth_col = "truth",
                    prob_cols = dplyr::any_of(outcome_levels)
                )
        } else {
            roc_curve <- NULL
            if (model_type %in% c("linear", "survival")) {
                confusion_matrix <- NULL
            }
        }

        if (model_type == "survival") {
            survival_curves <-
                tuning_data |>
                dplyr::mutate(
                    risk_group =
                        dplyr::if_else(
                            .data$relative_risk > tof_model$best_log_rank_threshold,
                            "high",
                            "low"
                        )
                )

            lrt_result <-
                tof_log_rank_test(
                    input_data = survival_curves,
                    relative_risk_col = .data$relative_risk,
                    time_col = .data$true_time_to_event,
                    event_col = .data$true_event,
                    threshold = tof_model$best_log_rank_threshold
                )

            lrt_tibble <-
                dplyr::tibble(metric = "log_rank_p_value", value = lrt_result)

            model_metrics <-
                model_metrics |>
                dplyr::bind_rows(lrt_tibble)

            survival_curves <-
                survival_curves |>
                dplyr::select(
                    "survival_curve",
                    "relative_risk",
                    time_to_event = "true_time_to_event",
                    event = "true_event",
                    "risk_group"
                )
        } else {
            survival_curves <- NULL
        }

        result <-
            list(
                model_metrics = model_metrics,
                roc_curve = roc_curve,
                confusion_matrix = confusion_matrix,
                survival_curves = survival_curves
            )

        return(result)
    }


#' Compute a trained elastic net model's performance metrics using new_data.
#'
#' @param tof_model A `tof_model` trained using \code{\link{tof_train_model}}
#'
#' @param new_data A tibble of new observations that should be used to evaluate
#' the `tof_model`'s performance.
#'
#' @return A list of performance metrics whose components depend on the model type.
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr as_tibble
#' @importFrom dplyr everything
#' @importFrom dplyr tibble
#'
#' @importFrom glmnet assess.glmnet
#' @importFrom glmnet confusion.glmnet
#' @importFrom glmnet roc.glmnet
#'
#' @importFrom survival survfit
#'
#' @importFrom stats pchisq
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#'
#'
tof_assess_model_new_data <-
    function(tof_model, new_data) {
        # set up
        # extract the recipe and the glmnet model
        model <- tof_model$model
        recipe <- tof_model$recipe
        lambda <- tof_model$penalty

        # extract the model type and outcome colnames
        model_type <- tof_get_model_type(tof_model)
        outcome_colnames <- tof_get_model_outcomes(tof_model)

        # preprocess the new data ------------------------------------------------
        new_xy <-
            new_data |>
            tof_setup_glmnet_xy(
                recipe = recipe,
                outcome_cols = dplyr::any_of(outcome_colnames),
                model_type = model_type
            )

        # calculate metrics ------------------------------------------------------
        model_metrics <-
            glmnet::assess.glmnet(
                object = model,
                newx = new_xy$x,
                newy = new_xy$y,
                family = tof_find_glmnet_family(model_type),
                s = lambda
            ) |>
            dplyr::as_tibble() |>
            tof_clean_metric_names(model_type = model_type) |>
            tidyr::pivot_longer(
                cols = dplyr::everything(),
                names_to = "metric",
                values_to = "value"
            )

        if (model_type == "two-class") {
            predictions <-
                tof_predict(
                    tof_model = tof_model,
                    new_data = new_data,
                    prediction_type = "response"
                )

            roc_tibble <-
                dplyr::tibble(
                    truth = new_data[[outcome_colnames]],
                    response = predictions$.pred
                )

            roc_curve <-
                tof_make_roc_curve(
                    input_data = roc_tibble,
                    truth_col = "truth",
                    prob_cols = "response"
                )
        } else if (model_type == "multiclass") {
            predictions <-
                tof_predict(
                    tof_model = tof_model,
                    new_data = new_data,
                    prediction_type = "response"
                ) |>
                dplyr::rename_with(
                    cols = dplyr::everything(),
                    .fn = ~ gsub(pattern = ".pred_", x = .x, replacement = "")
                )
            prediction_colnames <- colnames(predictions)

            roc_auc <-
                predictions |>
                dplyr::mutate(truth = new_data[[tof_model$outcome_colnames]]) |>
                yardstick::roc_auc(
                    truth = "truth",
                    dplyr::any_of(prediction_colnames)
                )

            roc_auc_tibble <-
                dplyr::tibble(metric = "roc_auc", value = roc_auc$.estimate)

            model_metrics <-
                model_metrics |>
                dplyr::bind_rows(roc_auc_tibble)

            roc_tibble <-
                dplyr::tibble(
                    truth = new_data[[outcome_colnames]]
                ) |>
                dplyr::bind_cols(predictions)

            outcome_levels <- unique(new_data[[outcome_colnames]])

            roc_curve <-
                tof_make_roc_curve(
                    input_data = roc_tibble,
                    truth_col = "truth",
                    prob_cols = dplyr::any_of(outcome_levels)
                )
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
                ) |>
                dplyr::as_tibble() |>
                dplyr::transmute(
                    true_outcome = .data$True,
                    predicted_outcome = .data$Predicted,
                    num_observations = .data$n
                )
        } else {
            confusion_matrix <- NULL
        }

        if (model_type == "survival") {
            survival_curves <-
                tof_model_plot_survival_curves(tof_model, new_x = new_xy$x)

            predictions <-
                tof_predict(
                    tof_model = tof_model,
                    new_data = new_data,
                    prediction_type = "response"
                )

            survival_curves <-
                survival_curves |>
                dplyr::mutate(
                    relative_risk = predictions$.pred,
                    time_to_event = new_xy$y[, 1],
                    event = new_xy$y[, 2],
                    risk_group =
                        dplyr::if_else(
                            .data$relative_risk > tof_model$best_log_rank_threshold,
                            "high",
                            "low"
                        )
                )

            lrt_result <-
                tof_log_rank_test(
                    input_data = survival_curves,
                    relative_risk_col = .data$relative_risk,
                    time_col = .data$time_to_event,
                    event_col = .data$event,
                    threshold = tof_model$best_log_rank_threshold
                )

            lrt_tibble <-
                dplyr::tibble(metric = "log_rank_p_value", value = lrt_result)

            model_metrics <-
                model_metrics |>
                dplyr::bind_rows(lrt_tibble)
        } else {
            survival_curves <- NULL
        }

        result <-
            list(
                model_metrics = model_metrics,
                roc_curve = roc_curve,
                confusion_matrix = confusion_matrix,
                survival_curves = survival_curves
            )

        return(result)
    }

#' Compute a receiver-operating curve (ROC) for a two-class or multiclass dataset
#'
#' @param input_data A tof_tbl, tbl_df, or data.frame in which each row is an
#' observation.
#' @param truth_col An unquoted column name indicating which column in `input_data`
#' contains the true class labels for each observation. Must be a factor.
#' @param prob_cols Unquoted column names indicating which columns in `input_data`
#' contain the probability estimates for each class in `truth_col`. These columns
#' must be specified in the same order as the factor levels in `truth_col`.
#'
#' @return A tibble that can be used to plot the ROC for a classification task.
#' For each candidate probability threshold, the following are reported:
#' specificity, sensitivity, true-positive rate (tpr), and false-positive rate
#' (fpr).
#'
#' @export
#'
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
#' @importFrom yardstick roc_curve
#'
#' @examples
#' feature_tibble <-
#'     dplyr::tibble(
#'         sample = as.character(1:100),
#'         cd45 = runif(n = 100),
#'         pstat5 = runif(n = 100),
#'         cd34 = runif(n = 100),
#'         outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
#'         class =
#'             as.factor(
#'                 dplyr::if_else(outcome > median(outcome), "class1", "class2")
#'             )
#'     )
#'
#' split_data <- tof_split_data(feature_tibble, split_method = "simple")
#'
#' # train a logistic regression classifier
#' log_model <-
#'     tof_train_model(
#'         split_data = split_data,
#'         predictor_cols = c(cd45, pstat5, cd34),
#'         response_col = class,
#'         model_type = "two-class"
#'     )
#'
#' # make predictions
#' predictions <-
#'     tof_predict(
#'         log_model,
#'         new_data = feature_tibble,
#'         prediction_type = "response"
#'     )
#' prediction_tibble <-
#'     dplyr::tibble(
#'         truth = feature_tibble$class,
#'         prediction = predictions$.pred
#'     )
#'
#' # make ROC curve
#' tof_make_roc_curve(
#'     input_data = prediction_tibble,
#'     truth_col = truth,
#'     prob_cols = prediction
#' )
#'
tof_make_roc_curve <- function(input_data, truth_col, prob_cols) {
    outcome_levels <-
        input_data |>
        dplyr::pull({{ truth_col }}) |>
        as.character() |>
        unique()

    num_prob_cols <-
        input_data |>
        dplyr::select({{ prob_cols }}) |>
        ncol()

    if (length(outcome_levels) >= 2) {
        roc_curve <-
            input_data |>
            dplyr::mutate(
                truth = dplyr::pull(input_data, {{ truth_col }})
            ) |>
            yardstick::roc_curve({{ prob_cols }}, truth = "truth", event_level = "second") |>
            dplyr::mutate(
                tpr = .data$sensitivity,
                fpr = 1 - .data$specificity
            )
    } else {
        warning("When only one outcome level is present in new_data, an ROC cannot be computed. Returning an empty roc_curve.")
        roc_curve <- dplyr::tibble()
    }
}


#' Compute the log-rank test p-value for the difference between the two survival
#' curves obtained by splitting a dataset into a "low" and "high" risk group
#' using a given relative-risk threshold.
#'
#' @param input_data A tbl_df or data.frame in which each observation is a row.
#'
#' @param relative_risk_col An unquote column name indicating which column contains
#' the relative-risk estimates for each observation.
#'
#' @param time_col An unquoted column name indicating which column contains the
#' true time-to-event information for each observation.
#'
#' @param event_col An unquoted column name indicating which column contains the
#' outcome (event or censorship). Must be a binary column - all values should be
#'  either 0 or 1 (with 1 indicating the adverse event and 0 indicating
#'  censorship) or FALSE and TRUE (with TRUE indicating the
#' adverse event and FALSE indicating censorship).
#'
#' @param threshold A numeric value indicating the relative-risk threshold that
#' should be used to split observations into low- and high-risk groups.
#'
#' @return A numeric value <1, the p-value of the log-rank test.
#'
#' @export
#'
#' @importFrom dplyr if_else
#' @importFrom dplyr transmute
#'
#' @importFrom survival Surv
#' @importFrom survival survdiff
#'
#' @importFrom stats pchisq
#'
#' @examples
#' NULL
#'
tof_log_rank_test <-
    function(
        input_data,
        relative_risk_col,
        time_col,
        event_col,
        threshold) {
        survival_df <-
            input_data |>
            dplyr::transmute(
                relative_risk = {{ relative_risk_col }},
                time = {{ time_col }},
                event = {{ event_col }},
                strata = dplyr::if_else(.data$relative_risk > threshold, "high", "low")
            )

        if (all(survival_df$strata == "high") | all(survival_df$strata == "low")) {
            p_val <- 1
        } else {
            log_rank_result <-
                survival::survdiff(survival::Surv(time, event) ~ strata, data = survival_df)

            df <- sum(1 * log_rank_result$exp > 0) - 1
            chi_sq <- log_rank_result$chisq
            p_val <- stats::pchisq(chi_sq, df = df, lower.tail = FALSE)
        }
        return(p_val)
    }

#' Compute the log-rank test p-value for the difference between the two survival
#' curves obtained by splitting a dataset into a "low" and "high" risk group
#' using all possible relative-risk thresholds.
#'
#' @param input_data A tbl_df or data.frame in which each observation is a row.
#'
#' @param relative_risk_col An unquote column name indicating which column contains
#' the relative-risk estimates for each observation.
#'
#' @param time_col An unquoted column name indicating which column contains the
#' true time-to-event information for each observation.
#'
#' @param event_col An unquoted column name indicating which column contains the
#' outcome (event or censorship). Must be a binary column - all values should be
#'  either 0 or 1 (with 1 indicating the adverse event and 0 indicating
#'  censorship) or FALSE and TRUE (with TRUE indicating the
#' adverse event and FALSE indicating censorship).
#'
#' @return A tibble with 3 columns: "candidate_thresholds" (the relative-risk threshold
#' used for the log-rank test), "log_rank_p_val" (the p-values of the log-rank
#' tests) and "is_best" (a logical value indicating which candidate threshold gave
#' the optimal, i.e. smallest, p-value).
#'
#'
#' @importFrom dplyr pull
#' @importFrom dplyr tibble
#'
#' @importFrom purrr map_dbl
#'
tof_find_log_rank_threshold <-
    function(
        input_data,
        relative_risk_col,
        time_col,
        event_col) {
        candidate_thresholds <-
            input_data |>
            dplyr::pull({{ relative_risk_col }}) |>
            sort()

        p_values <-
            purrr::map_dbl(
                .x = candidate_thresholds,
                .f = ~ tof_log_rank_test(
                    input_data = input_data,
                    relative_risk_col = {{ relative_risk_col }},
                    time_col = {{ time_col }},
                    event_col = {{ event_col }},
                    threshold = .x
                )
            )

        threshold_index <- which.min(p_values)
        best_threshold <- candidate_thresholds[[threshold_index]]

        result <-
            dplyr::tibble(
                candidate_thresholds = candidate_thresholds,
                log_rank_p_val = p_values,
                is_best = candidate_thresholds == best_threshold
            )
        return(result)
    }
