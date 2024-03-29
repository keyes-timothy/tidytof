% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patient-level_modeling.R
\name{tof_train_model}
\alias{tof_train_model}
\title{Train an elastic net model to predict sample-level phenomena using high-dimensional cytometry data.}
\usage{
tof_train_model(
  split_data,
  unsplit_data,
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
)
}
\arguments{
\item{split_data}{An `rsplit` or `rset` object from the \code{\link[rsample]{rsample}}
package containing the sample-level data to use for modeling.
The easiest way to generate this is to use \code{\link{tof_split_data}}.}

\item{unsplit_data}{A tibble containing sample-level data to use for modeling
without resampling. While using a resampling method is advised, this argument
provides an interface to fit a model without using cross-validation or
bootstrap resampling. Ignored if split_data is provided.}

\item{predictor_cols}{Unquoted column names indicating which columns in the
data contained in `split_data` should be used as predictors in the elastic net model.
Supports tidyselect helpers.}

\item{response_col}{Unquoted column name indicating which column in the data
contained in `split_data` should be used as the outcome in a "two-class", "multiclass",
or "linear" elastic net model. Must be a factor for "two-class" and "multiclass"
models and must be a numeric for "linear" models. Ignored if `model_type` is "survival".}

\item{time_col}{Unquoted column name indicating which column in the data
contained in `split_data` represents the time-to-event outcome in a "survival"
elastic net model. Must be numeric. Ignored if `model_type` is "two-class", "multiclass",
or "linear".}

\item{event_col}{Unquoted column name indicating which column in the data
contained in `split_data` represents the time-to-event outcome in a "survival"
elastic net model. Must be a binary column - all values should be either 0 or 1
(with 1 indicating the adverse event) or FALSE and TRUE (with TRUE indicating the
adverse event). Ignored if `model_type` is "two-class", "multiclass",
or "linear".}

\item{model_type}{A string indicating which kind of elastic net model to build.
If a continuous response is being predicted, use "linear" for linear regression;
if a categorical response with only 2 classes is being predicted, use
"two-class" for logistic regression; if a categorical response with more than 2
levels is being predicted, use "multiclass" for multinomial regression; and if
a time-to-event outcome is being predicted, use "survival" for Cox regression.}

\item{hyperparameter_grid}{A hyperparameter grid indicating which values of
the elastic net penalty (lambda) and the elastic net mixture (alpha) hyperparamters
should be used during model tuning. Generate this grid using \code{\link{tof_create_grid}}.}

\item{standardize_predictors}{A logical value indicating if numeric predictor columns
should be standardized (centered and scaled) before model fitting, as is
standard practice during elastic net regularization. Defaults to TRUE.}

\item{remove_zv_predictors}{A logical value indicating if predictor columns
with near-zero variance should be removed before model fitting using
\code{\link[recipes]{step_nzv}}. Defaults to FALSE.}

\item{impute_missing_predictors}{A logical value indicating if predictor columns
should have missing values imputed using k-nearest neighbors before model fitting (see
\code{\link[recipes]{step_impute_knn}}). Imputation is performed using an observation's
5 nearest-neighbors. Defaults to FALSE.}

\item{optimization_metric}{A string indicating which optimization metric
should be used for hyperparameter selection during model tuning. Valid values
depend on the model_type.

\itemize{
\item For "linear" models, choices are "mse" (the mean squared error
of the predictions; the default) and "mae" (the mean absolute error of the predictions).
\item For "two-class" models, choices are "roc_auc" (the area under the Receiver-Operating
Curve for the classification; the default), "misclassification error" (the proportion of
misclassified observations), "binomial_deviance" (see \code{\link[glmnet]{deviance.glmnet}}),
"mse" (the mean squared error of the logit function), and "mae" (the mean absolute error of the
logit function).
\item For "multiclass" models, choices are "roc_auc" (the area under the Receiver-Operating
Curve for the classification using the Hand-Till generalization of the ROC AUC
for multiclass models in \code{\link[yardstick]{roc_auc}}; the default), "misclassification error"
(the proportion of
misclassified observations), "multinomial_deviance" (see \code{\link[glmnet]{deviance.glmnet}}),
and "mse" and "mae" as above.
\item For "survival" models, choices are "concordance_index" (Harrel's C index;
see \code{\link[glmnet]{deviance.glmnet}}) and "partial_likelihood_deviance"
(see \code{\link[glmnet]{deviance.glmnet}}).
}}

\item{best_model_type}{Currently unused.}

\item{num_cores}{Integer indicating how many cores should be used for parallel
processing when fitting multiple models. Defaults to 1. Overhead to separate
models across multiple cores can be high, so significant speedup is unlikely
to be observed unless many large models are being fit.}
}
\value{
A `tof_model`, an S3 class that includes the elastic net model with
the best performance (assessed via cross-validation, bootstrapping, or simple splitting
depending on `split_data`) across all tested hyperparameter value combinations. `tof_models`
store the following information:

\describe{
\item{model}{The final elastic net ("glmnet") model, which is chosen by selecting
the elastic net hyperparameters with the best `optimization_metric` performance
on the validation sets of each resample used to train the model (on average)}
\item{recipe}{The \code{\link[recipes]{recipe}} used for data preprocessing}
\item{mixture}{The optimal mixture hyperparameter (alpha) for the glmnet model}
\item{penalty}{The optimal penalty hyperparameter (lambda) for the glmnet model}
\item{model_type}{A string indicating which type of glmnet model was fit}
\item{outcome_colnames}{A character vector representing the names of the columns in the training data modeled as outcome variables}
\item{training_data}{A tibble containing the (not preprocessed) data used to train the model}
\item{tuning_metrics}{A tibble containing the validation set performance metrics
(and model predictions) during for each resample fold during model tuning.}
\item{log_rank_thresholds}{For survival models only, a tibble containing information about the relative-risk
thresholds that can be used to split the training data into 2 risk groups (low-
and high-risk) based on the final model's predictions. For each relative-risk
threshold, the log-rank test p-value and an indicator of which threshold gives
the most significant separation is provided.}
\item{best_log_rank_threshold}{For survival models only, a numeric value
representing the relative-risk threshold that yields the most significant
log-rank test when separating the training data into low- and high-risk groups.}
}
}
\description{
This function uses a training set/test set paradigm to tune and fit an
elastic net model using a variety of user-specified details. Tuning can be
performed using either a simple training vs. test set split, k-fold cross-validation,
or bootstrapping, and multiple preprocessing options are available.
}
\examples{
feature_tibble <-
    dplyr::tibble(
        sample = as.character(1:100),
        cd45 = runif(n = 100),
        pstat5 = runif(n = 100),
        cd34 = runif(n = 100),
        outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
        class =
            as.factor(
                dplyr::if_else(outcome > median(outcome), "class1", "class2")
            ),
        multiclass =
            as.factor(
                c(rep("class1", 30), rep("class2", 30), rep("class3", 40))
            ),
        event = c(rep(0, times = 30), rep(1, times = 70)),
        time_to_event = rnorm(n = 100, mean = 10, sd = 2)
    )

split_data <- tof_split_data(feature_tibble, split_method = "simple")

# train a regression model
tof_train_model(
    split_data = split_data,
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear"
)

# train a logistic regression classifier
tof_train_model(
    split_data = split_data,
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class"
)

# train a cox regression survival model
tof_train_model(
    split_data = split_data,
    predictor_cols = c(cd45, pstat5, cd34),
    time_col = time_to_event,
    event_col = event,
    model_type = "survival"
)

}
\seealso{
Other modeling functions: 
\code{\link{tof_assess_model}()},
\code{\link{tof_create_grid}()},
\code{\link{tof_predict}()},
\code{\link{tof_split_data}()}
}
\concept{modeling functions}
