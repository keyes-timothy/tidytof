###    tof_classifier_apply.R 
# Description: 
### performs standard mass cytometry transformations on a tof_tibble.  
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or somehing that can be coerced into either 
#       of these
#     - classifier_fit = a classifier_fit object, i.e. the output from `tof_clasifier_build()`
#     - num_cores = number of cores to use for parallel processing (default is 1)
#     - parallel_var = an unquoted variable name indicating which variable should be used to split 
#                      the data for parallel processing.
#     - dist_fun = function for determining the distance between each cell and each cluster centroid
#                  (default is mahalanobis)
#
# Outputs: 
#     - classification_data = a tibble in which each of the columns represents the distance from 
#                             a given cell (each row) to a population in the `classifier_fit`. Will have
#                             the same number of rows as the original `tof_tibble`.
#
# Dependencies: 
#     - tidyverse library
#     - foreach library
#     - classify_cell function

tof_classifier_apply <- function(
  tof_tibble = NULL,
  classifier_fit = NULL, 
  num_cores = 1, 
  parallel_var = file_names, 
  dist_fun = "mahalanobis"
) {
  
  classifier_markers <- colnames(classifier_fit$covariance_matrix[[1]])
  `%my_do%` <- `%do%`
  my_cluster <- NULL
  
  # nest cancer data
  tof_tibble <- 
    tof_tibble %>% 
    group_by({{parallel_var}}) %>% 
    select(one_of(classifier_markers)) %>% 
    nest() %>% 
    ungroup()
  
  # set up parallel back-end
  if (num_cores != 1) {
    my_cluster <- makeCluster(num_cores)
    registerDoParallel(my_cluster)
    `%my_do%` <- `%dopar%`
  }
  
  classification_data <- 
    foreach(
      expression_matrix = tof_tibble$data,
      .combine = list, 
      .packages = c("dplyr", "purrr", "tidyr"), 
      .export = c("tof_classify_cell", "classifier_fit", "tof_cosine_dist"), 
      .multicombine = TRUE, 
      .maxcombine = nrow(tof_tibble)
    ) %my_do%
    tof_classify_cell(
      classifier_fit = classifier_fit, 
      expression_matrix = expression_matrix, 
      dist_fun = dist_fun
    )

  # stop cluster if it was set up
  if (!is.null(my_cluster)) {
    stopCluster(my_cluster)
  }
  
  classification_data <- 
    tibble(classification_data = classification_data) %>%
    unnest(cols = classification_data) %>% 
    rename_with(function(x) str_c(dist_fun, x, sep = "_"), .cols = everything())
  
  return(classification_data)
}
