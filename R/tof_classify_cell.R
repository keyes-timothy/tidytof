# Description:
### A helper function to assign cells to populations using the Davis Lab classifier.
#
# Inputs:
#     - classifier_fit = a classifier_fit object, i.e. the output from `tof_clasifier_build()`
#     - expression_matrix = ____
#     - dist_fun = a string indicating which distance metric should be used to calculated the distances
#                  between a centroid and each cell in `expression_matrix`. Options are as the following:
#                       * mahalanobis distance
#                       * cosine distance (also called "cosine similarity")
#                       * pearson correlation distance
#
# Outputs:
#     -
#
# Dependencies:
#     - tidyverse library
#     - tof_cosine_dist()
#

tof_classify_cell <-
  function(
    classifier_fit,
    expression_matrix,
    dist_fun = "mahalanobis"
  ){

    # Calculate distance to each population
    if (dist_fun == "mahalanobis") {
      classifications <-
        classifier_fit %>%
        transmute(
          classification_data =
            map2(
              .x = centroid,
              .y = covariance_matrix,
              .f = ~ mahalanobis(expression_matrix, center = .x, cov = .y)
            )
        ) %>%
        pull(classification_data) %>%
        quietly(bind_cols)() %>%
        pluck("result")

    } else if (dist_fun == "cosine") {
      classifications <-
        classifier_fit %>%
        transmute(
          classification_data =
            map(
              .x = centroid,
              .f =
                ~ tof_cosine_dist(
                  matrix = as.matrix(expression_matrix),
                  vector = .x
                )
            )
        ) %>%
        pull(classification_data) %>%
        quietly(bind_cols)() %>%
        pluck("result")

    } else if (dist_fun == "pearson") {
      classifications <-
        classifier_fit %>%
        mutate(
          classification_data =
            map(
              .x = centroid,
              .f = ~ 1 - cor(x = t(as.matrix(expression_matrix)), y = .x)
            )
        ) %>%
        pull(classification_data) %>%
        quietly(bind_cols)() %>%
        pluck("result")
    }

    # Is there a way to set colnames directly in the piped call above?
    colnames(classifications) <- classifier_fit$population

    #This has to be optimized
    population_names <-
      classifications %>%
      mutate(cell_id = 1:nrow(classifications)) %>%
      pivot_longer(
        cols = -cell_id,
        names_to = "cluster",
        values_to = "distance"
      ) %>%
      group_by(cell_id) %>%
      filter(distance == min(distance)) %>%
      select(-distance)

    classifications <-
      classifications %>%
      mutate(cell_id = 1:nrow(classifications)) %>%
      left_join(population_names, by = "cell_id") %>%
      select(-cell_id)

    return(classifications)
  }

