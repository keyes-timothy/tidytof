###   tof_classifier_build.R

# Description: 
### performs standard mass cytometry transformations on a tof_tibble.  
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#                   of these. Must contain information about each healthy population that 
#                   the classifier will be able to classify into. 
#     - population_vector = a character vector in which each entry represents the healthy 
#                           cell population that each row in  `tof_tibble` belongs to. 
#     - classifier_markers = a character vector of channel names that should be included in the classifier.
#
# Outputs: 
#     - classifier_fit = a tibble with three columns: 
#                           * population = name of the healthy cell population 
#                           * centroid = the centroid vector for that cell population
#                           * covariance_matrix = the covariance matrix for that cell population
#
# Dependencies: 
#     - tidyverse library

tof_classifier_build <- function(
  tof_tibble = NULL, 
  population_vector = NULL, 
  classifier_markers = NULL 
) { 
  message("Building classifier from healthy cells...")

  message("Reshaping healthy data...")
  tof_tibble <- 
    tof_tibble %>% 
    mutate(population = population_vector) %>% 
    select(one_of(classifier_markers), population) %>% 
    group_by(population) %>% 
    nest()
  
  # Calculate mean and covariance matrix for each population
  message("Calculating mean and covariance matrix for all healthy populations")
  classifier_fit <- 
    tof_tibble %>% 
    transmute(
      population,
      centroid = 
        map(
          data, 
          ~ summarize(.x, across(everything(), mean)) %>% 
            pivot_longer(everything()) %>% 
            deframe()
        ), 
      covariance_matrix = map(data, cov)
    )
  message("Done! Returning classifier_fit object")
  return(classifier_fit)
}
