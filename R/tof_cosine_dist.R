###   tof_cosine_dist.R

# Description: A function that finds the cosine distance between a given vector and each column 
#              of a matrix. 
###  
#
# Inputs: 
#     - matrix = a matrix or matrix-like object. Rows should represent independent observations (cells) 
#                and columns should represent features (measured variables). 
#     - vector = a vector or vector-like object from which cosine distances to each row in `matrix` 
#                should be calculated. 
#
# Outputs: 
#     - distances = a vector of distances from each cell in `matrix` to `vector`
#
# Dependencies: 
#     - tidyverse library
tof_cosine_dist <- function(matrix, vector) {
  diag_matrix_matrix_t <- rep(0, nrow(matrix))
  for (i in 1:nrow(matrix)) { 
    diag_matrix_matrix_t[[i]] <- crossprod(matrix[i , ])
  }
  distances <- 
    (matrix %*% vector) / 
    (as.vector(sqrt(crossprod(vector))) * sqrt(diag_matrix_matrix_t))
  
  distances <- 1 - distances
  
  return(distances)
}
