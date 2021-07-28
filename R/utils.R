
# tidytof_example_data ------------------

#' Get paths to tidytof example data
#'
#' tidytof comes bundled with a number of sample .fcs files in its
#' inst/extdata directory. This function makes them easy to access.
#'
#' @param dataset_name Name of the dataset you want to access. If NULL,
#' the names of the datasets (each of which is from a different study)
#' will be listed.
#'
#' @return A character vector of file paths where the requested .fcs
#' files are located. If `dataset_name` is NULL, a character vector of
#' dataset names (that can be used as values for `dataset_name`) is
#' returned instead.
#'
#' @export
#'
#' @examples
#' tidytof_example_data()
#' tidytof_example_data(dataset_name = "phenograph")
#'
tidytof_example_data <-
  function(dataset_name = NULL) {

    if (is.null(dataset_name)) {
      dir(system.file("extdata", package = "tidytof"))
    }
    else {
      dir(
        system.file("extdata", dataset_name, package = "tidytof", mustWork = TRUE),
        full.names = TRUE
      )
    }
  }


# get_extension ------------------

#' Find the extension for a file
#'
#' @param filename A string representing the name of a file in its local directory
#'
#' @return The the file extension of `filename`
#'
#' @examples
#' \dontrun{
#' # example file name
#' my_filename <- "my_file.txt"
#'
#' # find and print the extension
#' my_extension <- getExtension(my_filename)
#' print(my_extension)
#' }
get_extension <- function(filename) {
  ex <- strsplit(basename(filename), split="\\.")[[1]]
  return(ex[[length(ex)]])
}


#' Reverses arcsinh transformation with cofactor `scale_factor` and a shift of `shift_factor`.
#'
#' @param x A numeric vector.
#'
#' @param shift_factor The scalar value `a` in the following equation used to
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function:
#'    `new_x <- asinh(a + b * x)`.
#'
#' @param scale_factor The scalar value `b` in the following equation used to
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function:
#'    `new_x <- asinh(a + b * x)`.
#'
#' @return A numeric vector after undergoing reverse
#' arcsinh transformation
#'
#' @export
#'
#'
rev_asinh <- function(x, shift_factor, scale_factor) {

  new_x <- (sinh(x) - shift_factor) / scale_factor
  return(new_x)

}

#' Find if a vector is numeric
#'
#' This function takes an input vector `.vec` and checks if it is either an
#' integer or a double (i.e. is the type of vector that might encode CyTOF
#' measurements).
#'
#' @param .vec
#'
#' @return A boolean value indicating if .vec is of type integer or double.
#'
#' @examples
#' NULL
tof_is_numeric <- function(.vec) {
  return(purrr::is_integer(.vec) || purrr::is_double(.vec))
}

tof_knn_density <-
  function(
    neighbor_ids, # an N by K matrix representing each cell's knn IDs
    neighbor_distances, # an N by K matrix representing each cell's knn distances
    method = c("mean_distance", "sum_distance", "xshift"),
    d # optional argument, only if method = "xshift", represents the number of dimensions over which distances were calculated
  ) {
    # check method argument
    method <-
      match.arg(method, choices = c("mean_distance", "sum_distance", "xshift"))

    # extract needed values
    k <- ncol(neighbor_ids)
    n <- nrow(neighbor_ids)

    # find densities using one of 3 methods
    if (method == "mean_distance") {
      densities <- base::colSums(abs(neighbor_distances))
    } else if (method == "sum_distance") {
      densities <- base::colMeans(abs(neighbor_distances))
    } else if (method == "xshift") {
      # transform distances into actual cosine space
      neighbor_distances <- neighbor_distances

      # find longest distance for every row (each cell)
      largest_dist <- apply(X = neighbor_distances, MARGIN = 1, FUN = max)

      #
      densities <-
        # need to figure out what compile.dists means...
        (1/(n*(largest_dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d

    } else {
      stop("Not a valid method.")
    }

    # normalize densities?
    densities <-
      (densities - min(densities)) /
      ((max(densities) - min(densities)))

    return(densities) # a vector of length N (number of cells) with the ith
    # entry representing the KNN-estimated density of the ith cell.
  }



