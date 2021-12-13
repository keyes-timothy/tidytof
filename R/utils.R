
# tidytof_example_data ---------------------------------------------------------

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
      system.file("extdata", dataset_name, package = "tidytof", mustWork = TRUE)#,
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
#' @param .vec A vector.
#'
#' @return A boolean value indicating if .vec is of type integer or double.
#'
#' @importFrom purrr is_integer
#' @importFrom purrr is_double
#'
#'
tof_is_numeric <- function(.vec) {
  return(purrr::is_integer(.vec) | purrr::is_double(.vec))
}


#' Find the k-nearest neighbors of each cell in a CyTOF dataset.
#'
#' @param .data A `tof_tibble` or `tibble` in which each row represents a cell
#' and each column represents a CyTOF measurement.
#'
#' @param k An integer indicating the number of nearest neighbors to return for
#' each cell.
#'
#' @param distance_function A string indicating which distance function to use
#' for the nearest-neighbor calculation. Options include "euclidean"
#' (the default) and "cosine" distances.
#'
#' @param ... Optional additional arguments to pass to RANN::nn2
#'
#' @return A list with two elements: "neighbor_ids" and "neighbor_distances,"
#' both of which are n by k matrices (in which n is the number of cells in the
#' input `.data`. The [i,j]-th entry of "neighbor_ids" represents the row index
#' for the j-th nearest neighbor of the cell in the i-th row of `.data`.
#' The [i,j]-th entry of "neighbor_distances" represents the distance between
#' those two cells according to `distance_function`.
#'
#' @export
#'
#'
tof_find_knn <-
  function(
    .data,
    k = min(10, nrow(.data)),
    distance_function = c("euclidean", "cosine"),
    ...
  ) {
    # check distance function
    distance_function <-
      match.arg(distance_function, choices = c("euclidean", "cosine"))

    # if nns in cosine distance wanted, l2 normalize all rows, as the euclidean
    # distance between l2-normalized vectors will scale with
    # cosine distance (and this allows us to us to proceed as normal)
    if (distance_function == "cosine") {
      .data <- t(apply(X = .data, MARGIN = 1, FUN = l2_normalize))
    }

    # compute result
    # have found that eps up to 0.4 will provide NN accuracy above 95%
    nn_result <- RANN::nn2(data = .data, k = k + 1, ...)
    names(nn_result) <- c("neighbor_ids", "neighbor_distances")

    # remove the first-closest neighbor (column), which is always the point itself
    nn_result$neighbor_ids <- nn_result$neighbor_ids[, 2:(k + 1)]
    nn_result$neighbor_distances <- nn_result$neighbor_distances[, 2:(k + 1)]

    # return result
    return(nn_result)
  }

#' Find the KNN density estimate for each cell in a CyTOF dataset.
#'
#' @param neighbor_ids A n by k matrix returned by `tof_find_knn` representing the row indices
#' of the k nearest neighbors of each of the n cells in a CyTOF dataset.
#'
#' @param neighbor_distances A n by k matrix returned by `tof_find_knn` representing the pairwise distances
#' between a cell and each of its k nearest neighbors in a CyTOF dataset.
#'
#' @param method A string indicating how the relative density for each cell should be
#' calculated from the distances between it and each of its k nearest neighbors. Options are
#' "mean_distance" (the default; estimates the relative density for a cell's neighborhood by
#' taking the negative average of the distances to its nearest neighbors) and "sum_distance"
#' (estimates the relative density for a cell's neighborhood by taking the negative sum of the
#' distances to its nearest neighbors).
#'
#' @param normalize TO DO
#'
#' @return a vector of length N (number of cells) with the ith
# entry representing the KNN-estimated density of the ith cell.
#'
#' @export
#'
#' @examples
#' NULL
tof_knn_density <-
  function(
    neighbor_ids, # an N by K matrix representing each of N cell's knn IDs
    neighbor_distances, # an N by K matrix representing each of N cell's knn distances
    method = c("mean_distance", "sum_distance"),
    normalize = TRUE
  ) {

    # check method argument
    method <-
      match.arg(method, choices = c("mean_distance", "sum_distance"))

    # extract needed values
    k <- ncol(neighbor_ids)
    n <- nrow(neighbor_ids)

    # find densities using one of 3 methods
    if (method == "sum_distance") {
      densities <- -base::rowSums(abs(neighbor_distances))
    } else if (method == "mean_distance") {
      densities <- -base::rowMeans(abs(neighbor_distances))
    } else {
    stop("Not a valid method.")
    }

    if (normalize) {
      # normalize densities
      densities <-
        (densities - min(densities)) /
        ((max(densities) - min(densities)))
    }

    return(densities) # a vector of length N (number of cells) with the ith
    # entry representing the KNN-estimated density of the ith cell.
  }

make_binary_vector <- function(length, indices) {
  result <- rep.int(0, times = length)
  result[indices] <- 1

  return(result)
}


#' Find the dot product between two vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return The dot product between x and y.
#'
#'
#' @examples
#' NULL
#'
dot <- function(x, y) {
  return(as.numeric(t(x) %*% y))
}

#' Find the magnitude of a vector.
#'
#' @param x A numeric vector.
#'
#' @return A scalar value (the magnitude of x).
#'
#' @examples
#' NULL
magnitude <- function(x) {
  return(sqrt(dot(x, x)))
}

l2_normalize <- function(x) {
  return(x / magnitude(x))
}

cosine_similarity <- function(x, y) {
  result <- dot(x, y) / (magnitude(x) * magnitude(y))
  return(result)
}



#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#' @importFrom purrr map
#' @importFrom tidyr nest
#'
prepare_diffcyt_args <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    fixed_effect_cols,
    random_effect_cols,
    method = c("glmm", "edgeR", "voom"),
    include_observation_level_random_effects = FALSE
  ) {
    # initialize formula
    my_formula <- NULL

    # extract sample column name as a string
    sample_colname <-
      rlang::enquo(sample_col) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
      names()

    # extract cluster column name as a string
    cluster_colname <-
      rlang::enquo(cluster_col) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
      names()

    # extract fixed effect columns as a character vector
    # will return an empty character vector if the argument is missing
    fixed_effect_colnames <-
      rlang::enquo(fixed_effect_cols) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
      names()

    # extract random effect columns as a character vector
    # will return an empty character vector if the argument is missing
    random_effect_colnames <-
      rlang::enquo(random_effect_cols) %>%
      tidyselect::eval_select(expr = ., data = tof_tibble) %>%
      names()

    # find all marker names by process of elimination
    marker_names <-
      colnames(tof_tibble)[
        !(
          colnames(tof_tibble) %in%
            c(cluster_colname, sample_colname, fixed_effect_colnames, random_effect_colnames)
        )
      ]

    if (length(marker_names) < 2) {
      stop("At least 2 markers must be selected (due to diffcyt's underlying implementation).")
    }

    # create diffcyt experiment_info
    experiment_info <-
      tof_tibble %>%
      dplyr::select({{sample_col}}, {{fixed_effect_cols}}, {{random_effect_cols}}) %>%
      dplyr::distinct() %>%
      dplyr::rename(sample_id = {{sample_col}}) %>%
      dplyr::arrange(sample_id)

    # create diffcyt marker_info
    marker_info <-
      dplyr::tibble(marker_name = marker_names) %>%
      dplyr::mutate(marker_class = "state")

    # create formula or design matrix depending on which method is being used
    if (method %in% c("glmm", "lmm") & include_observation_level_random_effects) {
      random_effect_colnames <-
        c("sample_id", random_effect_colnames)
    }

    if (method %in% c("glmm", "lmm")) {
      # if using glmms, create formula

      if (length(random_effect_colnames) == 0) {
        # if there are no random effects
        my_formula <-
          diffcyt::createFormula(
            experiment_info = experiment_info,
            cols_fixed = fixed_effect_colnames
          )
      } else {
        # if there are random effects
        my_formula <-
          diffcyt::createFormula(
            experiment_info = experiment_info,
            cols_fixed = fixed_effect_colnames,
            cols_random = random_effect_colnames
          )
      }
    }

    # create design matrix
    my_design <-
      diffcyt::createDesignMatrix(
        experiment_info = experiment_info,
        cols_design = fixed_effect_colnames
      )

    ## make contrast matrix
    contrast_names <- colnames(my_design)

    # tibble of contrast matrices to test the null hypothesis that any given
    # fixed-effect coefficient is 0
    contrast_matrix_tibble <-
      dplyr::tibble(
        contrast_names = contrast_names,
        contrast_matrices =
          purrr::map(
            .x = 1:length(contrast_names),
            .f = ~
              make_binary_vector(length = length(contrast_names), indices = .x) %>%
              diffcyt::createContrast()
          )
      ) %>%
      dplyr::filter(contrast_names != "(Intercept)")

    # test against the null hypothesis that all fixed-effect coefficients are 0
    initial_contrast <-
      diffcyt::createContrast(
        make_binary_vector(length = length(contrast_names), indices = -1)
      )

    contrast_matrix_tibble <-
      dplyr::bind_rows(
        dplyr::tibble(contrast_names = "omnibus", contrast_matrices = list(initial_contrast)),
        contrast_matrix_tibble
      )

    # configure data into the format diffcyt likes
    data_list <-
      tof_tibble %>%
      dplyr::group_by({{sample_col}}) %>%
      tidyr::nest() %>%
      dplyr::arrange({{sample_col}}) %>%
      dplyr::pull(data)

    cols_to_include <-
      colnames(data_list[[1]]) %in%
      marker_info$marker_name

    data_diff <-
      diffcyt::prepareData(
        d_input = data_list,
        experiment_info = as.data.frame(experiment_info),
        marker_info = as.data.frame(marker_info),
        cols_to_include = cols_to_include
      )

    # add clusters to diffcyt object
    temp <-
      data_diff %>%
      SummarizedExperiment::rowData()

    temp[,"cluster_id"] <-
      tof_tibble %>%
      dplyr::pull({{cluster_col}}) %>%
      as.factor()

    SummarizedExperiment::rowData(data_diff) <- temp

    # fix type-related issues in the exprs component of the SummarizedExperiment
    data_exprs <-
      data_diff %>%
      SummarizedExperiment::assays() %>%
      `[[`("exprs")

    data_colnames <- colnames(data_exprs)

    data_exprs <-
      data_exprs %>%
      apply(MARGIN = 2, FUN = as.numeric)

    colnames(data_exprs) <- data_colnames

    SummarizedExperiment::assays(data_diff)[["exprs"]] <- data_exprs

    # return result
    diffcyt_args <-
      list(
        sample_colname = sample_colname,
        cluster_colname = cluster_colname,
        fixed_effect_colnames = fixed_effect_colnames,
        random_effect_colnames = random_effect_colnames,
        marker_names = marker_names,
        experiment_info = experiment_info,
        marker_info = marker_info,
        my_formula = my_formula,
        my_design = my_design,
        contrast_matrix_tibble = contrast_matrix_tibble,
        data_diff = data_diff
      )

    return(diffcyt_args)
  }



#'
#' @importFrom lme4 glmer
#' @importFrom stats glm
#'
fit_da_model <- function(data, formula, has_random_effects = TRUE) {

  if (has_random_effects) {
    model_fit <-
      lme4::glmer(formula, data, family = "binomial", weights = total_cells)
  } else {
    model_fit <-
      stats::glm(formula, data, family = "binomial", weights = total_cells)
  }

}

#'
#' @importFrom lmerTest lmer
#' @importFrom stats glm
#'
fit_de_model <- function(data, formula, has_random_effects = TRUE) {
  if (has_random_effects) {
    model_fit <-
      lmerTest::lmer(formula, data)
  } else {
      model_fit <-
        stats::glm(formula, data, family = "gaussian")
  }
}



#' Find the earth-mover's distance between two numeric vectors
#'
#' @param vec_1 A numeric vector.
#'
#' @param vec_2 A numeric vector.
#'
#' @param num_bins An integer number of bins to use when performing kernel
#' density estimation on the two vectors. Defaults to 100.
#'
#' @return A double (of length 1) representing the EMD between the two vectors.
#'
#' @importFrom emdist emd2d
#'
#'
tof_find_emd <- function(vec_1, vec_2, num_bins = 100) {

  # check that both vec_1 and vec_2 are numeric
  if (!is.numeric(vec_1) | ! is.numeric(vec_2)) {
    return(NA_real_)
  }

  if (length(vec_1) < 10 | length(vec_2) < 10) {
    return(NA_real_)
  }

  # set up bins
  set_max <- ceiling(max(vec_1, vec_2))
  set_min <- floor(min(vec_1, vec_2))
  bins <- seq(set_max, set_min, length.out = num_bins)

  # find histograms
  density_1 <-
    graphics::hist(vec_1, breaks = bins, plot = FALSE)$density %>%
    as.matrix()

  density_2 <-
    graphics::hist(vec_2, breaks = bins, plot = FALSE)$density %>%
    as.matrix()

  # return final result
  em_dist <- emdist::emd2d(density_1, density_2)
  return(em_dist)
}


#' Find the Jensen-Shannon Divergence (JSD) between two numeric vectors
#'
#' @param vec_1 A numeric vector.
#'
#' @param vec_2 A numeric vector.
#'
#' @param num_bins An integer number of bins to use when binning
#' across the two vectors' commbined range. Defaults to 100.
#'
#' @return A double (of length 1) representing the JSD between the two vectors.
#'
#' @importFrom philentropy JSD
#' @importFrom purrr quietly
#'
#'
tof_find_jsd <- function(vec_1, vec_2, num_bins = 100) {

  # check that both vectors are numeric
  if (!is.numeric(vec_1) | ! is.numeric(vec_2)) {
    return(NA_real_)
  }

  # check that both vectors have at least 10 cells
  if (length(vec_1) < 10 | length(vec_2) < 10) {
    return(NA_real_)
  }

  # set up bins
  set_max <- ceiling(max(vec_1, vec_2))
  set_min <- floor(min(vec_1, vec_2))
  bins <- seq(set_max, set_min, length.out = num_bins)

  # find histograms
  p1 <-
    graphics::hist(vec_1, breaks = bins, plot = FALSE)$counts %>%
    as.numeric()

  p1 <- p1 / sum(p1)

  p2 <-
    graphics::hist(vec_2, breaks = bins, plot = FALSE)$counts %>%
    as.numeric()

  p2 <- p2 / sum(p2)

  # return final result
  js_dist <-
    purrr::quietly(philentropy::JSD)(rbind(p1, p2))$result %>%
    as.numeric()

  return(js_dist)
}




pull_unless_null <- function(tib, uq_colname) {
  if (is.null(tib)) {
    return(NULL)
  } else {
    return(dplyr::pull(tib, {{uq_colname}}))
  }
}

# because rsample::nested_cv is a bit buggy when not used interactively
# (call handling is a bit weird in the inside/outside arguments)
#' @importFrom rsample vfold_cv
#' @importFrom rsample training
#' @importFrom rsample bootstraps
#' @importFrom purrr map
tof_nested_cv <-
  function(
    .data,
    strata,
    inner_split_method = c("k-fold", "bootstrap"),
    num_folds = 10
  ) {
    # perform outer resampling
    outside <-
      rsample::vfold_cv(data = .data, v = num_folds, strata = {{strata}})

    # perform inner resampling
    if (inner_split_method == "k-fold") {
      inner_resample_fn <- function(split, num_folds, strata) {
        inner_sample_result <-
          split %>%
          rsample::training() %>%
          rsample::vfold_cv(v = num_folds, strata = {{strata}})

        return(inner_sample_result)
      }

      inner_attr <- paste0(num_folds - 1, "-fold cross-validation")

    } else {
      inner_resample_fn <- function(split, num_folds, strata) {
        inner_sample_result <-
          split %>%
          rsample::training() %>%
          rsample::bootstraps(times = num_folds, strata = {{strata}})

        return(inner_sample_result)
      }

      inner_attr <- paste0("Bootstrapped resampling (", num_folds - 1, " times)")

    }

    inside <- purrr::map(
      .x = outside$splits,
      .f = inner_resample_fn,
      num_folds = num_folds - 1,
      strata = {{strata}}
    )

    split_result <-
      outside %>%
      dplyr::mutate(inner_resamples = inside)

    # add required class information
    attr(split_result, "outside") <- paste0(num_folds, "-fold cross-validation")
    attr(split_result, "inside") <- inner_attr

    class(split_result) <- c("nested_cv", class(split_result))

    # return result
    return(split_result)


  }

where <- tidyselect::vars_select_helpers$where


tof_ttest <-
  function(enough_samples, x, y, paired = FALSE) {
    if (enough_samples) {
      return(t.test(x = x, y = y, paired = paired))
    } else {
      return(list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_))
    }
  }

deframe <- function(x) {
  value <- x[[2L]]
  name <- x[[2L]]
  result <- setNames(value, nm = name)
  return(result)
}

#' Coerce flowFrames or flowSets into tof_tbl's.
#'
#' @param flow_data A flowFrame or flowSet
#'
#' @param sep A string indicating which symbol should be used to separate
#' antigen names and metal names in the columns of the output tof_tbl.
#'
#'
#' @export
as_tof_tbl <- function(flow_data, sep = "|") {
  UseMethod("as_tof_tbl")
}

#' @export
as_tof_tbl.flowSet <- function(flow_data, sep = "|") {
  # check if flowset is empty
  if (length(flow_data) < 1) {
    stop("This flowSet is empty.")
  }
  panel_info <-
    flow_data[[1]] %>%
    tof_find_panel_info()

  flowset_exprs <-
    flow_data %>%
    flowCore::fsApply(FUN = flowCore::exprs) %>%
    tibble::as_tibble()

  col_names <-
    stringr::str_c(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  colnames(flowset_exprs) <- col_names

  result <- new_tof_tibble(x = flowset_exprs, panel = panel_info)

  return(result)
}

#' @export
as_tof_tbl.flowFrame <- function(flow_data, sep = "|") {
  panel_info <-
    flow_data %>%
    tof_find_panel_info()

  col_names <-
    stringr::str_c(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  flowframe_exprs <-
    flow_data %>%
    {
      setNames(
        object = tibble::as_tibble(flowCore::exprs(.)),
        nm = col_names
      )
    }

  result <-
    new_tof_tibble(
      x = flowframe_exprs,
      panel = panel_info
    )

  return(result)

}
