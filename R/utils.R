
# reading data -----------------------------------------------------------------

## tidytof_example_data --------------------------------------------------------

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
      system.file("extdata", dataset_name, package = "tidytof", mustWork = TRUE)
    }
  }


## get_extension ---------------------------------------------------------------

#' Find the extension for a file
#'
#' @param filename A string representing the name of a file in its local directory
#'
#' @return The the file extension of `filename`
#'
#'
get_extension <- function(filename) {
  ex <- strsplit(basename(filename), split="\\.")[[1]]
  return(ex[[length(ex)]])
}

# preprocessing/postprocessing -------------------------------------------------

#' Reverses arcsinh transformation with cofactor `scale_factor` and a
#' shift of `shift_factor`.
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
#' @examples
#' shift_factor <- 0
#' scale_factor <- 1/5
#'
#' input_value <- 20
#' asinh_value <- asinh(shift_factor + input_value * scale_factor)
#'
#' restored_value <- rev_asinh(asinh_value, shift_factor, scale_factor)
#'
#'
rev_asinh <- function(x, shift_factor, scale_factor) {

  new_x <- (sinh(x) - shift_factor) / scale_factor
  return(new_x)

}

# downsampling -----------------------------------------------------------------

tof_spade_downsampling <-
  function(densities, cell_ids, target_num_cells, power = 1) {
    densities_sorted <- sort(densities)

    # a fast estimate of the expected number
    # of cells that will be sampled at each density threshold (starting
    # with the most inclusive threshold)
    cdf <- rev(cumsum(1.0 / rev(densities_sorted)))

    max_number_cells <- cdf[[1]]

    boundary <- target_num_cells / max_number_cells

    if (boundary > densities_sorted[1]) {
      targets <- (target_num_cells - (1:length(densities_sorted))) / cdf
      boundary <- targets[which.min(targets - densities_sorted > 0)]
    }

    result_indices <-
      which((boundary / densities)^power > stats::runif(length(densities)))

    result <-
      cell_ids$..cell_id[result_indices]

    return(result)

  }

# metaclustering ---------------------------------------------------------------

tof_summarize_clusters <-
  function(
    tof_tibble,
    cluster_col,
    metacluster_cols = where(tof_is_numeric),
    group_cols,
    central_tendency_function = stats::median
  ) {

    if (missing(group_cols)) {
      group_cols <- NULL
    }

    result <-
      tof_tibble |>
      tof_extract_central_tendency(
        cluster_col = {{cluster_col}},
        group_cols = {{group_cols}},
        marker_cols = {{metacluster_cols}},
        central_tendency_function = central_tendency_function,
        format = "long"
      ) |>
      tidyr::pivot_wider(
        names_from = "channel",
        values_from = "values"
      )

    return(result)
  }

tof_join_metacluster <-
  function(
    tof_tibble,
    meta_tibble,
    cluster_col,
    metacluster_vector
  ) {
    cluster_colname <- colnames(dplyr::select(tof_tibble, {{cluster_col}}))

    cluster_dictionary <-
      meta_tibble |>
      dplyr::select({{cluster_col}}) |>
      dplyr::mutate(.metacluster = metacluster_vector)

    result <-
      tof_tibble |>
      dplyr::select({{cluster_col}}) |>
      dplyr::left_join(cluster_dictionary, by = cluster_colname) |>
      dplyr::select(-{{cluster_col}})

    return(result)
  }

flowsom_consensus <-
  function(data_matrix, num_clusters) {
    tof_tibble <-
      t(data_matrix) |>
      dplyr::as_tibble()

    result <-
      tof_cluster_flowsom(
      tof_tibble = tof_tibble,
      som_xdim = num_clusters,
      som_ydim = 1,
      som_distance_function = "euclidean",
      perform_metaclustering = FALSE
    ) |>
      dplyr::pull(.data$flowsom_cluster)

    return(result)
  }

# mathematical operations ------------------------------------------------------

#' Find the dot product between two vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return The dot product between x and y.
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
magnitude <- function(x) {
  return(sqrt(dot(x, x)))
}

#' L2 normalize an input vector x to a length of 1
#'
#' @param x a numeric vector
#'
#' @return a vector of length length(x) with a magnitude of 1
#'
l2_normalize <- function(x) {
  return(x / magnitude(x))
}

#' Find the cosine similarity between two vectors
#'
#' @param x a numeric vector
#' @param y a numeric vector
#'
#' @return a scalar value representing the cosine similarity between x and y
#'
cosine_similarity <- function(x, y) {
  result <- dot(x, y) / (magnitude(x) * magnitude(y))
  return(result)
}

# differential discovery analysis ----------------------------------------------

#' @importFrom dplyr as_tibble
#'
tidy_lmer_test_glmm <- function(lmer_test) {
  result <-
    summary(lmer_test)$coefficients |>
    dplyr::as_tibble(rownames = "term") |>
    dplyr::rename(
      p.value = .data$`Pr(>|z|)`,
      estimate = .data$Estimate,
      statistic = .data$`z value`,
      std_error = .data$`Std. Error`
    )

  return(result)
}


#' @importFrom dplyr as_tibble
#'
tidy_lmer_test <- function(lmer_test) {
  result <-
    summary(lmer_test)$coefficients |>
    dplyr::as_tibble(rownames = "term") |>
    dplyr::rename(
      p.value = .data$`Pr(>|t|)`,
      estimate = .data$Estimate,
      statistic = .data$`t value`
    )

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
    diffcyt_method = c("glmm", "edgeR", "voom"),
    include_observation_level_random_effects = FALSE
  ) {
    # initialize formula
    my_formula <- NULL

    # extract sample column name as a string
    sample_colname <-
      tidyselect::eval_select(
        expr = rlang::enquo(sample_col),
        data = tof_tibble
      ) |>
      names()

    # extract cluster column name as a string
    cluster_colname <-
      tidyselect::eval_select(
        expr = rlang::enquo(cluster_col),
        data = tof_tibble
      ) |>
      names()

    # extract fixed effect columns as a character vector
    # will return an empty character vector if the argument is missing
    fixed_effect_colnames <-
      tidyselect::eval_select(
        expr = rlang::enquo(fixed_effect_cols),
        data = tof_tibble
      ) |>
      names()

    # extract random effect columns as a character vector
    # will return an empty character vector if the argument is missing
    random_effect_colnames <-
      tidyselect::eval_select(
        expr = rlang::enquo(random_effect_cols),
        data = tof_tibble
      ) |>
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
      tof_tibble |>
      dplyr::select({{sample_col}}, {{fixed_effect_cols}}, {{random_effect_cols}}) |>
      dplyr::distinct() |>
      dplyr::rename(sample_id = {{sample_col}}) |>
      dplyr::arrange(.data$sample_id)

    # create diffcyt marker_info
    marker_info <-
      dplyr::tibble(marker_name = marker_names) |>
      dplyr::mutate(marker_class = "state")

    # create formula or design matrix depending on which method is being used
    if (diffcyt_method %in% c("glmm", "lmm") & include_observation_level_random_effects) {
      random_effect_colnames <-
        c("sample_id", random_effect_colnames)
    }

    if (diffcyt_method %in% c("glmm", "lmm")) {
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
              make_binary_vector(length = length(contrast_names), indices = .x) |>
              diffcyt::createContrast()
          )
      ) |>
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
      tof_tibble |>
      dplyr::group_by({{sample_col}}) |>
      tidyr::nest() |>
      dplyr::arrange({{sample_col}}) |>
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
      data_diff |>
      SummarizedExperiment::rowData()

    temp[,"cluster_id"] <-
      tof_tibble |>
      dplyr::pull({{cluster_col}}) |>
      as.factor()

    SummarizedExperiment::rowData(data_diff) <- temp

    # fix type-related issues in the exprs component of the SummarizedExperiment
    data_exprs <-
      data_diff |>
      SummarizedExperiment::assays() |>
      `[[`("exprs")

    data_colnames <- colnames(data_exprs)

    data_exprs <-
      data_exprs |>
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
#' @importFrom stats glm
#'
fit_da_model <- function(data, formula, has_random_effects = TRUE) {

  # check to see if lme4 is installed
  rlang::check_installed(pkg = "lme4")

  if (!requireNamespace(package = "lme4")) {
    stop("tof_daa_glmm requires the lme4 package to be installed")
  }

  if(is.null(formula)) {
    total_cells <- NULL
  }

  if (has_random_effects) {
    model_fit <-
      lme4::glmer(formula, data, family = "binomial", weights = total_cells)
  } else {
    model_fit <-
      stats::glm(formula, data, family = "binomial", weights = total_cells)
  }

  return(model_fit)
}


#'
#' @importFrom stats glm
#'
fit_de_model <-
  function(data, formula, has_random_effects = TRUE) {

    # check to see if lmerTest is installed
    rlang::check_installed(pkg = "lmerTest")

    if (!requireNamespace(package = "lmerTest")) {
      stop("tof_dea_lmm requires the lmerTest package to be installed")
    }

    if (has_random_effects) {
      model_fit <-
        lmerTest::lmer(formula, data)
    } else {
      model_fit <-
        stats::glm(formula, data, family = "gaussian")
    }
  }



tof_ttest <-
  function(enough_samples, x, y, paired = FALSE) {
    if (enough_samples) {
      return(t.test(x = x, y = y, paired = paired))
    } else {
      return(list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_))
    }
  }


# feature_extraction -----------------------------------------------------------

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
#'
tof_find_emd <- function(vec_1, vec_2, num_bins = 100) {

  # check to see if emdist package is installed
  rlang::check_installed(
    pkg = "emdist",
    reason = "to compute the Earth-mover's distance."
  )

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
    graphics::hist(vec_1, breaks = bins, plot = FALSE)$density |>
    as.matrix()

  density_2 <-
    graphics::hist(vec_2, breaks = bins, plot = FALSE)$density |>
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
#' across the two vectors' combined range. Defaults to 100.
#'
#' @return A double (of length 1) representing the JSD between the two vectors.
#'
#'
#'
tof_find_jsd <- function(vec_1, vec_2, num_bins = 100) {

  # check to see if philentropy package is installed
  rlang::check_installed(
    pkg = "philentropy",
    reason = "to compute the Jensen-Shannon Divergence."
  )

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
    graphics::hist(vec_1, breaks = bins, plot = FALSE)$counts |>
    as.numeric()

  p1 <- p1 / sum(p1)

  p2 <-
    graphics::hist(vec_2, breaks = bins, plot = FALSE)$counts |>
    as.numeric()

  p2 <- p2 / sum(p2)

  # return final result
  js_dist <-
    purrr::quietly(philentropy::JSD)(rbind(p1, p2))$result |>
    as.numeric()

  return(js_dist)
}


#' @importFrom dplyr pull
pull_unless_null <- function(tib, uq_colname) {
  if (is.null(tib)) {
    return(NULL)
  } else {
    return(dplyr::pull(tib, {{uq_colname}}))
  }
}


# local density estimation -----------------------------------------------------


#' Estimate the local densities for all cells in a CyTOF dataset.
#'
#' This function is a wrapper around {tidytof}'s tof_*_density() function family.
#' It performs local density estimation on CyTOF data using a user-specified
#' method (of 3 choices) and each method's corresponding input parameters.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param distance_cols Unquoted names of the columns in `tof_tibble` to use in
#' calculating cell-to-cell distances during the local density estimation for
#' each cell. Defaults to all numeric columns in `tof_tibble`.
#'
#' @param distance_function A string indicating which distance function to use
#' for calculating cell-to-cell distances during local density estimation. Options
#' include "euclidean" (the default) and "cosine".
#'
#' @param normalize A boolean value indicating if the vector of local density
#' estimates should be normalized to values between 0 and 1. Defaults to TRUE.
#'
#' @param ... Additional arguments to pass to the `tof_*_density()` function family
#' member corresponding to the chosen `method`.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' local density estimates of each cell as a new column in `tof_tibble` (TRUE; the default) or if
#' a single-column tibble including only the local density estimates should be returned (FALSE).
#'
#' @param method  A string indicating which local density estimation method should be used.
#' Valid values include "mean_distance", "sum_distance", and "spade".
#'
#' @return A `tof_tbl` or `tibble` If augment = FALSE, it will have a single column encoding
#' the local density estimates for each cell in `tof_tibble`. If augment = TRUE, it will have
#' ncol(tof_tibble) + 1 columns: each of the (unaltered) columns in `tof_tibble`
#' plus an additional column encoding the local density estimates.
#'
#' @family local density estimation functions
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#'
#' # perform the density estimation
#' tof_estimate_density(tof_tibble = sim_data, method = "spade")
#'
#' # perform the density estimation with a smaller search radius around
#' # each cell
#' tof_estimate_density(
#'     tof_tibble = sim_data,
#'     alpha_multiplier = 2,
#'     method = "spade"
#' )
#'
tof_estimate_density <-
  function(
    tof_tibble,
    distance_cols = where(tof_is_numeric),
    distance_function = c("euclidean", "cosine"),
    normalize = TRUE,
    ...,
    augment = TRUE,
    method = c("mean_distance", "sum_distance", "spade")
  ) {
  # check method argument
    method <-
      rlang::arg_match(method, values = c("mean_distance", "sum_distance", "spade"))

    if (method %in% c("mean_distance", "sum_distance")) {
      densities <-
        tof_knn_density(
          tof_tibble = tof_tibble,
          distance_cols = {{distance_cols}},
          distance_function = distance_function,
          estimation_method = method,
          normalize = normalize,
          ...
        )
    } else if (method == "spade") {
      densities <-
        tof_spade_density(
          tof_tibble = tof_tibble,
          distance_cols = {{distance_cols}},
          distance_function = distance_function,
          normalize = normalize,
          ...
        )
    }

    if (augment) {
      result <-
        dplyr::bind_cols(tof_tibble, densities)
    } else {
      result <- densities
    }

    return(result)
}

#' Estimate cells' local densities as done in Spanning-tree Progression Analysis of Density-normalized Events (SPADE)
#'
#' This function uses the algorithm described in
#' \href{https://pubmed.ncbi.nlm.nih.gov/21964415/}{Qiu et al., (2011)} to estimate
#' the local density of each cell in a `tof_tbl` or `tibble` containing CyTOF data.
#' Briefly, this algorithm involves counting the number of neighboring cells
#' within  a sphere of radius alpha surrounding each cell. Here, we do so using
#' the \code{\link[RANN]{nn2}} function.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param distance_cols Unquoted names of the columns in `tof_tibble` to use in
#' calculating cell-to-cell distances during the local density estimation for
#' each cell. Defaults to all numeric columns in `tof_tibble`.
#'
#' @param distance_function A string indicating which distance function to use
#' for calculating cell-to-cell distances during local density estimation. Options
#' include "euclidean" (the default) and "cosine".
#'
#' @param num_alpha_cells An integer indicating how many cells from `tof_tibble`
#' should be randomly sampled from `tof_tibble` in order to estimate `alpha`, the
#' radius of the sphere constructed around each cell during local density
#' estimation. Alpha is calculated by taking the median nearest-neighbor distance
#' from the `num_alpha_cells` randomly-sampled cells and multiplying it by
#' `alpha_multiplier`. Defaults to 2000.
#'
#' @param alpha_multiplier An numeric value indicating the multiplier that should be used
#' when calculating `alpha`, the radius of the sphere constructed around each
#' cell during local density estimation. Alpha is calculated by taking
#' the median nearest-neighbor distance from the `num_alpha_cells` cells
#' randomly-sampled from `tof_tibble` and multiplying it by
#' `alpha_multiplier`. Defaults to 5.
#'
#' @param max_neighbors An integer indicating the maximum number of neighbors
#' that can be counted within the sphere surrounding any given cell. Implemented
#' to reduce the density estimation procedure's speed and memory requirements.
#' Defaults to 1\% of the number of rows in `tof_tibble`.
#'
#' @param normalize A boolean value indicating if the vector of local density
#' estimates should be normalized to values between 0 and 1. Defaults to TRUE.
#'
#' @param ... Additional optional arguments to pass to \code{\link{tof_find_knn}}.
#'
#' @return A tibble with a single column named ".spade_density" containing the
#' local density estimates for each input cell in `tof_tibble`.
#'
#' @family local density estimation functions
#'
#' @export
#'
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#'
#' # perform the density estimation
#' tof_spade_density(tof_tibble = sim_data)
#'
#' # perform the density estimation using cosine distance
#' tof_spade_density(
#'     tof_tibble = sim_data,
#'     distance_function = "cosine",
#'     alpha_multiplier = 2
#' )
#'
#' # perform the density estimation with a smaller search radius around
#' # each cell
#' tof_spade_density(
#'     tof_tibble = sim_data,
#'     alpha_multiplier = 2
#' )
#'
tof_spade_density <-
  function(
    tof_tibble,
    distance_cols = where(tof_is_numeric),
    distance_function = c("euclidean", "cosine"),
    num_alpha_cells = 2000L,
    alpha_multiplier = 5,
    max_neighbors = round(0.01 * nrow(tof_tibble)),
    normalize = TRUE,
    ...
  ) {
    # estimate alpha -----------------------------------------------------------

    # find median nearest-neighbor distances for num_alpha_cells sampled cells
    alpha_knns <-
      tof_tibble |>
      tof_downsample_constant(num_cells = num_alpha_cells) |>
      dplyr::select({{distance_cols}}) |>
      tof_find_knn(
        k = 1L,
        distance_function = distance_function,
        ...
      )
    alpha_median <- median(alpha_knns$neighbor_distances)
    alpha <- alpha_multiplier * alpha_median

    # estimate local densities -------------------------------------------------

    # compute number of neighbors within a sphere of radius alpha for each
    # input cell in tof_tibble
    spade_neighbors <-
      tof_tibble |>
      dplyr::select({{distance_cols}}) |>
      tof_find_knn(
        k = max_neighbors,
        distance_function = distance_function,
        searchtype = "radius",
        radius = alpha
      )

    # calculate densities equal to the number of neighbors
    # (up to max_neighbors) that each cell has
    num_zeros <-
      (spade_neighbors$neighbor_ids == 0) |>
      rowSums()
    densities <- (max_neighbors) - num_zeros

    if (normalize) {
      # normalize densities
      densities <-
        (densities - min(densities)) /
        ((max(densities) - min(densities)))
    }

    # return result
    result <-
      dplyr::tibble(.spade_density = densities)

    return(result)
  }


#' Estimate cells' local densities using K-nearest-neighbor density estimation
#'
#' This function uses the distances between a cell and each of its K nearest
#' neighbors to estimate local density of each cell in a
#' `tof_tbl` or `tibble` containing CyTOF data.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param distance_cols Unquoted names of the columns in `tof_tibble` to use in
#' calculating cell-to-cell distances during the local density estimation for
#' each cell. Defaults to all numeric columns in `tof_tibble`.
#'
#' @param num_neighbors An integer indicating the number of nearest neighbors
#' to use in estimating the local density of each cell. Defaults to the minimum
#' of 15 and the number of rows in `tof_tibble`.
#'
#' @param distance_function A string indicating which distance function to use
#' for calculating cell-to-cell distances during local density estimation. Options
#' include "euclidean" (the default) and "cosine".
#'
#' @param estimation_method A string indicating how the relative density for each cell should be
#' calculated from the distances between it and each of its k nearest neighbors. Options are
#' "mean_distance" (the default; estimates the relative density for a cell's neighborhood by
#' taking the negative average of the distances to its nearest neighbors) and "sum_distance"
#' (estimates the relative density for a cell's neighborhood by taking the negative sum of the
#' distances to its nearest neighbors).
#'
#' @param normalize A boolean value indicating if the vector of local density
#' estimates should be normalized to values between 0 and 1. Defaults to TRUE.
#'
#' @param ... Additional optional arguments to pass to
#' \code{\link{tof_find_knn}}.
#'
#' @return A tibble with a single column named ".knn_density" containing the
#' local density estimates for each input cell in `tof_tibble`.
#'
#' @family local density estimation functions
#'
#'
tof_knn_density <-
  function(
    tof_tibble,
    distance_cols = where(tof_is_numeric),
    num_neighbors = min(15L, nrow(tof_tibble)),
    distance_function = c("euclidean", "cosine"),
    estimation_method = c("mean_distance", "sum_distance"),
    normalize = TRUE,
    ...
  ) {

    # check method argument
    estimation_method <-
      rlang::arg_match(estimation_method, values = c("mean_distance", "sum_distance"))

    # compute knn information
    knn_result <-
      tof_tibble |>
      dplyr::select({{distance_cols}}) |>
      tof_find_knn(
        k = num_neighbors,
        distance_function = distance_function,
        ...
      )

    # find densities using one of 2 methods
    if (estimation_method == "sum_distance") {
      densities <- -base::rowSums(abs(knn_result$neighbor_distances))
    } else if (estimation_method == "mean_distance") {
      densities <- -base::rowMeans(abs(knn_result$neighbor_distances))
    } else {
      stop("Not a valid estimation_method.")
    }

    if (normalize) {
      # normalize densities
      densities <-
        (densities - min(densities)) /
        ((max(densities) - min(densities)))
    }

    result <-
      dplyr::tibble(.knn_density = densities)

    # a tibble with N (number of cells) rows with the ith
    # row representing the KNN-estimated density of the ith cell.
    return(result)
  }


# miscellaneous ----------------------------------------------------------------

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
#' @examples
#' sim_data <-
#'     dplyr::tibble(
#'         cd45 = rnorm(n = 1000),
#'         cd38 = rnorm(n = 1000),
#'         cd34 = rnorm(n = 1000),
#'         cd19 = rnorm(n = 1000)
#'     )
#'
#' # Find the 10 nearest neighbors of each cell in the dataset
#' tof_find_knn(
#'     .data = sim_data,
#'     k = 10,
#'     distance_function = "euclidean"
#' )
#'
#' # Find the 10 approximate nearest neighbors (see RANN::nn2)
#' tof_find_knn(
#'     .data = sim_data,
#'     k = 10,
#'     distance_function = "euclidean",
#'     eps = 0.3
#' )
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



tof_make_knn_graph <-
  function(
    tof_tibble, # single-cell data
    knn_cols, # unquoted column names of columns to use in the knn calculation
    num_neighbors, # number of knn's
    distance_function = c("euclidean", "cosine"),
    knn_error = 0, # eps value in RANN::nn2
    graph_type = c("weighted", "unweighted") # weighted or unweighted graph
  ) {
    knn_data <-
      tof_tibble |>
      # select only knn_cols
      dplyr::select({{knn_cols}}) |>
      tof_find_knn(
        k = num_neighbors,
        distance_function = distance_function,
        eps = knn_error
      )

    # extract knn_ids and put them into long format
    knn_ids <-
      knn_data |>
      purrr::pluck("neighbor_ids")
    colnames(knn_ids) <- 1:ncol(knn_ids)

    knn_ids <-
      knn_ids |>
      dplyr::as_tibble() |>
      dplyr::mutate(from = seq(from = 1, to = nrow(knn_ids), by = 1)) |>
      tidyr::pivot_longer(
        cols = -"from",
        names_to = "neighbor_index",
        values_to = "to"
      )

    if (graph_type == "weighted") {
      # extract knn distances and put them into long format
      knn_dists <-
        knn_data |>
        purrr::pluck("neighbor_distances")
      colnames(knn_dists) <- 1:ncol(knn_dists)

      knn_dists <-
        knn_dists |>
        dplyr::as_tibble() |>
        dplyr::mutate(from = seq(from = 1, to = nrow(knn_dists), by = 1)) |>
        tidyr::pivot_longer(
          cols = -"from",
          names_to = "neighbor_index",
          values_to = "distance"
        )

      # join knn distances with knn ids for final edge tibble
      edge_tibble <-
        knn_ids |>
        dplyr::left_join(knn_dists, by = (c("from", "neighbor_index")))

      if (distance_function == "euclidean") {
        edge_tibble <-
          edge_tibble |>
          dplyr::mutate(weight = 1 / (1 + .data$distance))
      } else {
        edge_tibble <-
          edge_tibble |>
          dplyr::mutate(weight = 1 - .data$distance)
      }

    } else {
      edge_tibble <-
        knn_ids
    }

    # create the knn_graph as a tidygraph object
    knn_graph <-
      tidygraph::tbl_graph(
        nodes = tof_tibble,
        edges = edge_tibble,
        directed = FALSE
      )

    return(knn_graph)

  }


make_binary_vector <- function(length, indices) {
  result <- rep.int(0, times = length)
  result[indices] <- 1

  return(result)
}

deframe <- function(x) {
  value <- x[[2L]]
  name <- x[[1L]]
  result <- setNames(value, nm = name)
  return(result)
}

tof_rescale <-
  function(.vec) {
    vec_max <- max(.vec)
    vec_min <- min(.vec)
    vec_diff <- vec_max - vec_min

    result <- (.vec - vec_min) / vec_diff
    return(result)
  }


#' Make a heatmap summarizing group marker expression patterns in CyTOF data
#'
#' This function makes a heatmap of group-to-group marker expression patterns
#' in single-cell data. Markers are plotted along the horizontal (x-) axis of
#' the heatmap and groups are plotted along the vertical (y-) axis of the
#' heatmap.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param y_col An unquoted column name indicating which column in `tof_tibble`
#' stores the ids for the group to which each cell belongs.
#'
#' @param marker_cols Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be plotted along the x-axis of the heatmap.
#' Supports tidyselect helpers.
#'
#' @param central_tendency_function A function to use for computing the
#' measure of central tendency that will be aggregated from each cluster in
#' cluster_col. Defaults to the median.
#'
#' @param scale_markerwise A boolean value indicating if the heatmap should
#' rescale the columns of the heatmap such that the maximum value for each
#' marker is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param scale_ywise A boolean value indicating if the heatmap should
#' rescale the rows of the heatmap such that the maximum value for each
#' group is 1 and the minimum value is 0. Defaults to FALSE.
#'
#' @param line_width A numeric value indicating how thick the lines separating
#' the tiles of the heatmap should be. Defaults to 0.25.
#'
#' @param theme A ggplot2 theme to apply to the heatmap.
#' Defaults to \code{\link[ggplot2]{theme_minimal}}
#'
#' @param cluster_markers A boolean value indicating if the heatmap should
#' order its columns (i.e. markers) using hierarchical clustering. Defaults to
#' TRUE.
#'
#' @param cluster_groups A boolean value indicating if the heatmap should
#' order its rows (i.e. groups) using hierarchical clustering. Defaults to
#' TRUE.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr as_tibble
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 scale_fill_viridis_c
#' @importFrom ggplot2 theme_minimal
#'
#' @importFrom purrr pluck
#'
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#' @importFrom tidyr pivot_longer
#'
tof_plot_heatmap <-
  function(
    tof_tibble,
    y_col,
    marker_cols = where(tof_is_numeric),
    central_tendency_function = stats::median,
    scale_markerwise = FALSE,
    scale_ywise = FALSE,
    cluster_markers = TRUE,
    cluster_groups = TRUE,
    line_width = 0.25,
    theme = ggplot2::theme_minimal()
  ) {
    # compute summary statistics for each group ------------------------------
    group_tibble <-
      tof_tibble |>
      dplyr::select(
        {{y_col}},
        {{marker_cols}}
      ) |>
      # compute one summary statistic for each group across all knn_cols
      tof_summarize_clusters(
        cluster_col = {{y_col}},
        metacluster_cols = c({{marker_cols}}),
        central_tendency_function = central_tendency_function
      )

    group_matrix <-
      group_tibble |>
      dplyr::select({{marker_cols}}) |>
      as.matrix()

    marker_names <- colnames(group_matrix)

    if (scale_markerwise) {
      group_matrix <-
        apply(X = group_matrix, MARGIN = 2, FUN = tof_rescale)

      # if there are NaN values, tell the user
      if (any(is.nan(group_matrix))) {
        stop("NaN values resulted from marker-wise scaling.
                Consider setting scale_markerwise to FALSE.")
      }
    }

    if (scale_ywise) {
      group_matrix <-
        apply(X = group_matrix, MARGIN = 1, FUN = tof_rescale) |>
        t()

      # if there are NaN values, tell the user
      if (any(is.nan(group_matrix))) {
        stop("NaN values resulted from group-wise scaling.
                Consider setting scale_* to FALSE.")
      }
    }

    colnames(group_matrix) <- marker_names

    if (cluster_markers) {

    marker_order <-
      group_matrix |>
      t() |>
      stats::dist(method = "euclidean") |>
      stats::hclust() |>
      purrr::pluck("order")
    marker_order <- marker_names[marker_order]
    } else {
      marker_order <- marker_names
    }

    if (cluster_groups) {
    group_order <-
      stats::dist(x = group_matrix, method = "euclidean") |>
      stats::hclust() |>
      purrr::pluck("order")
    group_order <- dplyr::pull(group_tibble, {{y_col}})[group_order]
    } else {
      group_order <- dplyr::pull(group_tibble, {{y_col}})
    }

    group_tibble_long <-
      group_matrix |>
      dplyr::as_tibble() |>
      dplyr::bind_cols(dplyr::select(group_tibble, {{y_col}})) |>
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "marker",
        values_to = "expression"
      ) |>
      dplyr::mutate(
        "{{y_col}}" := factor({{y_col}}, levels = group_order),
        marker = factor(.data$marker, levels = marker_order)
      )

    heatmap <-
      group_tibble_long |>
      ggplot2::ggplot(
        ggplot2::aes(x = .data$marker, y = {{y_col}}, fill = .data$expression)
      ) +
      ggplot2::geom_tile(color = "black", linewidth = line_width) +
      ggplot2::scale_fill_viridis_c()

    # rotate x axis labels
    theme <-
      theme +
      ggplot2::theme(
        axis.text.x =
          ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )

    return(heatmap + theme)

  }


flatten_model_predictions <- function(model_prediction, lambdas) {
  num_lambdas <- length(lambdas)
  prediction_list <- list()
  for (i in 1:num_lambdas) {
    lambda <- lambdas[[i]]
    prediction <- model_prediction[,,i]

    prediction_list[[i]] <- dplyr::as_tibble(prediction)
  }

  result <- dplyr::tibble(
    lambda = lambdas,
    predictions = prediction_list
  )

  return(result)
}

tof_model_plot_survival_curves <- function(tof_model, new_x) {
  cox_model <- tof_model$model
  lambda <- tof_model$penalty

  if(missing(new_x)) {
    new_x <-
      tof_setup_glmnet_xy(
        feature_tibble = tof_get_model_training_data(tof_model),
        recipe = tof_model$recipe,
        outcome_cols = dplyr::any_of(tof_model$outcome_colnames),
        model_type = tof_model$model_type
      )$x
  }

  survfit_result <-
    survival::survfit(
      cox_model,
      s = lambda,
      x = tof_get_model_x(tof_model),
      y = tof_get_model_y(tof_model),
      newx = new_x
    )

  times <-
    dplyr::tibble(
      time = survfit_result$time,
      .timepoint_index = 1:length(survfit_result$time)
    )

  survival_curves <-
    survfit_result$surv |>
    dplyr::as_tibble() |>
    dplyr::mutate(.timepoint_index = 1:nrow(survfit_result$surv)) |>
    tidyr::pivot_longer(
      cols = -".timepoint_index",
      names_to = "row_index",
      values_to = "probability"
    ) |>
    dplyr::left_join(times, by = ".timepoint_index") |>
    dplyr::select(-".timepoint_index") |>
    tidyr::nest(survival_curve = -"row_index")

  return(survival_curves)
}

tof_plot_survival_curves <-
  function(cox_model, lambda, recipe, train_x, train_y, new_x) {

    survfit_result <-
      survival::survfit(
        cox_model,
        s = lambda,
        x = train_x,
        y = train_y,
        newx = new_x
      )

    times <-
      dplyr::tibble(
        time = survfit_result$time,
        .timepoint_index = 1:length(survfit_result$time)
      )

    survival_curves <-
      survfit_result$surv |>
      dplyr::as_tibble() |>
      dplyr::mutate(.timepoint_index = 1:nrow(survfit_result$surv)) |>
      tidyr::pivot_longer(
        cols = -".timepoint_index",
        names_to = "row_index",
        values_to = "probability"
      ) |>
      dplyr::left_join(times, by = ".timepoint_index") |>
      dplyr::select(-".timepoint_index") |>
      tidyr::nest(survival_curve = -"row_index")

    return(survival_curves)
  }

#' Compute a Kaplan-Meier curve from sample-level survival data
#'
#' @param survival_curves A tibble from which the Kaplan-Meier curve will be
#' computed. Each row must represent an observation and must have two
#' columns named "time_to_event" and "event".
#'
#' @return A tibble with 3 columns: time_to_event, survival_probability, and
#' is_censored (whether or not an event was censored at that timepoint).
#'
#' @importFrom dplyr add_row
#' @importFrom dplyr arrange
#' @importFrom dplyr lag
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
tof_compute_km_curve <- function(survival_curves) {

  if (!is.factor(survival_curves$event)) {
    survival_curves$event <- factor(survival_curves$event, levels = c(0, 1))
  }

  censor_level <- levels(survival_curves$event)[1]
  death_level <- levels(survival_curves$event)[2]

  km_curve <-
    survival_curves |>
    dplyr::select("time_to_event", "event") |>
    dplyr::arrange(.data$time_to_event) |>
    dplyr::mutate(
      num_current_deaths = as.character(.data$event) == death_level,
      num_deaths = cumsum(.data$num_current_deaths),
      num_current_censored = as.character(.data$event) == censor_level,
      num_censored = cumsum(.data$num_current_censored)
    ) |>
    dplyr::add_row(
      dplyr::tibble(
        time_to_event = 0,
        num_current_deaths = 0,
        num_deaths = 0,
        num_current_censored = 0,
        num_censored = 0
      ),
      .before = 1L
    ) |>
    dplyr::mutate(
      num_at_risk = dplyr::n() - (.data$num_deaths + .data$num_censored + 1),
      num_at_risk = dplyr::lag(.data$num_at_risk, default = dplyr::n() - 1),
      multiplier = (.data$num_at_risk - .data$num_current_deaths) / .data$num_at_risk,
      survival_probability = 1,
      survival_probability = cumprod(.data$multiplier),
      is_censored = (.data$num_current_censored == 1)
    ) |>
   dplyr::select(
     "time_to_event",
     "survival_probability",
     "is_censored"
   )

  return(km_curve)
}

