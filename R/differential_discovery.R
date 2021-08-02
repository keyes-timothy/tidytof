# differential_discovery.R
# This file contains functions relevant to performing differential discovery
# analyses (differential abundance analysis and differential expression analysis)
# on tof_tibble objects containing CyTOF data.

# diffcyt -------------------------------

tof_daa_diffcyt <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    fixed_effect_cols,
    random_effect_cols,
    method = c("glmm", "edgeR", "voom"),
    include_observation_level_random_effects = FALSE,
    min_cells = 3,
    min_samples = 5,
    ...
  ) {

    # check to see if the diffcyt package is installed
    has_diffcyt <- requireNamespace(package = "diffcyt")
    if (!has_diffcyt) {
      stop(
        "This function requires the {diffcyt} package. Install it with this code:\n
           if (!requireNamespace(\"BiocManager\", quietly = TRUE))
           install.packages(\"BiocManager\")
           BiocManager::install(\"diffcyt\")"
      )
    }

    # check method argument
    method <- match.arg(method, choices = c("glmm", "edgeR", "voom"))

    # edgeR can't model random effects, so we throw an error for the user
    # if they are included
    if (method == "edgeR" & !missing(random_effect_cols)) {
      stop(
        "edgeR can't model random effects. Trying using another method or
           model everything as a fixed effect."
      )
    }

    # Only the "glmm" method supports observation-level random effects, so
    # provide a warning if include_observation_level_random_effects = TRUE
    # for any other method.
    if (include_observation_level_random_effects == TRUE & method != "glmm") {
      warning(
        "Warning: Only the \"glmm\" method can use observation-level random effects.
        Setting include_observation_level_random_effects to FALSE."
      )
    }

    # remove all columns from `tof_tibble` that aren't relevant
    tof_tibble <-
      tof_tibble %>%
      dplyr::select(
        {{sample_col}},
        {{marker_cols}},
        {{fixed_effect_cols}},
        {{random_effect_cols}}
      )

    diffcyt_args <-
      prepare_diffcyt_args(
        tof_tibble = tof_tibble,
        sample_col = {{sample_col}},
        cluster_col = {{cluster_col}},
        marker_cols = {{marker_cols}},
        fixed_effect_cols = {{fixed_effect_cols}},
        random_effect_cols = {{random_effect_cols}},
        method = method,
        include_observation_level_random_effects = include_observation_level_random_effects
      )

    # find counts of each cluster in all samples
    cell_counts <- diffcyt::calcCounts(diffcyt_args$data_diff)

    # perform difference abundance testing
    if (method == "glmm") {
      # if glmms are being used,

      result_tibble <-
        diffcyt_args$contrast_matrix_tibble %>%
        transmute(
          tested_effect = contrast_names,
          daa_results =
            map(
              .x = contrast_matrices,
              .f = function(x) {
                diffcyt::testDA_GLMM(
                  d_counts = cell_counts,
                  formula = diffcyt_args$my_formula,
                  contrast = x,
                  min_cells = min_cells,
                  min_samples = min_samples,
                  ...
                ) %>%
                  diffcyt::topTable(all = TRUE, show_all_cols = TRUE) %>%
                  tibble::as_tibble()
              }
            )
        )


    }

    else if (method == "voom") {
      # if limma/voom is being used,

      # We unite all random effect columns and treat them as
      # a single block ID. Note this occurs in the help file and that
      # it is not recommended to use more than 1 random effect variable with
      # the voom method.

      if (length(diffcyt_args$random_effect_colnames) != 0) {
        # if there are random effects, combine them into a single block_id
        block_id <-
          diffcyt_args$experiment_info %>%
          tidyr::unite(
            col = "block_id",
            tidyselect::any_of(diffcyt_args$random_effect_colnames)
          ) %>%
          dplyr::pull(block_id) %>%
          as.factor()

      } else {
        # otherwise, don't include a block_id
        block_id <- NULL
      }

      result_tibble <-
        diffcyt_args$contrast_matrix_tibble %>%
        dplyr::transmute(
          tested_effect = contrast_names,
          daa_results =
            purrr::map(
              .x = contrast_matrices,
              .f = function(x) {
                diffcyt::testDA_voom(
                  d_counts = cell_counts,
                  design = diffcyt_args$my_design,
                  contrast = x,
                  block_id = block_id,
                  min_cells = min_cells,
                  min_samples = min_samples,
                  ...
                ) %>%
                  diffcyt::topTable(all = TRUE, show_all_cols = TRUE) %>%
                  tibble::as_tibble()
              }
            )
        )

    }

    else {

      result_tibble <-
        diffcyt_args$contrast_matrix_tibble %>%
        dplyr::transmute(
          tested_effect = contrast_names,
          daa_results =
            purrr::map(
              .x = contrast_matrices,
              .f = function(x) {
                diffcyt::testDA_edgeR(
                  d_counts = cell_counts,
                  design = diffcyt_args$my_design,
                  contrast = x,
                  min_cells = min_cells,
                  min_samples = min_samples,
                  ...
                ) %>%
                  diffcyt::topTable(all = TRUE, show_all_cols = TRUE) %>%
                  tibble::as_tibble()
              }
            )
        )
    }

    return(result_tibble)

  }


tof_dea_diffcyt <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    fixed_effect_cols,
    random_effect_cols,
    method = c("lmm", "limma"),
    include_observation_level_random_effects = FALSE,
    min_cells = 3,
    min_samples = 5,
    ...
  ) {

    # check to see if the diffcyt package is installed
    has_diffcyt <- requireNamespace(package = "diffcyt")
    if (!has_diffcyt) {
      stop(
        "This function requires the {diffcyt} package. Install it with this code:\n
           if (!requireNamespace(\"BiocManager\", quietly = TRUE))
           install.packages(\"BiocManager\")
           BiocManager::install(\"diffcyt\")"
      )
    }

    # check method argument
    method <- match.arg(method, choices = c("lmm", "limma"))

    # Only the "lmm" method supports observation-level random effects, so
    # provide a warning if include_observation_level_random_effects = TRUE
    # for any other method.
    if (include_observation_level_random_effects == TRUE & method != "lmm") {
      warning(
        "Warning: Only the \"lmm\" method can use observation-level random effects.
        Setting include_observation_level_random_effects to FALSE."
      )
    }

    # remove all columns from `tof_tibble` that aren't relevant
    tof_tibble <-
      tof_tibble %>%
      dplyr::select(
        {{sample_col}},
        {{marker_cols}},
        {{fixed_effect_cols}},
        {{random_effect_cols}}
      )

    diffcyt_args <-
      prepare_diffcyt_args(
        tof_tibble = tof_tibble,
        sample_col = {{sample_col}},
        cluster_col = {{cluster_col}},
        marker_cols = {{marker_cols}},
        fixed_effect_cols = {{fixed_effect_cols}},
        random_effect_cols = {{random_effect_cols}},
        method = method,
        include_observation_level_random_effects = include_observation_level_random_effects
      )

    # find cluster counts and cluster medians
    cell_counts <- diffcyt::calcCounts(diffcyt_args$data_diff)
    cell_medians <- diffcyt::calcMedians(diffcyt_args$data_diff)

    # Perform the differential expression analysis

    my_contrast <-
      diffcyt_args$contrast_matrix_tibble %>%
      pull(contrast_matrices) %>%
      pluck(1)

    if (method == "lmm") {
      # if lmm's are being used,

      result_tibble <-
        diffcyt_args$contrast_matrix_tibble %>%
        dplyr::transmute(
          tested_effect = contrast_names,
          dea_results =
            purrr::map(
              .x = contrast_matrices,
              .f = function(x) {
                diffcyt::testDS_LMM(
                  d_counts = cell_counts,
                  d_medians = cell_medians,
                  formula = diffcyt_args$my_formula,
                  contrast = x,
                  markers_to_test = rep(TRUE, nrow(diffcyt_args$marker_info)),
                  min_cells = min_cells,
                  min_samples = min_samples,
                  ...
                ) %>%
                  diffcyt::topTable(all = TRUE, show_all_cols = TRUE) %>%
                  tibble::as_tibble()
              }
            )
        )
    }

    else if (method == "limma") {
      # if limma is being used,

      if (length(diffcyt_args$random_effect_colnames) != 0) {
        # if there are random effects, combine them into a single block_id
        block_id <-
          diffcyt_args$experiment_info %>%
          tidyr::unite(
            col = "block_id",
            tidyselect::any_of(diffcyt_args$random_effect_colnames)
          ) %>%
          dplyr::pull(block_id) %>%
          as.factor()

      } else {
        # otherwise, don't include a block_id
        block_id <- NULL
      }

      result_tibble <-
        diffcyt_args$contrast_matrix_tibble %>%
        dplyr::transmute(
          tested_effect = contrast_names,
          dea_results =
            purrr::map(
              .x = contrast_matrices,
              .f = function(x) {
                diffcyt::testDS_limma(
                  d_counts = cell_counts,
                  d_medians = cell_medians,
                  design = diffcyt_args$my_design,
                  contrast = x,
                  block_id = block_id,
                  min_cells = min_cells,
                  min_samples = min_samples,
                  markers_to_test = rep(TRUE, nrow(diffcyt_args$marker_info)),
                  ...
                ) %>%
                  diffcyt::topTable(all = TRUE, show_all_cols = TRUE) %>%
                  tibble::as_tibble() %>%
                  select(-ID)
              }
            )
        )
    }

    return(result_tibble)

  }



# GLMs and GLMMs -----------------------

tof_daa_glmm <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    #marker_cols = where(tof_is_numeric),
    fixed_effect_cols,
    random_effect_cols,
    min_cells = 3,
    min_samples = 5,
    alpha = 0.05
  ) {

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


    # count cells in all samples
    cell_counts <-
      tof_tibble %>%
      dplyr::group_by(
        {{sample_col}},
        {{cluster_col}},
        dplyr::across(tidyselect::any_of(c(fixed_effect_colnames, random_effect_colnames)))
      ) %>%
      dplyr::count(name = "num_cells") %>%
      dplyr::group_by({{sample_col}}) %>%
      dplyr::mutate(
        total_cells = sum(num_cells),
        prop = num_cells / total_cells
      ) %>%
      dplyr::ungroup()

    # find the clusters that don't have over the threshold of minimum cells
    # in over the threshold of minimum samples
    clusters_to_remove <-
      cell_counts %>%
      dplyr::count({{sample_col}}, {{cluster_col}}, wt = num_cells) %>%
      dplyr::mutate(has_over_min_cells = n > min_cells) %>%
      dplyr::count({{cluster_col}}, has_over_min_cells) %>%
      dplyr::filter(has_over_min_cells) %>%
      dplyr::filter(n < min_samples) %>%
      pull({{cluster_col}})

    cell_counts <-
      cell_counts %>%
      dplyr::filter(!({{cluster_col}} %in% clusters_to_remove))

    # nest the count data so we can fit one model per cluster
    fit_data <-
      cell_counts %>%
      dplyr::group_by({{cluster_col}}) %>%
      tidyr::nest() %>%
      dplyr::ungroup()

    # specify if there are random effects
    if (length(random_effect_colnames) == 0) {
      has_random_effects <- FALSE
    } else {
      has_random_effects <- TRUE
    }

    # construct formula for each model
    if (has_random_effects) {
      formula_string <-
        stringr::str_c(
          "prop ~ ",
          stringr::str_c(fixed_effect_colnames, sep = "+", collapse = " + "),
          "+",
          stringr::str_c(paste0("(1 | ", random_effect_colnames, ")"), sep = "+")
        )
    } else {
      formula_string <-
        stringr::str_c(
          "prop ~ ",
          stringr::str_c(fixed_effect_colnames, sep = "+", collapse = " + "),
          sep = ""
        )
    }

    formula <- stats::as.formula(formula_string)

    # fit one model per cluster
    fit_data <-
      fit_data %>%
      dplyr::mutate(
        results =
          purrr::map(
            .x = data,
            .f = fit_da_model,
            formula = formula,
            has_random_effects = has_random_effects
          ),
        results = purrr::map(.x = results, .f = broomExtra::tidy)
      ) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = results) %>%
      dplyr::filter(term != "(Intercept)", !is.na(p.value)) %>%
      dplyr::mutate(p_adj = stats::p.adjust(p.value, method = "fdr")) %>%
      dplyr::arrange(p_adj) %>%
      dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
      dplyr::rename(tested_effect = term) %>%
      dplyr::rename_with(stringr::str_replace_all, pattern = "\\.", replacement = "_") %>%
      tidyr::nest(daa_results = c(-tested_effect))

    return(fit_data)
  }

tof_dea_lmm <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    fixed_effect_cols,
    random_effect_cols,
    central_tendency_function = median,
    min_cells = 3,
    min_samples = 5,
    alpha = 0.05
  ) {

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

    # count cells in all samples
    cell_counts <-
      tof_tibble %>%
      dplyr::count({{sample_col}}, {{cluster_col}}, name = "num_cells")

    # find the clusters that don't have over the threshold of minimum cells
    # in over the threshold of minimum samples
    clusters_to_remove <-
      cell_counts %>%
      dplyr::mutate(has_over_min_cells = num_cells > min_cells) %>%
      dplyr::count({{cluster_col}}, has_over_min_cells) %>%
      dplyr::filter(has_over_min_cells) %>%
      dplyr::filter(n < min_samples) %>%
      pull({{cluster_col}})

    # find the median for each cluster in each sample
    expression_data <-
      tof_tibble %>%
      dplyr::select(
        {{cluster_col}},
        {{sample_col}},
        tidyselect::any_of(c(fixed_effect_colnames, random_effect_colnames)),
        {{marker_cols}}
      ) %>%
      # remove clusters that don't fit the minimum criteria
      dplyr::filter(!({{cluster_col}} %in% clusters_to_remove)) %>%
      dplyr::group_by(
        {{sample_col}},
        {{cluster_col}},
        dplyr::across(tidyselect::any_of(c(fixed_effect_colnames, random_effect_colnames)))
      ) %>%
      dplyr::summarize(
        dplyr::across(
          tidyselect::everything(),
          .fns = central_tendency_function
        )
      ) %>%
      dplyr::ungroup()

    # nest the expression data so we can fit one model per cluster per channel
    fit_data <-
      expression_data %>%
      tidyr::pivot_longer(
        cols = {{marker_cols}},
        names_to = "marker",
        values_to = "expression"
      ) %>%
      dplyr::group_by({{cluster_col}}, marker) %>%
      tidyr::nest() %>%
      dplyr::ungroup()

    # specify if there are random effects
    if (length(random_effect_colnames) == 0) {
      has_random_effects <- FALSE
    } else {
      has_random_effects <- TRUE
    }

    # construct formula for each model
    if (has_random_effects) {
      formula_string <-
        stringr::str_c(
          "expression ~ ",
          stringr::str_c(fixed_effect_colnames, sep = "+", collapse = " + "),
          "+",
          stringr::str_c(paste0("(1 | ", random_effect_colnames, ")"), sep = "+")
        )
    } else {
      formula_string <-
        stringr::str_c(
          "expression ~ ",
          stringr::str_c(fixed_effect_colnames, sep = "+", collapse = " + "),
          sep = ""
        )
    }

    formula <- stats::as.formula(formula_string)

    # fit one model per cluster
    fit_data <-
      fit_data %>%
      dplyr::mutate(
        results =
          purrr::map(
            .x = data,
            .f = fit_de_model,
            formula = formula,
            has_random_effects = has_random_effects
          ),
        results = purrr::map(.x = results, .f = broomExtra::tidy)
      ) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = results) %>%
      dplyr::filter(term != "(Intercept)", !is.na(p.value)) %>%
      dplyr::mutate(p_adj = stats::p.adjust(p.value, method = "fdr")) %>%
      dplyr::arrange(p_adj) %>%
      dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
      dplyr::rename(tested_effect = term) %>%
      dplyr::rename_with(stringr::str_replace_all, pattern = "\\.", replacement = "_") %>%
      tidyr::nest(daa_results = c(-tested_effect))

    return(fit_data)
  }


# cydar -----------------------

tof_daa_cydar <- function(tof_tibble) {
  # TO DO
  return(NULL)
}


