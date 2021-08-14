# differential_discovery.R
# This file contains functions relevant to performing differential discovery
# analyses (differential abundance analysis and differential expression analysis)
# on tof_tibble objects containing CyTOF data.

# diffcyt -------------------------------

#' Differential Abundance Analysis (DAA) with diffcyt
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tibble` using one of three
#' methods implemented in the \href{https://www.bioconductor.org/packages/release/bioc/html/diffcyt.html}{diffcyt}
#' package for differential discovery analysis in high-dimensional cytometry
#' data.
#'
#' The three methods are based on generalized linear mixed models ("glmm"),
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/}{edgeR} ("edgeR"), and
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29}{voom} ("voom").
#' While both the "glmm" and "voom" methods can model both fixed effects and random
#' effects, the "edgeR" method can only model fixed effects.
#'
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the sample id of the sample from which each cell was collected.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the differential
#' discovery analysis. Defaults to all numeric (integer or double) columns.
#' Supports tidyselection.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' abundance analysis. Generally speaking, fixed effects should represent the
#' comparisons of biological interest (often the the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' abundance analysis. Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age.
#'
#' Note that without many samples at each level of each of the
#' random effect variables, it can be easy to overfit mixed models. For most CyTOF
#' experiments, 2 or fewer (and often 0) random effect variables are appropriate.
#'
#' @param method A string indicating which diffcyt method should be used for the
#' differential abundance analysis. Valid methods include "glmm" (the default),
#' "edgeR", and "voom".
#'
#' @param include_observation_level_random_effects A boolean value indicating
#' if "observation-level random effects" (OLREs) should be included as random effect
#' terms in a "glmm" differential abundance model. For details about what OLREs are, see
#' \href{https://www.nature.com/articles/s42003-019-0415-5}{the diffcyt paper}.
#' Defaults to FALSE.
#'
#' @param min_cells An integer value used to filter clusters out of the differential
#' abundance analysis. Clusters are not included in the differential abundance testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 3.
#'
#' @param min_samples An integer value used to filter clusters out of the differential
#' abundance analysis. Clusters are not included in the differential abundance testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 5.
#'
#' @param ... Optional additional arguments to pass to the under-the-hood diffcyt
#' function being used to perform the differential abundance analysis. See
#' \code{\link[diffcyt]{testDA_GLMM}}, \code{\link[diffcyt]{testDA_edgeR}}, and
#' \code{\link[diffcyt]{testDA_voom}} for details.
#'
#' @return A nested tibble with two columns: `tested_effect` and `daa_results`.
#'
#' The first column, `tested_effect`
#' is a character vector indicating which term in the differential abundance model
#' was used for significance testing. The values in this row are obtained
#' by pasting together the column names for each fixed effect variable and each
#' of its values. For example, a fixed effect column named fixed_effect with
#' levels "a", "b", and "c" have two terms in `tested_effect`: "fixed_effectb" and
#' "fixed_effectc" (note that level "a" of fixed_effect is set as the reference
#' level during dummy coding). These values correspond to the terms in the
#' differential abundance model that represent the difference in cluster abundances
#' between samples with fixed_effect = "b" and fixed_effect = "a" and between
#' samples with fixed_effect = "c" and fixed_effect = "a", respectively. In addition,
#' note that the first row in `tested_effect` will always represent the "omnibus"
#' test, or the test that there were significant differences between any levels of
#' any fixed effect variable in the model.
#'
#' The second column, `daa_results` is a list of tibbles in which each entry gives
#' the differential abundance results for each tested_effect. Within each entry
#' of `daa_results`, you will find `p_val`, the p-value associated with each
#' tested effect in each input cluster; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential abundance analysis (such as the log-fold change of cluster
#' abundance between the levels being compared).
#'
#' @export
#'
#' @examples
#' NULL
#'
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
    method <- rlang::arg_match(method)

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
        Setting include_observation_level_random_effects to FALSE.\n"
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
                  dplyr::as_tibble()
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
                  dplyr::as_tibble()
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
                  dplyr::as_tibble()
              }
            )
        )
    }

    return(result_tibble)

  }


#' Differential Expression Analysis (DEA) with diffcyt
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tibble` using one of two
#' methods implemented in the \href{https://www.bioconductor.org/packages/release/bioc/html/diffcyt.html}{diffcyt}
#' package for differential discovery analysis in high-dimensional cytometry
#' data.
#'
#' The two methods are based on linear mixed models ("glmm") and
#' \href{https://academic.oup.com/nar/article/43/7/e47/2414268}{limma} ("limma").
#' Both the "lmm" and "limma" methods can model both fixed effects and random
#' effects.
#'
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the sample id of the sample from which each cell was collected.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the differential
#' discovery analysis. Defaults to all numeric (integer or double) columns.
#' Supports tidyselection.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' expression analysis. Generally speaking, fixed effects should represent the
#' comparisons of biological interest (often the the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' expression analysis. Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age.
#'
#' Note that without many samples at each level of each of the
#' random effect variables, it can be easy to overfit mixed models. For most CyTOF
#' experiments, 2 or fewer (and often 0) random effect variables are appropriate.
#'
#' @param method A string indicating which diffcyt method should be used for the
#' differential expression analysis. Valid methods include "lmm" (the default)
#' and "limma".
#'
#' @param include_observation_level_random_effects A boolean value indicating
#' if "observation-level random effects" (OLREs) should be included as random effect
#' terms in a "lmm" differential expression model. For details about what OLREs are, see
#' \href{https://www.nature.com/articles/s42003-019-0415-5}{the diffcyt paper}.
#' Defaults to FALSE.
#'
#' @param min_cells An integer value used to filter clusters out of the differential
#' expression analysis. Clusters are not included in the differential expression testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 3.
#'
#' @param min_samples An integer value used to filter clusters out of the differential
#' expression analysis. Clusters are not included in the differential expression testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 5.
#'
#' @param ... Optional additional arguments to pass to the under-the-hood diffcyt
#' function being used to perform the differential expression analysis. See
#' \code{\link[diffcyt]{testDS_LMM}} and \code{\link[diffcyt]{testDS_limma}}
#' for details.
#'
#' @return A nested tibble with two columns: `tested_effect` and `dea_results`.
#'
#' The first column, `tested_effect`
#' is a character vector indicating which term in the differential expression model
#' was used for significance testing. The values in this row are obtained
#' by pasting together the column names for each fixed effect variable and each
#' of its values. For example, a fixed effect column named fixed_effect with
#' levels "a", "b", and "c" have two terms in `tested_effect`: "fixed_effectb" and
#' "fixed_effectc" (note that level "a" of fixed_effect is set as the reference
#' level during dummy coding). These values correspond to the terms in the
#' differential expression model that represent the difference in cluster median
#' expression values of each marker between samples with fixed_effect = "b" and
#' fixed_effect = "a" and between samples with fixed_effect = "c" and
#' fixed_effect = "a", respectively. In addition,
#' note that the first row in `tested_effect` will always represent the "omnibus"
#' test, or the test that there were significant differences between any levels of
#' any fixed effect variable in the model.
#'
#' The second column, `dea_results` is a list of tibbles in which each entry gives
#' the differential expression results for each tested_effect. Within each entry
#' of `daa_results`, you will find `p_val`, the p-value associated with each
#' tested effect in each input cluster/marker pair; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential expression analysis (such as the log-fold change of clusters' median
#' marker expression values between the levels being compared).
#'
#' @export
#'
#' @examples
#' NULL
#'
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
        Setting include_observation_level_random_effects to FALSE.\n"
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
      dplyr::pull(contrast_matrices) %>%
      purrr::pluck(1)

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
                  dplyr::as_tibble()
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
                  dplyr::as_tibble() %>%
                  select(-ID)
              }
            )
        )
    }

    return(result_tibble)

  }



# GLMs and GLMMs -----------------------

#' Differential Abundance Analysis (DAA) with generalized linear mixed-models (GLMMs)
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tibble` using generalized linear mixed-models. Users
#' specify which columns represent sample, cluster, fixed effect, and random effect
#' information, and a (mixed) binomial regression model is fit using either
#' \code{\link[lme4]{glmer}} or \code{\link[stats]{glm}}.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the sample id of the sample from which each cell was collected.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
#'
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' abundance analysis. Supports tidyselection.
#'
#' Generally speaking, fixed effects should represent the
#' comparisons of biological interest (often the the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' abundance analysis. Supports tidyselection.
#'
#' Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age.
#'
#' Note that without many samples at each level of each of the
#' random effect variables, it can be easy to overfit mixed models. For most CyTOF
#' experiments, 2 or fewer (and often 0) random effect variables are appropriate.
#'
#' @param min_cells An integer value used to filter clusters out of the differential
#' abundance analysis. Clusters are not included in the differential abundance testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 3.
#'
#' @param min_samples An integer value used to filter clusters out of the differential
#' abundance analysis. Clusters are not included in the differential abundance testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 5.
#'
#' @param alpha A numeric value between 0 and 1 indicating which significance
#' level should be applied to multiple-comparison adjusted p-values during the
#' differential abundance analysis. Defaults to 0.05.
#'
#' @return A nested tibble with two columns: `tested_effect` and `daa_results`.
#'
#' The first column, `tested_effect`
#' is a character vector indicating which term in the differential abundance model
#' was used for significance testing. The values in this row are obtained
#' by pasting together the column names for each fixed effect variable and each
#' of its values. For example, a fixed effect column named fixed_effect with
#' levels "a", "b", and "c" have two terms in `tested_effect`: "fixed_effectb" and
#' "fixed_effectc" (note that level "a" of fixed_effect is set as the reference
#' level during dummy coding). These values correspond to the terms in the
#' differential abundance model that represent the difference in cluster abundances
#' between samples with fixed_effect = "b" and fixed_effect = "a" and between
#' samples with fixed_effect = "c" and fixed_effect = "a", respectively. In addition,
#' note that the first row in `tested_effect` will always represent the "omnibus"
#' test, or the test that there were significant differences between any levels of
#' any fixed effect variable in the model.
#'
#' The second column, `daa_results` is a list of tibbles in which each entry gives
#' the differential abundance results for each tested_effect. Within each entry
#' of `daa_results`, you will find `p_value`, the p-value associated with each
#' tested effect in each input cluster; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential abundance analysis (such as the log-fold change of cluster
#' abundance between the levels being compared).
#'
#' @export
#'
#' @examples
#' NULL
#'
tof_daa_glmm <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
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




#' Differential Expression Analysis (DEA) with linear mixed-models (LMMs)
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tibble` using linear mixed-models. Users
#' specify which columns represent sample, cluster, marker, fixed effect, and random effect
#' information, and a (mixed) linear regression model is fit using either
#' \code{\link[lmerTest]{lmer}} or \code{\link[stats]{glm}}.
#'
#' Specifically, one linear model is fit for each cluster/marker pair. For each cluster/marker
#' pair, a user-supplied measurement of central tendency (`central_tendency_function`), such
#' as mean or median, is calculated across all cells in the cluster on a sample-by-sample
#' basis. Then, this central tendency value is used as the dependent variable in a
#' linear model with `fixed_effect_cols` as fixed effects predictors and `random_effect_cols`
#' as random effects predictors. Once all models (one per each cluster/marker pair) are fit,
#' p-values for each coefficient in each model are multiple-comparisons adjusted using the
#' \code{\link[stats]{p.adjust}} function.
#'
#' @param tof_tibble A `tof_tibble` or a `tibble` in which each row represents a
#' single cell and each column represents a CyTOF measurement or a piece of metadata
#' (i.e. cluster id, patient id, etc.) about each cell.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the sample id of the sample from which each cell was collected.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of each cell. These cluster columns can be produced via
#' any method the user chooses, such as manual gating, any of the functions in the
#' `tof_cluster_*` function family, or another method.
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be included in the differential
#' discovery analysis. Defaults to all numeric (integer or double) columns.
#' Supports tidyselection.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' expression analysis. Supports tidyselection.
#'
#' Generally speaking, fixed effects should represent the
#' comparisons of biological interest (often the the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' expression analysis. Supports tidyselection.
#'
#' Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age.
#'
#' @param central_tendency_function The function that will be used to calculate
#' the measurement of central tendency for each cluster/marker pair (to be used
#' as the dependent variable in the linear model). Defaults to \code{\link[stats]{median}}.
#'
#' @param min_cells An integer value used to filter clusters out of the differential
#' expression analysis. Clusters are not included in the differential expression testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 3.
#'
#' @param min_samples An integer value used to filter clusters out of the differential
#' expression analysis. Clusters are not included in the differential expression testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 5.
#'
#' @param alpha A numeric value between 0 and 1 indicating which significance
#' level should be applied to multiple-comparison adjusted p-values during the
#' differential abundance analysis. Defaults to 0.05.
#'
#' @return A nested tibble with two columns: `tested_effect` and `dea_results`.
#'
#' The first column, `tested_effect`
#' is a character vector indicating which term in the differential expression model
#' was used for significance testing. The values in this row are obtained
#' by pasting together the column names for each fixed effect variable and each
#' of its values. For example, a fixed effect column named fixed_effect with
#' levels "a", "b", and "c" have two terms in `tested_effect`: "fixed_effectb" and
#' "fixed_effectc" (note that level "a" of fixed_effect is set as the reference
#' level during dummy coding). These values correspond to the terms in the
#' differential expression model that represent the difference in cluster median
#' expression values of each marker between samples with fixed_effect = "b" and
#' fixed_effect = "a" and between samples with fixed_effect = "c" and
#' fixed_effect = "a", respectively. In addition,
#' note that the first row in `tested_effect` will always represent the "omnibus"
#' test, or the test that there were significant differences between any levels of
#' any fixed effect variable in the model.
#'
#' The second column, `dea_results` is a list of tibbles in which each entry gives
#' the differential expression results for each tested_effect. Within each entry
#' of `daa_results`, you will find `p_val`, the p-value associated with each
#' tested effect in each input cluster/marker pair; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential expression analysis (such as the log-fold change of clusters' median
#' marker expression values between the levels being compared).#'
#' @export
#'
#' @examples
#' NULL
#'
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


