# differential_discovery.R
# This file contains functions relevant to performing differential discovery
# analyses (differential abundance analysis and differential expression analysis)
# on tof_tbl objects containing CyTOF data.

# diffcyt ----------------------------------------------------------------------

#' Differential Abundance Analysis (DAA) with diffcyt
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tbl` using one of three
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
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the id of the sample from which each cell was collected. `sample_col`
#' should serve as a unique identifier for each sample collected during data acquisition -
#' all cells with the same value for `sample_col` will be treated as a part of the same
#' observational unit.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' abundance analysis. Generally speaking, fixed effects represent the
#' comparisons of biological interest (often the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Optional. Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' abundance analysis. Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age.
#'
#' Note that without multiple samples at each level of each of the
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
#' \href{https://www.nature.com/articles/s42003-019-0415-5}{the diffcyt paper}. Only the
#' "glmm" method can model observation-level random effects, and all other values will ignore
#' this argument (and throw a warning if it is set to TRUE).
#' Defaults to FALSE.
#'
#' @param min_cells An integer value used to filter clusters out of the differential
#' abundance analysis. Clusters are not included in the differential abundance testing
#' if they do not have at least `min_cells` in at least `min_samples` samples.
#' Defaults to 3.
#'
#' @param alpha A numeric value between 0 and 1 indicating which significance
#' level should be applied to multiple-comparison adjusted p-values during the
#' differential abundance analysis. Defaults to 0.05.
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
#' of its values. For example, a fixed effect column named `fixed_effect` with
#' levels "a", "b", and "c" have two terms in `tested_effect`: "fixed_effectb" and
#' "fixed_effectc" (note that level "a" of fixed_effect is set as the reference
#' level during dummy coding). These values correspond to the terms in the
#' differential abundance model that represent the difference in cluster abundances
#' between samples with fixed_effect = "b" and fixed_effect = "a" and between
#' samples with fixed_effect = "c" and fixed_effect = "a", respectively. In addition,
#' the first row in `tested_effect` will always represent the "omnibus"
#' test, or the test that there were significant differences between \emph{any} levels of
#' \emph{any} fixed effect variable in the model.
#'
#' The second column, `daa_results` is a list of tibbles in which each entry gives
#' the differential abundance results for each tested_effect. Within each entry
#' of `daa_results`, you will find several columns including the following:
#' * `p_val`, the p-value associated with each
#' tested effect in each input cluster
#' * `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function)
#' * Other values associated with the underlying method used to perform the
#' differential abundance analysis (such as the log-fold change of cluster
#' abundance between the levels being compared). For details, see
#' \code{\link[edgeR]{glmFit}}, \code{\link[limma]{voom}}, \code{\link[limma]{topTable}},
#' and \code{\link[diffcyt]{testDA_GLMM}}.
#'
#' @family differential abundance analysis functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#' @importFrom tidyselect eval_select
#' @importFrom tidyr unite
#' @importFrom purrr map
#' @importFrom rlang arg_match
#'
tof_daa_diffcyt <-
  function(
    tof_tibble,
    sample_col,
    cluster_col,
    fixed_effect_cols,
    random_effect_cols,
    method = c("glmm", "edgeR", "voom"),
    include_observation_level_random_effects = FALSE,
    min_cells = 3,
    min_samples = 5,
    alpha = 0.05,
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

      include_observation_level_random_effects <- FALSE
    }

    # a hack-y approach for dealing with the diffcyt software - we can pick 2 random
    # columns corresponding to CyTOF measurements to fill the SummarizedExperiment
    # that diffcyt requires later, but because in DAA these measurements are not used,
    # it doesn't matter that we ignore all the other protein measurements.
    #
    # This will make the implementation faster for DAA as well because fewer values will
    # need to by copied into the SummarizedExperiment data structure.
    marker_colnames <-
      tidyselect::eval_select(expr = where(tof_is_numeric), data = tof_tibble) %>%
      names()
    marker_colnames <- marker_colnames[1:2]

    # remove all columns from `tof_tibble` that aren't relevant
    tof_tibble <-
      tof_tibble %>%
      dplyr::select(
        {{sample_col}},
        any_of(marker_colnames),
        {{cluster_col}},
        {{fixed_effect_cols}},
        {{random_effect_cols}}
      )

    diffcyt_args <-
      prepare_diffcyt_args(
        tof_tibble = tof_tibble,
        sample_col = {{sample_col}},
        cluster_col = {{cluster_col}},
        marker_cols = any_of(marker_colnames),
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
        dplyr::transmute(
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
                  dplyr::as_tibble() %>%
                  dplyr::arrange(p_adj) %>%
                  dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
                  dplyr::rename("{{cluster_col}}" := cluster_id)
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
        suppressWarnings(suppressMessages(
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
                      dplyr::as_tibble() %>%
                      dplyr::arrange(p_adj) %>%
                      dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
                      dplyr::select(
                        "{{cluster_col}}" := cluster_id,
                        p_val,
                        p_adj,
                        significant,
                        tidyselect::everything()
                      )
                  }
                )
            )
        ))
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
                  dplyr::as_tibble() %>%
                  dplyr::arrange(p_adj) %>%
                  dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
                  dplyr::select(
                    "{{cluster_col}}" := cluster_id,
                    p_val,
                    p_adj,
                    significant,
                    tidyselect::everything()
                  )
              }
            )
        )
    }

    # remove the omnibus test information (and unnest the results tibble)
    # if there are only 2 levels to the fixed effects being tested (because
    # this means the omnibus test and the individual effect will be identical)
    if (nrow(result_tibble) == 2) {
      result_tibble <-
        result_tibble %>%
        dplyr::filter(tested_effect != "omnibus") %>%
        tidyr::unnest(cols = daa_results)
    }

    return(result_tibble)

  }


#' Differential Expression Analysis (DEA) with diffcyt
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tbl` using one of two
#' methods implemented in the \href{https://www.bioconductor.org/packages/release/bioc/html/diffcyt.html}{diffcyt}
#' package for differential discovery analysis in high-dimensional cytometry
#' data.
#'
#' The two methods are based on linear mixed models ("lmm") and
#' \href{https://academic.oup.com/nar/article/43/7/e47/2414268}{limma} ("limma").
#' Both the "lmm" and "limma" methods can model both fixed effects and random
#' effects.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the id of the sample from which each cell was collected. `sample_col`
#' should serve as a unique identifier for each sample collected during data acquisition -
#' all cells with the same value for `sample_col` will be treated as a part of the same
#' observational unit.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be tested for differential expression between
#' levels of the `fixed_effect_cols`. Defaults to all numeric (integer or double) columns.
#' Supports tidyselect helpers.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' expression analysis. Generally speaking, fixed effects represent the
#' comparisons of biological interest (often the the variables manipulated during
#' experiments), such as treated vs. non-treated, before-treatment vs. after-treatment,
#' or healthy vs. non-healthy.
#'
#' @param random_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' expression analysis. Generally speaking, random effects represent variables
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
#' @param alpha A numeric value between 0 and 1 indicating which significance
#' level should be applied to multiple-comparison adjusted p-values during the
#' differential abundance analysis. Defaults to 0.05.
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
#' test, or the test that there are significant differences between \emph{any} levels of
#' \emph{any} fixed effect variable in the model.
#'
#' The second column, `dea_results` is a list of tibbles in which each entry gives
#' the differential expression results for each tested_effect. Within each entry
#' of `dea_results`, you will find `p_val`, the p-value associated with each
#' tested effect in each input cluster/marker pair; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential expression analysis (such as the log-fold change of clusters' median
#' marker expression values between the conditions being compared). Each tibble in `dea_results`
#' will also have two columns representing the cluster and marker corresponding to the
#' p-value in each row.
#'
#' @family differential expression analysis functions
#'
#' @export
#'
#' @importFrom purrr pluck
#' @importFrom tidyr unite
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
    alpha = 0.05,
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
        {{cluster_col}},
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
        suppressWarnings(suppressMessages(
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
                      dplyr::as_tibble() %>%
                      dplyr::rename(
                        marker = marker_id,
                        "{{cluster_col}}" := cluster_id
                      ) %>%
                      dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
                      dplyr::select(
                        {{cluster_col}},
                        marker,
                        p_val,
                        p_adj,
                        significant,
                        tidyselect::everything()
                      )
                  }
                )
            )
        ))
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
        suppressWarnings(suppressMessages(
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
                      select(-ID) %>%
                      dplyr::rename(
                        marker = marker_id,
                        "{{cluster_col}}" := cluster_id
                      ) %>%
                      dplyr::mutate(significant = dplyr::if_else(p_adj < alpha, "*", "")) %>%
                      dplyr::select(
                        {{cluster_col}},
                        marker,
                        p_val,
                        p_adj,
                        significant,
                        tidyselect::everything()
                      )
                  }
                )
            )
        ))
    }

    # remove the omnibus test information (and unnest the results tibble)
    # if there are only 2 levels to the fixed effects being tested (because
    # this means the omnibus test and the individual effect will be identical)
    if (nrow(result_tibble) == 2) {
      result_tibble <-
        result_tibble %>%
        dplyr::filter(tested_effect != "omnibus") %>%
        tidyr::unnest(cols = dea_results)
    }

    return(result_tibble)

  }



# GLMs and GLMMs ---------------------------------------------------------------

#' Differential Abundance Analysis (DAA) with generalized linear mixed-models (GLMMs)
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tbl` using generalized linear mixed-models. Users
#' specify which columns represent sample, cluster, fixed effect, and random effect
#' information, and a (mixed) binomial regression model is fit using either
#' \code{\link[lme4]{glmer}} or \code{\link[stats]{glm}}.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the id of the sample from which each cell was collected. `sample_col`
#' should serve as a unique identifier for each sample collected during data acquisition -
#' all cells with the same value for `sample_col` will be treated as a part of the same
#' observational unit.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param fixed_effect_cols Unquoted column names representing which columns in
#' `tof_tibble` should be used to model fixed effects during the differential
#' abundance analysis. Supports tidyselect helpers.
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
#' The first column, `tested_effect`,
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
#' The second column, `daa_results`, is a list of tibbles in which each entry gives
#' the differential abundance results for each tested_effect. Within each entry
#' of `daa_results`, you will find `p_value`, the p-value associated with each
#' tested effect in each input cluster; `p_adj`, the multiple-comparison
#' adjusted p-value (using the \code{\link[stats]{p.adjust}} function), and
#' other values associated with the underlying method used to perform the
#' differential abundance analysis (such as the log-fold change of cluster
#' abundance between the levels being compared).
#'
#' @family differential abundance analysis functions
#'
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom stringr str_c
#' @importFrom stats as.formula
#' @importFrom purrr map
#' @importFrom broomExtra tidy
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

    if (length(fixed_effect_colnames) == 0) {
      stop("Fixed effects must be specified. Did you forget to set the `fixed_effect_cols` argument?")
    }

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
      dplyr::pull({{cluster_col}})

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
      suppressMessages(suppressWarnings(
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
          )
      )) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = results) %>%
      dplyr::filter(term != "(Intercept)", !is.na(p.value)) %>%
      dplyr::mutate(p_adj = stats::p.adjust(p.value, method = "fdr")) %>%
      dplyr::arrange(p_adj) %>%
      dplyr::mutate(
        significant = dplyr::if_else(p_adj < alpha, "*", ""),
        mean_fc = exp(estimate)
      ) %>%
      dplyr::rename(
        tested_effect = term,
        p_val = p.value,
        f_statistic = statistic
      )

    if (has_random_effects) {
      fit_data <-
        fit_data %>%
        dplyr::select(-group, -effect)
    }

    fit_data <-
      fit_data %>%
      dplyr::rename_with(stringr::str_replace_all, pattern = "(?<=.)\\.", replacement = "_") %>%
      dplyr::select({{cluster_col}}, p_val, p_adj, significant, tidyselect::everything()) %>%
      tidyr::nest(daa_results = c(-tested_effect))

    # if result tibble only has 1 row (only 2 levels of fixed_effect_cols), \
    # unnest it (which is more intuitive)
    if (nrow(fit_data) == 1) {
      fit_data <-
        tidyr::unnest(fit_data, cols = daa_results)
    }

    return(fit_data)
  }


#' Differential Expression Analysis (DEA) with linear mixed-models (LMMs)
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tbl` using linear mixed-models. Users
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
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param sample_col An unquoted column name indicating which column in `tof_tibble`
#' represents the id of the sample from which each cell was collected. `sample_col`
#' should serve as a unique identifier for each sample collected during data acquisition -
#' all cells with the same value for `sample_col` will be treated as a part of the same
#' observational unit.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
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
#' @param random_effect_cols Optional. Unquoted column names representing which columns in
#' `tof_tibble` should be used to model random effects during the differential
#' expression analysis. Supports tidyselection.
#'
#' Generally speaking, random effects should represent variables
#' that a researcher wants to control/account for, but that are not necessarily
#' of biological interest. Example random effect variables might include batch id,
#' patient id (in a paired design), or patient age. Most analyses will not include random effects.
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
#' marker expression values between the levels being compared).
#'
#' @export
#'
#' @importFrom rlang enquo
#' @importFrom tidyselect eval_select
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom stringr str_c
#' @importFrom stats as.formula
#' @importFrom purrr map
#' @importFrom broomExtra tidy
#'
#' @family differential expression analysis functions
#'
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
      dplyr::pull({{cluster_col}})

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
      suppressWarnings(suppressMessages(
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
          dplyr::mutate(
            significant = dplyr::if_else(p_adj < alpha, "*", ""),
            mean_diff = estimate
          ) %>%
          dplyr::rename(
            tested_effect = term,
            p_val = p.value,
            f_statistic = statistic
          ) %>%
          dplyr::rename_with(stringr::str_replace_all, pattern = "(?<=.)\\.", replacement = "_") %>%
          dplyr::select(
            {{cluster_col}},
            marker,
            p_val,
            p_adj,
            significant,
            tidyselect::everything()
          )
      ))

    if (has_random_effects) {
      fit_data <-
        fit_data %>%
        dplyr::select(-group, -effect)
    }

    fit_data <-
      fit_data %>%
      tidyr::nest(dea_results = c(-tested_effect))

    # if result tibble only has 1 row (only 2 levels of fixed_effect_cols), \
    # unnest it (which is more intuitive)
    if (nrow(fit_data) == 1) {
      fit_data <-
        tidyr::unnest(fit_data, cols = dea_results)
    }

    return(fit_data)
  }




# t-tests ----------------------------------------------------------------------

#' Differential Abundance Analysis (DAA) with t-tests
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tbl` using simple t-tests. Users
#' specify which columns represent sample, cluster, and effect
#' information, and either a paired or unpaired t-test (one per cluster) is used to detect significant
#' differences between sample types.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param effect_col Unquoted column name representing which column in
#' `tof_tibble` should be used to break samples into groups for the t-test. Should
#' only have 2 unique values.
#'
#' @param group_cols Unquoted names of the columns other than `effect_col`
#' that should be used to group cells into independent observations. Fills a similar role
#' to `sample_col` in other `tof_daa_*` functions. For example, if an experiment involves
#' analyzing samples taken from multiple patients at two timepoints (with `effect_col = timepoint`),
#' then group_cols should be the name of the column representing patient IDs.
#'
#' @param test_type A string indicating whether the t-test should be "unpaired"
#' (the default) or "paired".
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
#' @param quiet A boolean value indicating whether warnings should be printed.
#' Defaults to `TRUE`.
#'
#' @return A tibble with 7 columns:
#'
#' \describe{
#'   \item{\strong{\{cluster_col\}}}{The name/ID of the cluster being tested.
#'   Each entry in this column will match a unique value in the input
#'   \{cluster_col\}.}
#'   \item{\strong{t}}{The t-statistic computed for each cluster.}
#'   \item{\strong{df}}{The degrees of freedom used for the t-test for each cluster.}
#'   \item{\strong{p_val}}{The (unadjusted) p-value for the t-test for each cluster.}
#'   \item{\strong{p_adj}}{The \code{\link[stats]{p.adjust}}-adjusted p-value for the t-test for each cluster.}
#'   \item{\strong{significant}}{A character vector that will be "*" for clusters
#'    for which p_adj < alpha and "" otherwise.}
#'   \item{\strong{mean_diff}}{For an unpaired t-test, the difference between the average
#'   proportions of each cluster in the two levels of `effect_col`. For a paired
#'   t-test, the average difference between the proportions of each cluster in
#'   the two levels of `effect_col` within a given patient.}
#'   \item{\strong{mean_fc}}{For an unpaired t-test, the ratio between the average
#'   proportions of each cluster in the two levels of `effect_col`. For a paired
#'   t-test, the average ratio between the proportions of each cluster in
#'   the two levels of `effect_col` within a given patient. 0.001 is added to
#'   the denominator of the ratio to avoid divide-by-zero errors.}
#' }
#'
#' The "levels" attribute of the result indicates the order in which the different
#' levels of the `effect_col` were considered. The `mean_diff` value for each
#' row of the output is computed by subtracting the second level from the first level,
#' and the `mean_fc` value for each row is computed by dividing the first level
#' by the second level.
#'
#' @family differential abundance analysis functions
#'
#' @export
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr nest
#' @importFrom purrr map_if
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom purrr map2_dbl
#' @importFrom stats p.adjust
#'
tof_daa_ttest <-
  function(
    tof_tibble,
    cluster_col,
    effect_col,
    group_cols,
    test_type = c("unpaired", "paired"),
    min_cells = 3,
    min_samples = 5,
    alpha = 0.05,
    quiet = FALSE
  ) {

    # check the test_type argument
    test_type <- rlang::arg_match(test_type, values = c("unpaired", "paired"))

    # check for the number of unique levels in effect_col
    if (missing(group_cols)) {
      stop("The `group_cols` argument must be specified.")
    }

    # check for the number of unique levels in effect_col
    if (missing(effect_col)) {
      stop("The `effect_col` argument must be specified.")
    }

    effect_levels <-
      tof_tibble %>%
      dplyr::pull({{effect_col}}) %>%
      unique()

    if (length(effect_levels) != 2L) {
      stop("`effect_col` must have 2 distinct levels`.")
    }

    # count the cells in each cluster within each sample and effect_col level
    count_df <-
      tof_tibble %>%
      dplyr::mutate("{{cluster_col}}" := as.factor({{cluster_col}})) %>%
      dplyr::count(across({{group_cols}}), {{effect_col}}, {{cluster_col}}, .drop = FALSE) %>%
      dplyr::group_by(across({{group_cols}}), {{effect_col}}) %>%
      dplyr::mutate(prop = n / sum(n)) %>%
      dplyr::ungroup()

    # find the clusters that should be removed due to not having enough cells in
    # enough samples
    clusters_to_keep <-
      count_df %>%
      dplyr::filter(n > min_cells) %>%
      dplyr::count({{cluster_col}}) %>%
      dplyr::filter(n > min_samples) %>%
      dplyr::pull({{cluster_col}})

    count_df <-
      count_df %>%
      dplyr::filter({{cluster_col}} %in% clusters_to_keep)

    # perform the t-tests
    if (test_type == "paired") {
      t_df <-
        count_df %>%
        dplyr::select({{group_cols}}, {{effect_col}}, {{cluster_col}}, prop) %>%
        tidyr::pivot_wider(
          names_from = {{effect_col}},
          values_from = prop
        ) %>%
        tidyr::nest(data = -{{cluster_col}}) %>%
        dplyr::mutate(
          t_test =
            purrr::map_if(
              .x = data,
              .p = ~ nrow(.x) > 1,
              .f = ~ stats::t.test(.x[[effect_levels[[1]]]], .x[[effect_levels[[2]]]], paired = TRUE),
              .else = ~return(list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_))
            ),
          t = purrr::map_dbl(.x = t_test, .f = ~ .x$statistic),
          df = purrr::map_dbl(.x = t_test, .f = ~.x$parameter),
          p_val = purrr::map_dbl(.x = t_test, .f = ~ .x$p.value),
          p_adj = stats::p.adjust(p_val, "fdr"),
          significant = dplyr::if_else(p_adj < alpha, "*", ""),
          mean_diff =
            purrr::map_dbl(.x = data, .f = ~mean(.x[[effect_levels[[1]]]] - .x[[effect_levels[[2]]]])),
          mean_fc =
            # note the correction factor to avoid divide-by-zero errors in calculating the fold-change
            purrr::map_dbl(.x = data, .f = ~mean(.x[[effect_levels[[1]]]] / (.x[[effect_levels[[2]]]] + 0.001)))
        ) %>%
        dplyr::select(-t_test, -data)

    } else {
      # extract cluster_col names as a string
      cluster_colname <-
        count_df %>%
        dplyr::select({{cluster_col}}) %>%
        colnames()

      effect_1_tibble <-
        count_df %>%
        dplyr::filter({{effect_col}} == effect_levels[[1]]) %>%
        tidyr::nest(data = -{{cluster_col}}) %>%
        dplyr::transmute(
          {{cluster_col}},
          effect_1_vector = purrr::map(.x = data, .f = ~ dplyr::pull(.x, prop))
        )

      effect_2_tibble <-
        count_df %>%
        dplyr::filter({{effect_col}} == effect_levels[[2]]) %>%
        tidyr::nest(data = -{{cluster_col}}) %>%
        dplyr::transmute(
          {{cluster_col}},
          effect_2_vector = purrr::map(.x = data, .f = ~ dplyr::pull(.x, prop))
        )

      t_df <-
        dplyr::left_join(effect_1_tibble, effect_2_tibble, by = cluster_colname) %>%
        dplyr::mutate(
          enough_samples =
            purrr::map2_lgl(
              .x = effect_1_vector,
              .y = effect_2_vector,
              .f = ~ length(.x) > 1 & length(.y) > 1
            ),
          t_test =
            purrr::pmap(
              .l = list(enough_samples, effect_1_vector, effect_2_vector),
              .f = tof_ttest
            ),
          t = purrr::map_dbl(.x = t_test, .f = ~ .x$statistic),
          df = purrr::map_dbl(.x = t_test, .f = ~.x$parameter),
          p_val = purrr::map_dbl(.x = t_test, .f = ~ .x$p.value),
          p_adj = stats::p.adjust(p_val, "fdr"),
          significant = dplyr::if_else(p_adj < alpha, "*", ""),
          mean_diff =
            purrr::map2_dbl(
              .x = effect_1_vector,
              .y = effect_2_vector,
              .f = ~ mean(.x, na.rm = TRUE) - mean(.y, na.rm = TRUE),
            ),
          mean_fc =
            purrr::map2_dbl(
              .x = effect_1_vector,
              .y = effect_2_vector,
              .f = ~ mean(.x, na.rm = TRUE) / mean(.y, na.rm = TRUE),
            )
        ) %>%
        dplyr::select(-t_test, -effect_1_vector, -effect_2_vector, -enough_samples) %>%
        dplyr::arrange(p_adj)
    }

    # convert cluster column back to a character vector for consistency with
    # other tidytof functions and arrange columns into standard order
    t_df <-
      dplyr::mutate(t_df, "{{cluster_col}}" := as.character({{cluster_col}})) %>%
      dplyr::select({{cluster_col}}, p_val, p_adj, significant, tidyselect::everything())

    if (!quiet) {
      if (any(is.na(t_df$t))) {
        warning("Some conditions did not have at least 2 replicates. Some NA values returned as a result.\nBe sure to check the `group_cols` argument.")
      }
    }

    # save the level in which the attributes were analyzed in case the user
    # wants to interpret the mean difference or fold change values.
    attr(t_df, "levels") <- effect_levels

    return(t_df)
}


#' Differential Expression Analysis (DEA) with t-tests
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tbl` using simple t-tests. Specifically, either an unpaired
#' or paired t-test will compare samples' marker expression distributions (between two conditions)
#' within each cluster using a user-specified summary function (i.e. mean or median).
#' One t-test is conducted per cluster/marker pair and significant differences
#' between sample types are detected after multiple-hypothesis correction.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids of the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param effect_col Unquoted column name representing which column in
#' `tof_tibble` should be used to break samples into groups for the t-test. Should
#' only have 2 unique values.
#'
#' @param marker_cols Unquoted column names representing which columns in `tof_tibble`
#' (i.e. which CyTOF protein measurements) should be tested for differential expression between
#' levels of the `effect_col`. Defaults to all numeric (integer or double) columns.
#' Supports tidyselect helpers.
#'
#' @param group_cols Unquoted names of the columns other than `effect_col`
#' that should be used to group cells into independent observations. Fills a similar role
#' to `sample_col` in other `tof_daa_*` functions. For example, if an experiment involves
#' analyzing samples taken from multiple patients at two timepoints (with `effect_col = timepoint`),
#' then group_cols should be the name of the column representing patient IDs.
#'
#' @param test_type A string indicating whether the t-test should be "unpaired"
#' (the default) or "paired".
#'
#' @param summary_function The vector-valued function that should be used to
#' summarize the distribution of each marker in each cluster (within each sample, as
#' grouped by `group_cols`). Defaults to `mean`.
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
#' @param quiet A boolean value indicating whether warnings should be printed.
#' Defaults to `TRUE`.
#'
#' @return A tibble with 7 columns:
#'
#' \describe{
#'   \item{\strong{\{cluster_col\}}}{The name/ID of the cluster in the cluster/marker pair
#'   being tested. Each entry in this column will match a unique value in the input
#'   \{cluster_col\}.}
#'   \item{\strong{marker}}{The name of the marker in the cluster/marker pair being tested.}
#'   \item{\strong{t}}{The t-statistic computed for each cluster.}
#'   \item{\strong{df}}{The degrees of freedom used for the t-test for each cluster.}
#'   \item{\strong{p_val}}{The (unadjusted) p-value for the t-test for each cluster.}
#'   \item{\strong{p_adj}}{The \code{\link[stats]{p.adjust}}-adjusted p-value for the t-test for each cluster.}
#'   \item{\strong{significant}}{A character vector that will be "*" for clusters
#'    for which p_adj < alpha and "" otherwise.}
#'   \item{\strong{mean_diff}}{For an unpaired t-test, the difference between the average
#'   proportions of each cluster in the two levels of `effect_col`. For a paired
#'   t-test, the average difference between the proportions of each cluster in
#'   the two levels of `effect_col` within a given patient.}
#'   \item{\strong{mean_fc}}{For an unpaired t-test, the ratio between the average
#'   proportions of each cluster in the two levels of `effect_col`. For a paired
#'   t-test, the average ratio between the proportions of each cluster in
#'   the two levels of `effect_col` within a given patient. 0.001 is added to
#'   the denominator of the ratio to avoid divide-by-zero errors.}
#' }
#'
#' The "levels" attribute of the result indicates the order in which the different
#' levels of the `effect_col` were considered. The `mean_diff` value for each
#' row of the output is computed subtracting the second level from the first level,
#' and the `mean_fc` value for each row is computed by dividing the first level
#' by the second level.
#'
#' @family differential expression analysis functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr nest
#' @importFrom purrr map_if
#' @importFrom purrr map_dbl
#' @importFrom purrr map
#' @importFrom purrr map2_dbl
#' @importFrom purrr map2_lgl
#' @importFrom purrr pmap
#' @importFrom stats p.adjust
#' @importFrom stats t.test
#'
#'
tof_dea_ttest <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols,
    effect_col,
    group_cols,
    test_type = c("unpaired", "paired"),
    summary_function = mean,
    min_cells = 3,
    min_samples = 5,
    alpha = 0.05,
    quiet = FALSE
  ) {
    # check the test_type argument
    test_type <- rlang::arg_match(test_type, values = c("unpaired", "paired"))

    # check for the number of unique levels in effect_col
    if (missing(effect_col)) {
      stop("The `effect_col` argument must be specified.")
    }

    # check for the number of unique levels in effect_col
    if (missing(group_cols)) {
      stop("The `group_cols` argument must be specified.")
    }

    effect_levels <-
      tof_tibble %>%
      dplyr::pull({{effect_col}}) %>%
      unique()

    if (length(effect_levels) != 2L) {
      stop("`effect_col` must have 2 distinct levels`.")
    }

    # summarize marker expression in each cluster within each sample and effect_col level
    expression_df <-
      tof_tibble %>%
      dplyr::mutate("{{cluster_col}}" := as.factor({{cluster_col}})) %>%
      dplyr::group_by(across({{group_cols}}), {{effect_col}}, {{cluster_col}}) %>%
      dplyr::summarize(across({{marker_cols}}, summary_function)) %>%
      dplyr::ungroup()

    # find the clusters that should be removed due to not having enough cells in
    # enough samples
    count_df <-
      tof_tibble %>%
      dplyr::mutate("{{cluster_col}}" := as.factor({{cluster_col}})) %>%
      dplyr::count(across({{group_cols}}), {{effect_col}}, {{cluster_col}}, .drop = FALSE) %>%
      dplyr::group_by(across({{group_cols}}), {{effect_col}}) %>%
      dplyr::ungroup()

    clusters_to_keep <-
      count_df %>%
      dplyr::filter(n > min_cells) %>%
      dplyr::count({{cluster_col}}) %>%
      dplyr::filter(n > min_samples) %>%
      dplyr::pull({{cluster_col}})

    expression_df <-
      expression_df %>%
      dplyr::filter({{cluster_col}} %in% clusters_to_keep)

    # throw an error if you didn't keep any clusters
    if (nrow(expression_df) == 0) {
      stop("No clusters were retained. Try changing `min_cells` and `min_samples`.")
    }

    # perform the t-tests
    if (test_type == "paired") {
      t_df <-
        expression_df %>%
        tidyr::pivot_longer(cols = {{marker_cols}}, names_to = "marker", values_to = "summary") %>%
        tidyr::nest(data = c(-marker, -{{cluster_col}})) %>%
        dplyr::mutate(
          data =
            purrr::map(
              data,
              ~tidyr::pivot_wider(.x, names_from = {{effect_col}}, values_from = summary)
            )
        ) %>%
        dplyr::mutate(
          t_test =
            purrr::map_if(
              .x = data,
              .p = ~ nrow(.x) > 1,
              .f = ~
                stats::t.test(
                  x = .x[[effect_levels[[1]]]],
                  y = .x[[effect_levels[[2]]]],
                  paired = TRUE
                ),
              .else = ~return(list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_))
            ),
          t = purrr::map_dbl(.x = t_test, .f = ~ .x$statistic),
          df = purrr::map_dbl(.x = t_test, .f = ~.x$parameter),
          p_val = purrr::map_dbl(.x = t_test, .f = ~ .x$p.value),
          p_adj = stats::p.adjust(p_val, "fdr"),
          significant = dplyr::if_else(p_adj < alpha, "*", ""),
          mean_diff =
            purrr::map_dbl(
              .x = data,
              .f = ~ mean(.x[[effect_levels[[1]]]] - .x[[effect_levels[[2]]]], na.rm = TRUE)
            ),
          mean_fc =
            purrr::map_dbl(
              .x = data,
              .f = ~ mean(.x[[effect_levels[[1]]]] / (.x[[effect_levels[[2]]]] + 0.001), na.rm = TRUE)
            )
        ) %>%
        dplyr::select(-t_test, -data) %>%
        dplyr::select({{cluster_col}}, marker, p_val, p_adj, significant, tidyselect::everything()) %>%
        dplyr::arrange(p_adj)

    } else {
      t_df <-
        expression_df %>%
        tidyr::pivot_longer(cols = {{marker_cols}}, names_to = "marker", values_to = "summary") %>%
        tidyr::nest(data = c(-{{effect_col}}, -{{cluster_col}}, -marker)) %>%
        tidyr::pivot_wider(names_from = {{effect_col}}, values_from = data) %>%
        dplyr::mutate(
          enough_samples =
            purrr::map2_lgl(
              .x = .data[[effect_levels[[1]]]],
              .y = .data[[effect_levels[[2]]]],
              .f = ~ if (!is.null(.x) & !is.null(.y)) {return(nrow(.x) > 1 & nrow(.y) > 1)} else {return(FALSE)}
            ),
          t_test =
            purrr::pmap(
              .l = list(enough_samples, .data[[effect_levels[[1]]]], .data[[effect_levels[[2]]]]),
              .f = function(x, y, z) tof_ttest(x, y$summary, z$summary)
            ),
          t = purrr::map_dbl(.x = t_test, .f = ~ .x$statistic),
          df = purrr::map_dbl(.x = t_test, .f = ~.x$parameter),
          p_val = purrr::map_dbl(.x = t_test, .f = ~ .x$p.value),
          p_adj = stats::p.adjust(p_val, "fdr"),
          significant = dplyr::if_else(p_adj < alpha, "*", ""),
          mean_diff =
            purrr::map2_dbl(
              .x = .data[[effect_levels[[1]]]],
              .y = .data[[effect_levels[[2]]]],
              .f = ~ if (!is.null(.x) & !is.null(.y)) {return(mean(.x$summary) - mean(.y$summary))} else {return(NA_real_)}
            ),
          mean_fc =
            # note the correction factor to avoid divide-by-zero errors in calculating the fold-change
            purrr::map2_dbl(
              .x = .data[[effect_levels[[1]]]],
              .y = .data[[effect_levels[[2]]]],
              .f = ~ if (!is.null(.x) & !is.null(.y)) {return(mean(.x$summary) / (mean(.y$summary) + 0.001))} else {return(NA_real_)}
            )
        ) %>%
        dplyr::select(-t_test, -any_of(effect_levels), -enough_samples) %>%
        dplyr::select({{cluster_col}}, marker, p_val, p_adj, significant, tidyselect::everything()) %>%
        dplyr::arrange(p_adj)
    }

    # convert cluster column back to a character vector for consistency with
    # other tidytof functions
    t_df <- dplyr::mutate(t_df, "{{cluster_col}}" := as.character({{cluster_col}}))

    if (!quiet) {
      if (any(is.na(t_df$t))) {
        warning("Some conditions did not have at least 2 replicates. Some NA values returned as a result.\nBe sure to check the `group_cols` argument.")
      }
    }

    # save the level in which the attributes were analyzed in case the user
    # wants to interpret the mean difference or fold change values.
    attr(t_df, "levels") <- effect_levels

    return(t_df)
  }

# wrapper functions ------------------------------------------------------------

#' Perform Differential Abundance Analysis (DAA) on CyTOF data
#'
#' This function performs differential abundance analysis on the cell clusters
#' contained within a `tof_tbl` using one of three methods
#' ("diffcyt", "glmm", and "ttest"). It wraps the members of the `tof_daa_*`
#' function family: \code{\link{tof_daa_diffcyt}},
#' \code{\link{tof_daa_glmm}}, and \code{\link{tof_daa_ttest}}.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param daa_method A string indicating which statistical method should be used. Valid
#' values include "diffcyt", "glmm", and "ttest".
#'
#' @param ... Additional arguments to pass onto the `tof_daa_*`
#' function family member corresponding to the chosen method.
#'
#' @return A tibble or nested tibble containing the differential abundance results
#' from the chosen method. See \code{\link{tof_daa_diffcyt}},
#' \code{\link{tof_daa_glmm}}, and \code{\link{tof_daa_ttest}} for details.
#'
#' @family differential abundance analysis functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#'
tof_daa <-
  function(
    tof_tibble,
    daa_method = c("diffcyt", "glmm", "ttest"),
    ...
  ) {
    # check method argument
    method <-
      rlang::arg_match(arg = daa_method, values = c("diffcyt", "glmm", "ttest"))

    if (daa_method == "diffcyt") {
      result <-
        tof_tibble %>%
        tof_daa_diffcyt(...)
    } else if (daa_method == "glmm") {
      result <-
        tof_tibble %>%
        tof_daa_glmm(...)
    } else {
      result <-
        tof_tibble %>%
        tof_daa_ttest(...)
    }
    # return result
    return(result)
  }



#' Perform Differential Expression Analysis (DEA) on CyTOF data
#'
#' This function performs differential expression analysis on the cell clusters
#' contained within a `tof_tbl` using one of three methods
#' ("diffcyt", "glmm", and "ttest"). It wraps the members of the `tof_dea_*`
#' function family: \code{\link{tof_dea_diffcyt}},
#' \code{\link{tof_dea_lmm}}, and \code{\link{tof_dea_ttest}}.
#'
#' @param tof_tibble A `tof_tbl` or a `tibble`.
#'
#' @param dea_method A string indicating which statistical method should be used. Valid
#' values include "diffcyt", "lmm", and "ttest".
#'
#' @param ... Additional arguments to pass onto the `tof_dea_*`
#' function family member corresponding to the chosen method.
#'
#' @return A tibble or nested tibble containing the differential abundance results
#' from the chosen method. See \code{\link{tof_dea_diffcyt}},
#' \code{\link{tof_dea_lmm}}, and \code{\link{tof_dea_ttest}} for details.
#'
#' @family differential expression analysis functions
#'
#' @export
#'
#' @importFrom rlang arg_match
#'
tof_dea <-
  function(
    tof_tibble,
    dea_method = c("diffcyt", "glmm", "ttest"),
    ...
  ) {
    # check method argument
    method <-
      rlang::arg_match(arg = dea_method, values = c("diffcyt", "glmm", "ttest"))

    if (dea_method == "diffcyt") {
      result <-
        tof_tibble %>%
        tof_dea_diffcyt(...)
    } else if (dea_method == "lmm") {
      result <-
        tof_tibble %>%
        tof_dea_lmm(...)
    } else {
      result <-
        tof_tibble %>%
        tof_dea_ttest(...)
    }
    # return result
    return(result)
  }





