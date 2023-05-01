# quality_control.R
# This file contains functions relevant to performing quality control analyses
# for high-dimensional cytometry data.

# single-cell ------------------------------------------------------------------

#' Detect low-expression (i.e. potentially failed) channels in high-dimensional cytometry data
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param channel_cols A vector of unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#' If nothing is specified, the default is to analyze all numeric columns.
#'
#' @param negative_threshold A scalar indicating the threshold below
#' which a measurement should be considered negative. Defaults to the hyperbolic
#' arcsine transformation of 10 counts.
#'
#' @param negative_proportion_flag A scalar between 0 and 1 indicating the proportion of
#' cells in tof_tibble that need to be below `negative_threshold` for a given marker
#' in order for that marker to be flagged. Defaults to 0.95.
#'
#'
#' @return A tibble 3 columns and a number of rows equal to the number of
#' columns in `tof_tibble` chosen by `channel_cols`. The three columns are "channel",
#' a character vector of channel names, "negative_proportion", a numeric vector with values
#' between 0 and 1 indicating how many cells in `tof_tibble` below `negative_threshold` for
#' each channel, and `flagged_channel`, a boolean vector indicating whether or not a channel
#' has been flagged as potentially failed (TRUE means that the channel had a large number of
#' cells below `negative_threshold`).
#'
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr arrange
#' @importFrom dplyr as_tibble
#' @importFrom dplyr everything
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#'
#'
#' @examples
#' # simulate some data
#' sim_data <-
#'    data.frame(
#'    cd4 = rnorm(n = 100, mean = 5, sd = 0.5),
#'    cd8 = rnorm(n = 100, mean = 0, sd = 0.1),
#'    cd33 = rnorm(n = 100, mean = 10, sd = 0.1)
#' )
#'
#' tof_assess_channels(tof_tibble = sim_data)
#'
#' tof_assess_channels(tof_tibble = sim_data, channel_cols = c(cd4, cd8))
#'
#' tof_assess_channels(tof_tibble = sim_data, negative_threshold = 2)
#'
tof_assess_channels <-
  function(
    tof_tibble,
    channel_cols = where(tof_is_numeric),
    negative_threshold = asinh(10 / 5),
    negative_proportion_flag = 0.95
  ) {
    proportions <-
      tof_tibble |>
      dplyr::select({{channel_cols}}) |>
      dplyr::summarize(
        dplyr::across(dplyr::everything(), \(x) mean(x < negative_threshold))
      ) |>
      t()

    colnames(proportions) <- "negative_proportion"

    majority_negative_channels <-
      proportions |>
      dplyr::as_tibble(rownames = "channel") |>
      dplyr::mutate(flagged_channel = .data$negative_proportion > negative_proportion_flag) |>
      dplyr::arrange(-.data$negative_proportion)

    return(majority_negative_channels)

  }


#' Calculate the relative flow rates of different timepoints throughout a flow
#' or mass cytometry run.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param time_col An unquoted column name indicating which column in `tof_tibble`
#' contains the time at which each cell was collected.
#'
#' @param num_timesteps The number of bins into which `time_col` should be split.
#' to define "timesteps" of the data collection process. The number of cells analyzed
#' by the cytometer will be counted in each bin separately and will represent
#' the relative average flow rate for that timestep in data collection.
#'
#' @return A tibble with 3 columns and num_timesteps rows. Each row will represent a single
#' timestep (and an error will be thrown if `num_timesteps` is larger than the number of rows in
#' `tof_tibble`). The three columns are as follows: "timestep", a numeric vector indicating which
#' timestep is represented by a given row; "time_window", a factor showing the interval in `time_col`
#' over which "timestep" is defined; and "num_cells", the number of cells that were collected during
#' each timestep.
#'
#' @export
#'
#' @importFrom dplyr count
#' @importFrom dplyr transmute
#'
#' @examples
#'
#' # simulate some data
#' sim_data <-
#'    data.frame(
#'    cd4 = rnorm(n = 100, mean = 5, sd = 0.5),
#'    cd8 = rnorm(n = 100, mean = 0, sd = 0.1),
#'    cd33 = rnorm(n = 100, mean = 10, sd = 0.1),
#'    time = sample(1:300, size = 100)
#' )
#'
#' tof_calculate_flow_rate(tof_tibble = sim_data, time_col = time, num_timesteps = 20L)
#'
#'
#'
tof_calculate_flow_rate <-
  function(tof_tibble, time_col, num_timesteps = nrow(tof_tibble) / 1000) {

    if (num_timesteps > nrow(tof_tibble)) {

      stop("num_timesteps must be smaller than the number of rows in tof_tibble.")
    }

    tof_tibble <-
      tof_tibble |>
      dplyr::transmute(
        time_window = cut({{time_col}}, breaks = num_timesteps),
        timestep = as.numeric(.data$time_window)
      )

    time_counts <-
      tof_tibble |>
      dplyr::count(.data$timestep, .data$time_window, name = "num_cells")

    return(time_counts)
  }



#' Detect flow rate abnormalities in high-dimensional cytometry data (stored in a
#' single data.frame)
#'
#' This function performs a simplified version of
#' \href{https://academic.oup.com/bioinformatics/article/32/16/2473/2240408}{flowAI's}
#' statistical test to detect time periods with abnormal flow rates over the
#' course of a flow cytometry experiment. Briefly, the relative flow rates for each timestep
#' throughout data acquisition are calculated (see \link{tof_calculate_flow_rate}), and
#' outlier timepoints with particularly high or low flow rates (i.e. those beyond
#' extreme values of the t-distribution across timesteps) are flagged.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param time_col An unquoted column name indicating which column in `tof_tibble`
#' contains the time at which each cell was collected.
#'
#' @param num_timesteps The number of bins into which `time_col` should be split.
#' to define "timesteps" of the data collection process. The number of cells analyzed
#' by the cytometer will be counted in each bin separately and will represent
#' the relative average flow rate for that timestep in data collection.
#'
#' @param alpha_threshold A scalar between 0 and 1 indicating the two-tailed significance level
#' at which to draw outlier thresholds in the t-distribution with `num_timesteps` - 1
#' degrees of freedom. Defaults to 0.01.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' computed flags for each cell (see below) as new columns in `tof_tibble` (TRUE) or if
#' a tibble including only the computed flags should be returned (FALSE, the default).
#'
#' @return A tibble with the same number of rows as `tof_tibble`. If augment = FALSE
#' (the default), it will have 3 columns: "{time_col}" (the same column as `time_col`),
#' "timestep" (the numeric timestep to which each cell was assigned based on its
#' value for `time_col`), and "flagged_window" (a boolean vector indicating if
#' each cell was collecting during a timestep flagged for having a high or low
#' flow rate). If augment = TRUE, these 3 columns will be column-bound to `tof_tibble`
#' to return an augmented version of the input dataset. (Note that in this case, time_col will
#' not be duplicated).
#'
#' @importFrom dplyr any_of
#' @importFrom dplyr everything
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr relocate
#' @importFrom dplyr right_join
#'
#' @importFrom stats mad
#' @importFrom stats qt
#' @importFrom stats sd
#'
#' @examples
#'
#' set.seed(1000L)
#' sim_data <-
#'   data.frame(
#'     cd4 = rnorm(n = 1000, mean = 5, sd = 0.5),
#'     cd8 = rnorm(n = 1000, mean = 0, sd = 0.1),
#'     cd33 = rnorm(n = 1000, mean = 10, sd = 0.1),
#'     time =
#'       c(
#'         sample(1:100, size = 200, replace = TRUE),
#'         sample(100:400, size = 300, replace = TRUE),
#'         sample(1:150, size = 400, replace = TRUE),
#'         sample(1:500, size = 100, replace = TRUE)
#'       )
#'   )
#'
#' sim_data |>
#'   tof_assess_flow_rate(
#'     time_col = time,
#'     num_timesteps = 20,
#'     visualize = TRUE
#'   )
#'
tof_assess_flow_rate_tibble <-
  function(
    tof_tibble,
    time_col,
    num_timesteps = nrow(tof_tibble) / 1000,
    alpha_threshold = 0.01,
    augment = FALSE
  ) {

    if (!augment) {
      tof_tibble <-
        tof_tibble |>
        dplyr::select({{time_col}})
    }

    cut_tibble <-
      tof_tibble |>
      dplyr::mutate(
        time_window = cut({{time_col}}, breaks = num_timesteps),
        timestep = as.numeric(.data$time_window)
      ) |>
      dplyr::select(-"time_window")

    time_counts <-
      tof_tibble |>
      tof_calculate_flow_rate(
        time_col = {{time_col}},
        num_timesteps = num_timesteps
      )

    scores <-
      time_counts |>
      dplyr::mutate(
        r_score = (.data$num_cells - median(.data$num_cells)) / stats::mad(.data$num_cells),
        flagged_window =
          .data$r_score > stats::qt(p = 1 - alpha_threshold / 2, df = nrow(time_counts)) |
          .data$r_score < stats::qt(p = alpha_threshold / 2, df = nrow(time_counts))
      )

    result <-
      scores |>
      dplyr::select("flagged_window", "timestep") |>
      dplyr::right_join(cut_tibble, by = "timestep") |>
      dplyr::relocate(dplyr::any_of(colnames(scores)), .after = dplyr::everything())

    return(result)
  }






#' Detect flow rate abnormalities in high-dimensional cytometry data
#'
#' This function performs a simplified version of
#' \href{https://academic.oup.com/bioinformatics/article/32/16/2473/2240408}{flowAI's}
#' statistical test to detect time periods with abnormal flow rates over the
#' course of a flow cytometry experiment. Briefly, the relative flow rates for each timestep
#' throughout data acquisition are calculated (see \link{tof_calculate_flow_rate}), and
#' outlier timepoints with particularly high or low flow rates (i.e. those beyond
#' extreme values of the t-distribution across timesteps) are flagged.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param time_col An unquoted column name indicating which column in `tof_tibble`
#' contains the time at which each cell was collected.
#'
#' @param group_cols Optional. Unquoted column names indicating which columns
#' should be used to group cells before analysis. Flow rate calculation is then
#' performed independently within each group. Supports tidyselect helpers.
#'
#' @param num_timesteps The number of bins into which `time_col` should be split.
#' to define "timesteps" of the data collection process. The number of cells analyzed
#' by the cytometer will be counted in each bin separately and will represent
#' the relative average flow rate for that timestep in data collection.
#'
#' @param alpha_threshold A scalar between 0 and 1 indicating the two-tailed significance level
#' at which to draw outlier thresholds in the t-distribution with `num_timesteps` - 1
#' degrees of freedom. Defaults to 0.01.
#'
#' @param visualize A boolean value indicating if a plot should be generated to
#' visualize each timestep's relative flow rate (by group). This plot is printed
#' to the console as a side-effect of this function; for non-interactive R
#' sessions, this argument should always be FALSE (the default).
#'
#' @param ... Optional additional arguments to pass to \code{\link[ggplot2]{facet_wrap}}.
#' Ignored if visualize = FALSE.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' computed flags for each cell (see below) as new columns in `tof_tibble` (TRUE) or if
#' a tibble including only the computed flags should be returned (FALSE, the default).
#'
#' @return A tibble with the same number of rows as `tof_tibble`. If augment = FALSE
#' (the default), it will have 3 columns: "{time_col}" (the same column as `time_col`),
#' "timestep" (the numeric timestep to which each cell was assigned based on its
#' value for `time_col`), and "flagged_window" (a boolean vector indicating if
#' each cell was collecting during a timestep flagged for having a high or low
#' flow rate). If augment = TRUE, these 3 columns will be column-bound to `tof_tibble`
#' to return an augmented version of the input dataset. (Note that in this case, time_col will
#' not be duplicated).
#'
#' @export
#'
#' @importFrom dplyr count
#' @importFrom dplyr select
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 vars
#'
#' @importFrom purrr map
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#'
#' @examples
#' set.seed(1000L)
#' sim_data <-
#'   data.frame(
#'     cd4 = rnorm(n = 1000, mean = 5, sd = 0.5),
#'     cd8 = rnorm(n = 1000, mean = 0, sd = 0.1),
#'     cd33 = rnorm(n = 1000, mean = 10, sd = 0.1),
#'     file_name = c(rep("a", times = 500), rep("b", times = 500)),
#'     time =
#'       c(
#'         sample(1:100, size = 200, replace = TRUE),
#'         sample(100:400, size = 300, replace = TRUE),
#'         sample(1:150, size = 400, replace = TRUE),
#'         sample(1:500, size = 100, replace = TRUE)
#'       )
#'   )
#'
#' sim_data |>
#'   tof_assess_flow_rate(
#'     time_col = time,
#'     num_timesteps = 20,
#'     visualize = TRUE
#'   )
#'
#' sim_data |>
#'   tof_assess_flow_rate(
#'     time_col = time,
#'     group_cols = file_name,
#'     num_timesteps = 20,
#'     visualize = TRUE
#'   )
#'
#'
tof_assess_flow_rate <-
  function(
    tof_tibble,
    time_col,
    group_cols,
    num_timesteps = nrow(tof_tibble) / 1000,
    alpha_threshold = 0.01,
    visualize = FALSE,
    ...,
    augment = FALSE
  ) {

    if (!augment) {
      tof_tibble <-
        tof_tibble |>
        dplyr::select({{time_col}}, {{group_cols}})
    }

    result <-
      tof_tibble |>
      tidyr::nest(.by = {{group_cols}})

    assessments <-
      purrr::map(
        .x = result$data,
        .f = tof_assess_flow_rate_tibble,
        time_col = {{time_col}},
        num_timesteps = num_timesteps,
        alpha_threshold = alpha_threshold,
        augment = augment
      )

    result$assessment <- assessments

    result <-
      result |>
      dplyr::select(-"data") |>
      tidyr::unnest(cols = "assessment")

    if (visualize) {
      rate_plot <-
        result |>
        dplyr::count(
          {{group_cols}},
          .data$flagged_window,
          .data$timestep,
          name = "num_cells"
        ) |>
        ggplot2::ggplot(ggplot2::aes(x = .data$timestep, y = .data$num_cells, fill = .data$flagged_window)) +
        ggplot2::geom_point(shape = 21) +
        ggplot2::theme_bw()

      if (missing(group_cols)) {
        print(rate_plot)

      } else {

        print(
          rate_plot + ggplot2::facet_wrap(facets = ggplot2::vars({{group_cols}}), ...)
        )
      }
    }

    return(result)

  }

# cluster ----------------------------------------------------------------------

#' Assess a clustering result by calculating the z-score of each cell's
#' mahalanobis distance to its cluster centroid and flagging outliers.
#'
#' This function evaluates the result of a clustering procedure by comparing
#' the mahalanobis distance between each cell and the centroid of the cluster
#' to which it was assigned among all cells in a given cluster. All cells with
#' a mahalanobis-distance z-score above a user-specified threshold are flagged
#' as potentially anomalous. Note that the z-score is calculated using a modified
#' formula to minimize the effect of outliers (Z = x - median(x) / mad(x)).
#'
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param marker_cols  Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be used in the mahalanobis distance calculation.
#' Defaults to all numeric columns. Supports tidyselection.
#'
#' @param z_threshold A scalar indicating the distance z-score threshold above
#' which a cell should be considered anomalous. Defaults to 3.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' computed flags for each cell (see below) as new columns in `tof_tibble` (TRUE) or if
#' a tibble including only the computed flags should be returned (FALSE, the default).
#'
#' @return If augment = FALSE (the default), a tibble with 3 columns:
#' ".mahalanobis_distance" (the mahalanobis distance from each cell to the centroid of
#' tits assigned cluster), "z_score" (the modified z-score of each cell's mahalanobis distance
#' relative to all other cells in the dataset), and "flagged_cell" (a boolean
#' indicating whether or not each cell was flagged as having a z-score above
#' z_threshold). If augment = TRUE, the same 3 columns will be column-bound to
#' tof_tibble, and the resulting tibble will be returned.
#'
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#'
#' @importFrom purrr map
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#'
#' @examples
#'
#' # simulate data
#' sim_data_inner <-
#'     dplyr::tibble(
#'         cd45 = c(rnorm(n = 600), rnorm(n = 500, mean = -4)),
#'         cd38 = c(rnorm(n = 100, sd = 0.5), rnorm(n = 500, mean = -3), rnorm(n = 500, mean = 8)),
#'         cd34 = c(rnorm(n = 100, sd = 0.2, mean = -10), rnorm(n = 500, mean = 4), rnorm(n = 500, mean = 60)),
#'         cd19 = c(rnorm(n = 100, sd = 0.3, mean = 10), rnorm(n = 1000)),
#'         cluster_id = c(rep("a", 100), rep("b", 500), rep("c", 500)),
#'         dataset = "inner"
#'     )
#'
#' sim_data_outer <-
#'   dplyr::tibble(
#'     cd45 = c(rnorm(n = 10), rnorm(50, mean = 3), rnorm(n = 50, mean = -12)),
#'     cd38 =
#'       c(
#'         rnorm(n = 10, sd = 0.5),
#'         rnorm(n = 50, mean = -10),
#'         rnorm(n = 50, mean = 10)
#'       ),
#'     cd34 =
#'       c(
#'         rnorm(n = 10, sd = 0.2, mean = -15),
#'         rnorm(n = 50, mean = 15),
#'         rnorm(n = 50, mean = 70)
#'       ),
#'     cd19 = c(rnorm(n = 10, sd = 0.3, mean = 19), rnorm(n = 100)),
#'     cluster_id = c(rep("a", 10), rep("b", 50), rep("c", 50)),
#'     dataset = "outer"
#'   )
#'
#' sim_data <- rbind(sim_data_inner, sim_data_outer)
#'
#' # detect anomalous cells (in this case, the "outer" dataset contains small
#' # clusters that get lumped into the larger clusters in the "inner" dataset)
#' z_result <-
#'   sim_data |>
#'   tof_assess_clusters_distance(cluster_col = cluster_id, z_threshold = 2.5)
#'
#'
tof_assess_clusters_distance <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    z_threshold = 3,
    augment = FALSE
  ) {
    result <-
      tof_tibble |>
      tidyr::nest(data = -{{cluster_col}})

    distances <-
      purrr::map(
        .x = result$data,
        .f = \(x) {
          tof_cluster_ddpr(
            tof_tibble = x,
            healthy_tibble = dplyr::mutate(x, placeholder = "distance"),
            healthy_label_col = "placeholder",
            cluster_cols = {{marker_cols}},
            distance_function = "mahalanobis",
            num_cores = 1L,
            return_distances = TRUE
          ) |>
            dplyr::select(-".mahalanobis_cluster")
        }
      )

    result$distances <- distances

    result <-
      result |>
      tidyr::unnest(cols = c("data", "distances")) |>
      dplyr::mutate(
        # z_score =
        #   (.data$.mahalanobis_distance - mean(.data$.mahalanobis_distance)) / sd(.data$.mahalanobis_distance),
        z_score = (.data$.mahalanobis_distance - median(.data$.mahalanobis_distance)) / stats::mad(.data$.mahalanobis_distance),
        flagged_cell = .data$z_score > z_threshold
      ) |>
      dplyr::ungroup()


    if (!augment) {
      result <-
        result |>
        dplyr::select(".mahalanobis_distance", "z_score", "flagged_cell")
    }

    return(result)

  }

tof_softmax <- function(vec) {
  vec <- vec - max(vec)
  result <- exp(vec) / sum(exp(vec))

  return(result)
}

tof_shannon_entropy <- function(vec) {
  vec <- vec + 1e-20
  result <- -sum(vec * log(vec))

  which_are_0 <- which(abs(max(vec) - 1) < 1e-20)

  result[which_are_0] <- 0

  return(result)
}

tof_sum_rescale <- function(vec) {
  result <- vec / sum(vec)
  return(result)
}


#' Assess a clustering result by calculating the shannon entropy of each cell's
#' mahalanobis distance to all cluster centroids and flagging outliers.
#'
#' This function evaluates the result of a clustering procedure by calculating
#' the mahalanobis distance between each cell and the centroids of all clusters
#' in the dataset and finding the shannon entropy of the resulting vector of distances.
#' All cells with an entropy threshold above a user-specified threshold are flagged
#' as potentially anomalous. Entropy is minimized (to 0) when a cell is close to
#' one (or a small number) of clusters, but far from the rest of them. If a cell is
#' close to multiple cluster centroids (i.e. has an ambiguous phenotype),
#' its entropy will be large.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param marker_cols Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be used in the mahalanobis distance calculation.
#' Defaults to all numeric columns. Supports tidyselection.
#'
#' @param entropy_threshold A scalar indicating the entropy threshold above
#' which a cell should be considered anomalous. If unspecified, a threshold will
#' be computed using `entropy_quantile` (see below). (Note: Entropy is often between
#' 0 and 1, but can be larger with many classes/clusters).
#'
#' @param entropy_quantile A scalar between 0 and 1 indicating the entropy quantile
#' above which a cell should be considered anomalous. Defaults to 0.9, which means
#' that cells with an entropy above the 90th percentile will be flagged. Ignored
#' if entropy_threshold is specified directly.
#'
#' @param num_closest_clusters An integer indicating how many of a cell's closest
#' cluster centroids should have their mahalanobis distance included in the entropy
#' calculation. Playing with this argument will allow you to ignore distances to
#' clusters that are far away from each cell (and thus may distort the result, as
#' many distant centroids with large distances can artificially inflate a cells'
#' entropy value; that being said, this is rarely an issue empirically).
#' Defaults to all clusters in tof_tibble.
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' computed flags for each cell (see below) as new columns in `tof_tibble` (TRUE) or if
#' a tibble including only the computed flags should be returned (FALSE, the default).
#'
#' @return If augment = FALSE (the default), a tibble with 2 + NUM_CLUSTERS columns.
#' where NUM_CLUSTERS is the number of unique clusters in cluster_col.
#' Two of the columns will be "entropy" (the entropy value for each cell) and "flagged_cell"
#' (a boolean value indicating if each cell had an entropy value above entropy_threshold).
#' The other NUM_CLUSTERS columns will contain the mahalanobis distances from each cell
#' to each of the clusters in cluster_col (named ".mahalanobis_{cluster_name}").
#' If augment = TRUE, the same 2 + NUM_CLUSTERS columns will be column-bound to
#' tof_tibble, and the resulting tibble will be returned.
#'
#' @export
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr tibble
#'
#' @examples
#'
#' # simulate data
#' sim_data <-
#'   dplyr::tibble(
#'     cd45 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd38 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd34 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd19 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cluster_id = c(rep("a", 1000), rep("b", 1000), rep("c", 1000))
#'   )
#'
#' # imagine a "reference" dataset in which "cluster a" isn't present
#' sim_data_reference <-
#'   sim_data |>
#'   dplyr::filter(cluster_id %in% c("b", "c"))
#'
#' # if we cluster into the reference dataset, we will force all cells in
#' # cluster a into a population where they don't fit very well
#' sim_data <-
#'   sim_data |>
#'   tof_cluster(
#'     healthy_tibble = sim_data_reference,
#'     healthy_label_col = cluster_id,
#'     method = "ddpr"
#'   )
#'
#' # we can evaluate the clustering quality by calculating by the entropy of the
#' # mahalanobis distance vector for each cell to all cluster centroids
#' entropy_result <-
#'   sim_data |>
#'   tof_assess_clusters_entropy(
#'     cluster_col = .mahalanobis_cluster,
#'     marker_cols = starts_with("cd"),
#'     entropy_quantile = 0.8,
#'     augment = TRUE
#'   )
#'
#' # most cells in "cluster a" are flagged, and few cells in the other clusters are
#' flagged_cluster_proportions <-
#'   entropy_result |>
#'   group_by(cluster_id) |>
#'   dplyr::summarize(
#'     prop_flagged = mean(flagged_cell)
#'   )
#'
tof_assess_clusters_entropy <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    entropy_threshold,
    entropy_quantile = 0.9,
    num_closest_clusters,
    augment = FALSE
  ) {
    distance_tibble <-
      tof_cluster_ddpr(
        tof_tibble = tof_tibble,
        healthy_tibble = tof_tibble,
        healthy_label_col = {{cluster_col}},
        cluster_cols = {{marker_cols}},
        return_distances = TRUE
      ) |>
      dplyr::select(-".mahalanobis_cluster")

    if (missing(num_closest_clusters)) {
      num_closest_clusters <- ncol(distance_tibble)
    } else if (num_closest_clusters > ncol(distance_tibble)) {
      num_closest_clusters <- ncol(distance_tibble)
    }

    distances <-
      distance_tibble |>
      as.matrix() |>
      apply(MARGIN = 1, FUN = sort, simplify = TRUE) |>
      t()

    distances <- distances[, 1:num_closest_clusters]
    colnames(distances) <- paste0(".mahalanobis_", 1:num_closest_clusters)

    entropies <-
      apply(distances, MARGIN = 1, FUN = tof_sum_rescale)

    entropies <-
      (1 - entropies) |>
      apply(MARGIN = 2, FUN = tof_shannon_entropy) |>
      as.vector()

    if (missing(entropy_threshold)) {
      entropy_threshold <- quantile(entropies, prob = entropy_quantile)
    }

    result <-
      dplyr::tibble(entropy = entropies) |>
      dplyr::mutate(flagged_cell = .data$entropy > entropy_threshold) |>
      dplyr::bind_cols(distance_tibble)

    if (augment) {
      result <- dplyr::bind_cols(tof_tibble, result)
    }

    return(result)

  }

#' Assess a clustering result by calculating a cell's cluster assignment to that
#' of its K nearest neighbors.
#'
#' This function evaluates the result of a clustering procedure by finding the cell's
#' K nearest neighbors, determining which cluster the majority of them are assigned to,
#' and checking if this matches the cell's own cluster assignment. If the cluster
#' assignment of the majority of a cell's nearest neighbors does not match with the
#' cell's own cluster assignment, the cell is flagged as potentially anomalous.
#'
#' @param tof_tibble A `tof_tbl` or `tibble`.
#'
#' @param cluster_col An unquoted column name indicating which column in `tof_tibble`
#' stores the cluster ids for the cluster to which each cell belongs.
#' Cluster labels can be produced via any method the user chooses - including manual gating,
#' any of the functions in the `tof_cluster_*` function family, or any other method.
#'
#' @param marker_cols Unquoted column names indicating which column in `tof_tibble`
#' should be interpreted as markers to be used in the mahalanobis distance calculation.
#' Defaults to all numeric columns. Supports tidyselection.
#'
#' @param num_neighbors An integer indicating how many neighbors should be found
#' during the nearest neighbor calculation.
#'
#' @param distance_function A string indicating which distance function should
#' be used to perform the k nearest neighbor calculation.
#'  Options are "euclidean" (the default) and "cosine".
#'
#' @param augment A boolean value indicating if the output should column-bind the
#' computed flags for each cell (see below) as new columns in `tof_tibble` (TRUE) or if
#' a tibble including only the computed flags should be returned (FALSE, the default).
#'
#' @return If augment = FALSE (the default), a tibble with 2 columns: ".knn_cluster"
#' (a character vector indicating which cluster received the majority vote of each
#' cell's k nearest neighbors) and "flagged_cell" (a boolean value indicating if
#' the cell's cluster assignment matched the majority vote (TRUE) or not (FALSE)).
#' If augment = TRUE, the same 2 columns will be column-bound to
#' tof_tibble, and the resulting tibble will be returned.
#'
#' @export
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#'
#' @examples
#'sim_data <-
#'   dplyr::tibble(
#'     cd45 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd38 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd34 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cd19 = c(rnorm(n = 1000, sd = 1.5), rnorm(n = 1000, mean = 2), rnorm(n = 1000, mean = -2)),
#'     cluster_id = c(rep("a", 1000), rep("b", 1000), rep("c", 1000))
#'   )
#'
#' knn_result <-
#'   sim_data |>
#'   tof_assess_clusters_knn(
#'     cluster_col = cluster_id,
#'     num_neighbors = 10
#'   )
#'
tof_assess_clusters_knn <-
  function(
    tof_tibble,
    cluster_col,
    marker_cols = where(tof_is_numeric),
    num_neighbors = min(10, nrow(tof_tibble)),
    distance_function = c("euclidean", "cosine"),
    augment = FALSE
  ) {
    result <-
      tof_tibble |>
      tof_upsample_neighbor(
        reference_tibble = tof_tibble,
        reference_cluster_col = {{cluster_col}},
        upsample_cols = {{marker_cols}},
        num_neighbors = num_neighbors,
        distance_function = distance_function
      ) |>
      dplyr::rename(.knn_cluster = ".upsample_cluster") |>
      dplyr::bind_cols(dplyr::select(tof_tibble, {{cluster_col}})) |>
      dplyr::mutate(flagged_cell = .data$.knn_cluster != {{cluster_col}}) |>
      dplyr::select(-{{cluster_col}})

    if (augment) {
      result <- dplyr::bind_cols(tof_tibble, result)
    }

    return(result)
  }


