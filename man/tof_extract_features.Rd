% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_extraction.R
\name{tof_extract_features}
\alias{tof_extract_features}
\title{Extract aggregated, sample-level features from CyTOF data.}
\usage{
tof_extract_features(
  tof_tibble,
  cluster_col,
  group_cols = NULL,
  stimulation_col = NULL,
  lineage_cols,
  signaling_cols,
  central_tendency_function = stats::median,
  signaling_method = c("threshold", "emd", "jsd", "central tendency"),
  basal_level = NULL,
  ...
)
}
\arguments{
\item{tof_tibble}{A `tof_tbl` or a `tibble`.}

\item{cluster_col}{An unquoted column name indicating which column in `tof_tibble`
stores the cluster ids of the cluster to which each cell belongs.
Cluster labels can be produced via any method the user chooses - including manual gating,
any of the functions in the `tof_cluster_*` function family, or any other method.}

\item{group_cols}{Unquoted column names representing which columns in `tof_tibble`
should be used to break the rows of `tof_tibble` into subgroups for the feature
extraction calculation. Defaults to NULL (i.e. performing the extraction without subgroups).}

\item{stimulation_col}{Optional. An unquoted column name that indicates which
column in `tof_tibble` contains information about which stimulation condition each cell
was exposed to during data acquisition. If provided, the feature extraction will be
further broken down into subgroups by stimulation condition (and features from each stimulation
condition will be included as their own features in wide format).}

\item{lineage_cols}{Unquoted column names representing which columns in `tof_tibble`
(i.e. which CyTOF protein measurements) should be considered lineage markers in the
feature extraction calculation. Supports tidyselect helpers.}

\item{signaling_cols}{Unquoted column names representing which columns in `tof_tibble`
(i.e. which CyTOF protein measurements) should be considered signaling markers in the
feature extraction calculation. Supports tidyselect helpers.}

\item{central_tendency_function}{The function that will be used to calculate
the measurement of central tendency for each cluster (to be used
as the dependent variable in the linear model). Defaults to \code{\link[stats]{median}}.}

\item{signaling_method}{A string indicating which feature extraction method to use
for signaling markers (as identified by the `signaling_cols` argument). Options are
"threshold" (the default), "emd", "jsd", and "central tendency".}

\item{basal_level}{A string indicating what the value in `stimulation_col`
corresponds to the basal stimulation condition (i.e. "basal" or "unstimulated").}

\item{...}{Optional additional arguments to be passed to tof_extract_threshold,
\code{\link{tof_extract_emd}}, or \code{\link{tof_extract_jsd}}.}
}
\value{
A tibble.

The output tibble will have 1 row for each combination of
the grouping variables provided in `group_cols` (thus, each row will represent
what is considered a single "sample" based on the grouping provided).
It will have one column for each grouping variable and one column for each
extracted feature ("wide" format).
}
\description{
This function wraps other members of the `tof_extract_*` function family to extract
sample-level features from both lineage (i.e. cell surface antigen) CyTOF channels
assumed to be stable across stimulation conditions and signaling CyTOF channels
assumed to change across stimulation conditions. Features are extracted for
each cluster within each independent sample (as defined with the `group_cols` argument).
}
\details{
Lineage channels are specified using the `lineage_cols` argument, and their
extracted features will be measurements of central tendency (as computed by the
user-supplied `central_tendency_function`).

Signaling channels are specified
using the `signaling_cols` argument, and their extracted features will depend
on the user's chosen `signaling_method`. If `signaling method` == "threshold"
(the default), \code{\link{tof_extract_threshold}} will be used to calculate the proportion of
cells in each cluster with signaling marker expression over `threshold` in each
stimulation condition. If `signaling_method` == "emd" or `signaling_method` == "jsd",
\code{\link{tof_extract_emd}} or \code{\link{tof_extract_jsd}} will be used to calculate the earth-mover's
distance (EMD) or Jensen-Shannon Distance (JSD), respectively, between the basal
condition and each of the stimulated conditions in each cluster for each sample.
Finally, if none of these options are chosen, \code{\link{tof_extract_central_tendency}}
will be used to calculate measurements of central tendency.

In addition, \code{\link{tof_extract_proportion}} will be used to extract
the proportion of cells in each cluster will be computed for each sample.

These calculations can be performed either
overall (across all cells in the dataset) or after breaking down the cells into
subgroups using `group_cols`.
}
\examples{
sim_data <-
    dplyr::tibble(
        cd45 = rnorm(n = 1000),
        cd38 = rnorm(n = 1000),
        cd34 = rnorm(n = 1000),
        cd19 = rnorm(n = 1000),
        cluster_id = sample(letters, size = 1000, replace = TRUE),
        patient = sample(c("kirby", "mario"), size = 1000, replace = TRUE),
        stim = sample(c("basal", "stim"), size = 1000, replace = TRUE)
    )

# extract the following features from each cluster in each
# patient/stimulation:
#    - proportion of each cluster
#    - central tendency (median) of cd45 and cd38 in each cluster
#    - the proportion of cells in each cluster with cd34 expression over
#      the default threshold (asinh(10 / 5))
tof_extract_features(
    tof_tibble = sim_data,
    cluster_col = cluster_id,
    group_cols = patient,
    lineage_cols = c(cd45, cd38),
    signaling_cols = cd34,
    stimulation_col = stim
)

# extract the following features from each cluster in each
# patient/stimulation:
#    - proportion of each cluster
#    - central tendency (mean) of cd45 and cd38 in each cluster
#    - the earth mover's distance between each cluster's cd34 histogram in
#      the "basal" and "stim" conditions
tof_extract_features(
    tof_tibble = sim_data,
    cluster_col = cluster_id,
    group_cols = patient,
    lineage_cols = c(cd45, cd38),
    signaling_cols = cd34,
    central_tendency_function = mean,
    stimulation_col = stim,
    signaling_method = "emd",
    basal_level = "basal"
)

}
\seealso{
Other feature extraction functions: 
\code{\link{tof_extract_central_tendency}()},
\code{\link{tof_extract_emd}()},
\code{\link{tof_extract_jsd}()},
\code{\link{tof_extract_proportion}()},
\code{\link{tof_extract_threshold}()}
}
\concept{feature extraction functions}
