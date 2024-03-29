% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_extraction.R
\name{tof_extract_emd}
\alias{tof_extract_emd}
\title{Extract aggregated features from CyTOF data using earth-mover's distance (EMD)}
\usage{
tof_extract_emd(
  tof_tibble,
  cluster_col,
  group_cols = NULL,
  marker_cols = where(tof_is_numeric),
  emd_col,
  reference_level,
  format = c("wide", "long"),
  num_bins = 100
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

\item{marker_cols}{Unquoted column names representing which columns in `tof_tibble`
(i.e. which CyTOF protein measurements) should be included in the earth-mover's distance
calculation. Defaults to all numeric (integer or double) columns.
Supports tidyselect helpers.}

\item{emd_col}{An unquoted column name that indicates which
column in `tof_tibble` should be used to group cells into different distributions
to be compared with one another during the EMD calculation. For example,
if you want to compare marker expression distributions across stimulation
conditions, `emd_col` should be the column in `tof_tibble` containing
information about which stimulation condition each cell
was exposed to during data acquisition.

If provided, the feature extraction will be
further broken down into subgroups by stimulation condition (and features from each stimulation
condition will be included as their own features in wide format).}

\item{reference_level}{A string indicating what the value in `emd_col`
corresponds to the "reference" value to which all other values in `emd_col`
should be compared. For example, if `emd_col` represents the stimulation
condition for a cell, reference_level might take the value of "basal" or
"unstimulated" if you want to compare each stimulation to the basal state.}

\item{format}{A string indicating if the data should be returned in "wide" format
(the default; each cluster feature is given its own column) or in "long" format
(each cluster feature is provided as its own row).}

\item{num_bins}{Optional. The number of bins to use in dividing one-dimensional
marker distributions into discrete segments for the EMD calculation. Defaults to 100.}
}
\value{
A tibble.

If format == "wide", the tibble will have 1 row for each combination of
the grouping variables provided in `group_cols` and one column for each grouping variable,
one column for each extracted feature (the EMD between the distribution of a
given marker in a given cluster in the basal condition and the distribution of
that marker in a given cluster in a stimulated condition).
The names of each column containing cluster features is obtained using the following pattern:
"\{stimulation_id\}_\{marker_id\}@\{cluster_id\}_emd".

If format == "long", the tibble will have 1 row for each combination of the grouping variables
in `group_cols`, each cluster id (i.e. level) in `cluster_col`, and each marker in `marker_cols`.
It will have one column for each grouping variable, one column for the cluster ids, one
column for the CyTOF channel names, and one column (`value`) containing the features.
}
\description{
This feature extraction function calculates the earth-mover's distance (EMD) between
the stimulated and unstimulated ("basal") experimental conditions of samples in a
CyTOF experiment. This calculation is performed across a
user-specified selection of CyTOF antigens and can be performed either
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

# extract emd of each cluster in each patient (using the "basal" stim
# condition as a reference) in wide format
tof_extract_emd(
    tof_tibble = sim_data,
    cluster_col = cluster_id,
    group_cols = patient,
    emd_col = stim,
    reference_level = "basal"
)

# extract emd of each cluster (using the "basal" stim
# condition as a reference) in long format
tof_extract_emd(
    tof_tibble = sim_data,
    cluster_col = cluster_id,
    emd_col = stim,
    reference_level = "basal",
    format = "long"
)

}
\seealso{
Other feature extraction functions: 
\code{\link{tof_extract_central_tendency}()},
\code{\link{tof_extract_features}()},
\code{\link{tof_extract_jsd}()},
\code{\link{tof_extract_proportion}()},
\code{\link{tof_extract_threshold}()}
}
\concept{feature extraction functions}
