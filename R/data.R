
#' CyTOF data from two samples: 5,000 B-cell lineage cells from a healthy
#' patient and 5,000 B-cell lineage cells from a B-cell precursor Acute
#' Lymphoblastic Leukemia (BCP-ALL) patient.
#'
#' A dataset containing CyTOF measurements from immune cells originally studied
#' in the following paper:  \cr \cr
#' Good Z, Sarno J, et al.
#' Single-cell developmental classification of B cell precursor acute
#' lymphoblastic leukemia at diagnosis reveals predictors of relapse.
#' Nat Med. 2018 May;24(4):474-483. doi: 10.1038/nm.4505. Epub 2018 Mar 5.
#' PMID: 29505032; PMCID: PMC5953207.
#'
#'
#' @format A data frame with 10000 rows and 24 variables:
#' \describe{
#'   \item{sample_name}{name of the sample from which the data was read}
#'   \item{cd45}{A CyTOF measurement in raw ion counts}
#'   \item{cd19}{A CyTOF measurement in raw ion counts}
#'   \item{cd22}{A CyTOF measurement in raw ion counts}
#'   \item{cd79b}{A CyTOF measurement in raw ion counts}
#'   \item{cd20}{A CyTOF measurement in raw ion counts}
#'   \item{cd34}{A CyTOF measurement in raw ion counts}
#'   \item{cd179a}{A CyTOF measurement in raw ion counts}
#'   \item{cd123}{A CyTOF measurement in raw ion counts}
#'   \item{cd10}{A CyTOF measurement in raw ion counts}
#'   \item{cd179b}{A CyTOF measurement in raw ion counts}
#'   \item{cd24}{A CyTOF measurement in raw ion counts}
#'   \item{cd127}{A CyTOF measurement in raw ion counts}
#'   \item{cd43}{A CyTOF measurement in raw ion counts}
#'   \item{cd38}{A CyTOF measurement in raw ion counts}
#'   \item{cd58}{A CyTOF measurement in raw ion counts}
#'   \item{cd3}{A CyTOF measurement in raw ion counts}
#'   \item{psyk}{A CyTOF measurement in raw ion counts}
#'   \item{p4ebp1}{A CyTOF measurement in raw ion counts}
#'   \item{pstat5}{A CyTOF measurement in raw ion counts}
#'   \item{pakt}{A CyTOF measurement in raw ion counts}
#'   \item{ps6}{A CyTOF measurement in raw ion counts}
#'   \item{perk}{A CyTOF measurement in raw ion counts}
#'   \item{pcreb}{A CyTOF measurement in raw ion counts}
#' }
#' @source \url{https://github.com/kara-davis-lab/DDPR}
"ddpr_data"


#' Clinical metadata for each patient sample in Good & Sarno et al. (2018).
#'
#' A dataset containing patient-level clinical metadata for samples originally studied
#' in the following paper: \cr \cr
#' Good Z, Sarno J, et al.
#' Single-cell developmental classification of B cell precursor acute
#' lymphoblastic leukemia at diagnosis reveals predictors of relapse.
#' Nat Med. 2018 May;24(4):474-483. doi: 10.1038/nm.4505. Epub 2018 Mar 5.
#' PMID: 29505032; PMCID: PMC5953207.
#'
#'
#' @format A data frame with 10000 rows and 12 variables:
#' \describe{
#'   \item{patient_id}{Name of the sample from which the data was read}
#'   \item{gender}{Gender of the patient from which each sample was collected}
#'   \item{age_at_diagnosis}{Age (in years) of the patient from which each sample was collected}
#'   \item{wbc_count}{The diagnostic White Blood Cell (WBC) count of the patient from which each sample was collected}
#'   \item{mrd_risk}{Risk stratification category for each patient using minimal residual disease (MRD) criteria}
#'   \item{nci_rome_risk}{Risk stratification category for each patient using National Cancer Institute (NCI) criteria}
#'   \item{relapse_status}{A string representing whether or not a patient relapsed}
#'   \item{time_to_relapse}{The time (in days) it took each patient to relapse. Patients who did not relapse will have the value of NA}
#'   \item{type_of_relapse}{
#'   A string representing the timing of relapse for each patient.
#'   "Very early" relapses occurred less than 18 months after diagnosis;
#'   "Early" relapses occurred between 18 months and 32 months after diagnosis;
#'   "Late" relapses occurred later than 32 months after diagnosis.
#'   }
#'   \item{ccr}{The number of documented days of continuous complete remission (CCR) for patients who did not relapse. All patients who relapsed will have a value of NA.}
#'   \item{cohort}{A string representing if each sample was used in the "Training" or "Validation" cohort in the original study}
#'   \item{ddpr_risk}{The risk category ("Low" or "High") assigned to each sample using the original paper's risk-stratification algorithm}
#' }
#' @source Good Z, Sarno J, et al.
#' Single-cell developmental classification of B cell precursor acute
#' lymphoblastic leukemia at diagnosis reveals predictors of relapse.
#' Nat Med. 2018 May;24(4):474-483. doi: 10.1038/nm.4505. Epub 2018 Mar 5.
#' PMID: 29505032; PMCID: PMC5953207. Supplementary Table 1.
"ddpr_metadata"

#' CyTOF data from 6,000 healthy immune cells from a single patient.
#'
#' A dataset containing CyTOF measurements from healthy control cells originally studied
#' in the following paper:  \cr \cr
#' Levine JH, Simonds EF, et al.
#' Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that
#' Correlate with Prognosis. Cell. 2015 Jul 2;162(1):184-97.
#' doi: 10.1016/j.cell.2015.05.047. Epub 2015 Jun 18. PMID: 26095251;
#' PMCID: PMC4508757.
#'
#' 2000 cells from 3 clusters identified in the original paper have been
#' sampled.
#'
#' @format A data frame with 6000 rows and 26 variables:
#' \describe{
#'   \item{sample_name}{Name of the sample from which the data was read}
#'   \item{phenograph_cluster}{Numeric ID of the cluster assignment of each row}
#'   \item{cd19}{A CyTOF measurement in raw ion counts}
#'   \item{cd11b}{A CyTOF measurement in raw ion counts}
#'   \item{cd34}{A CyTOF measurement in raw ion counts}
#'   \item{cd45}{A CyTOF measurement in raw ion counts}
#'   \item{cd123}{A CyTOF measurement in raw ion counts}
#'   \item{cd33}{A CyTOF measurement in raw ion counts}
#'   \item{cd47}{A CyTOF measurement in raw ion counts}
#'   \item{cd7}{A CyTOF measurement in raw ion counts}
#'   \item{cd15}{A CyTOF measurement in raw ion counts}
#'   \item{cd44}{A CyTOF measurement in raw ion counts}
#'   \item{cd38}{A CyTOF measurement in raw ion counts}
#'   \item{cd3}{A CyTOF measurement in raw ion counts}
#'   \item{cd117}{A CyTOF measurement in raw ion counts}
#'   \item{cd64}{A CyTOF measurement in raw ion counts}
#'   \item{cd41}{A CyTOF measurement in raw ion counts}
#'   \item{pstat3}{A CyTOF measurement in raw ion counts}
#'   \item{pstat5}{A CyTOF measurement in raw ion counts}
#'   \item{pampk}{A CyTOF measurement in raw ion counts}
#'   \item{p4ebp1}{A CyTOF measurement in raw ion counts}
#'   \item{ps6}{A CyTOF measurement in raw ion counts}
#'   \item{pcreb}{A CyTOF measurement in raw ion counts}
#'   \item{pzap70-syk}{A CyTOF measurement in raw ion counts}
#'   \item{prb}{A CyTOF measurement in raw ion counts}
#'   \item{perk1-2}{A CyTOF measurement in raw ion counts}
#' }
#' @source \url{https://cytobank.org/nolanlab/reports/Levine2015.html}
"phenograph_data"

#' A character vector of metal name patterns supported by tidytof.
#'
#' A character vector used by `tof_read_fcs` and `tof_read_data` to detect and
#' parse which CyTOF metals correspond to each channel in an input .fcs file.
#'
#' @format A character vector in which each entry is a pattern that tidytof searches
#' for in every CyTOF channel in input .fcs files. These patterns are an amalgamate
#' of example .fcs files sampled from the studies linked below.
#'
#' @source \url{https://github.com/kara-davis-lab/DDPR}
#' \url{https://cytobank.org/nolanlab/reports/Levine2015.html}
#' \url{https://cytobank.org/nolanlab/reports/Spitzer2015.html}
#' \url{https://cytobank.org/nolanlab/reports/Spitzer2017.html}
#' \url{https://community.cytobank.org/cytobank/projects/609}
"metal_masterlist"
