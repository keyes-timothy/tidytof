###    tof_clean_marker_names.R

# Description: 
### Replaces the metals in the colnames() of the output from tof_read_fcs() using an input file
#
# Inputs: 
#     - tof_tibble = expression matrix in tibble form (m cells by n proteins)
#     - lookup_table_path = path to a .csv or excel file containing a lookup 
#        table mapping metals to marker names 
#     - Others?
#
# Outputs: 
#     - tof_tibble = an [m x (n+1)] tibble in which each row represents a single cell 
#                    (of m total in the dataset) and each column represents a marker 
#                    measurement (of n total in the dataset). In addition, the first column 
#                    of the tibble will represent the filename from which the cell was read.
#
# Dependencies: 
#     - tidyverse library
#     - readxl library

# Author: Timothy Keyes
# Version: 2020-06-11

tof_clean_marker_names <- function(tof_tibble = NULL, lookup_table_path = NULL) { 
  NULL
} 



