###    tof_clean_marker_names.R

# Description: 
### Replaces the metals in the colnames() of the output from tof_read_fcs() using an input file
#
# Inputs: 
#     - tof_tibble = expression matrix in tibble form (m cells by n proteins)
#     - sep = regular expression representing the separator between metal and marker names 
#     - Others?
#
# Outputs: 
#     - tof_tibble = an [m x (n+1)] tibble in which each row represents a single cell 
#                    (of m total in the dataset) and each column represents a marker 
#                    measurement (of n total in the dataset). In addition, the last column 
#                    of the tibble will represent the filename from which the cell was read.
#
# Dependencies: 
#     - tidyverse library

# Author: Timothy Keyes
# Version: 2020-06-12

# Notes: Currently, will only perform if the column names are in the format "metal{sep}protein", 
#        where `sep` cannot also be used in the protein name. 

tof_clean_marker_names <- function(tof_tibble = NULL, sep = "_") { 
  tof_tibble %>% 
    rename_with(
      .fn = ~ str_extract(string = .x, pattern = str_glue("(?<={sep}).+")), #this use of regular expressions should be optimized
      .cols = contains(sep)
    ) %>%
    return()
} 

