###    tof_read_fcs.R

# Description: 
### Reads an .fcs file (or a director of .fcs files) into a tibble. 
#
# Inputs: 
#     - file_path = file path to a single .fcs file to be read into the R session. 
#     - folder_path = path to a folder of multiple .fcs files to read into the R session. 
#     - Others?
#
# Outputs: 
#     - tof_tibble = an [m x (n+1)] tibble in which each row represents a single cell 
#                    (of m total in the dataset) and each column represents a metal 
#                    measurement (of n total in the dataset). In addition, the last column 
#                    of the tibble will represent the filename from which the cell was read.
#
# Dependencies: 
#     - flowCore library
#     - tidyverse library
#     - data.table library

# Author: Timothy Keyes
# Version: 2020-09-17

#optimization: -use file_test from built-in {utils} package to remove redundant arguments
#               (just implemented but needs to be tested)
#              -currently only works if all .fcs files have the same panel

tof_read_fcs <- function(file_path = NULL) {
  
  if (file_test(op = "-f", file_path)) {
  tof_tibble <- 
    file_path %>% 
    read.FCS(transformation = FALSE, truncate_max_range = FALSE)
  
  col_names <- 
    if_else(
      are_na(as.character(tof_tibble@parameters@data$desc)), 
      tof_tibble@parameters@data$name, 
      tof_tibble@parameters@data$desc
    ) %>% 
    as.character()
  
  tof_tibble <- 
    tof_tibble %>%
    {
      setNames(
        object = as_tibble(flowCore::exprs(.)), 
        nm = col_names
      )
    }
  
  } else {
    file_names <- list.files(file_path, full.names = FALSE)
    
    tof_tibble <-
      file_path %>% 
      list.files(full.names = TRUE) %>% 
      map(read.FCS, transformation = FALSE, truncate_max_range = FALSE) 
    
    col_names <- 
      if_else(
        are_na(as.character(tof_tibble[[1]]@parameters@data$desc)), 
        tof_tibble[[1]]@parameters@data$name, 
        tof_tibble[[1]]@parameters@data$desc
      ) %>% 
      as.character()
    
    tof_tibble <- 
      tof_tibble %>% 
      map(
        .f =
          ~ setNames(
            object = as_tibble((flowCore::exprs(.x))), 
            nm = col_names
          )
      ) %>%
      map2(.y = file_names, .f = ~ mutate(.x, file_names = .y)) %>% 
      #use data.table for the row binding for speed
      data.table::rbindlist() %>%
      as_tibble()
    
  } 
  return(tof_tibble)
}

