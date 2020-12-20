###    tof_write_out.R

# Description: 
### Reads an .fcs file (or a director of .fcs files) into a tibble. 
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or somehing that can be coerced into either 
#       of these
#     - group_vars = variables that should be used to break cells into individual files
#     - out_folder_path = file path to the folder where the output files should be saved.
#     - format = format for the files being written. Currently supports .csv and .fcs files. 
#     - Others?
#
# Outputs: 
#     - none
#
# Side-effects: 
#       - Save .csv or .fcs files in the designated folder. 
#
# Dependencies: 
#     - flowCore library
#     - tidyverse library
#     - data.table library

# Author: Timothy Keyes
# Version: 2020-06-12

tof_write_out <- 
  function(
    tof_tibble = NULL, 
    group_vars = "file_names", 
    out_folder_path = NULL, 
    format = NULL
  ) { 
    if ("csv" %in% format) { 
      tof_tibble %>% 
        group_by(across(one_of(group_vars))) %>% 
        nest(data = everything()) %>%
        ungroup() %>% 
        mutate(prefix = str_c(setdiff(colnames(.), "data"))) %>% 
        map2(
          .x =
        )
        pmap(
          .l = ., 
          .f = 
            ~ 
            write_csv(
              x = ..1, 
              path = file.path(out_folder_path, str_c(!!! group_vars, ".csv"))
            )
        )
    }
    if ("fcs" %in% format) { 
      tof_tibble %>% 
        group_by({group_vars}) %>% 
        nest(data = everything()) %>%
        map(
          .x = data,
          .f = write.fcs, 
          path = file.path(out_folder_path, str_c({group_vars}, ".csv"))
        )
    }
    
  }
