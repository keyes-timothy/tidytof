###    tof_preprocess.R

# Description: 
### performs standard mass cytometry transformations on a tof_tibble.  
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or somehing that can be coerced into either 
#       of these
#     - metadata_vars = variables that contain information about cells that should not
#                       computed over, i.e. file names, patient names, stimulation names, etc.
#                       Supports tidy selection using tidy_select helpers. 
#                       Not currently used. 
#     - channel_vars = A vector of non-quoted variables representing columns that contain
#                      single-cell protein measurements. Anything that works in the 
#                      first argument of dplyr::across will work. See ?across. 
#                      Supports tidy selection using tidy_select helpers. The default is to 
#                      transform all numeric columns. 
#     - undo_noise = logical indicating if you'd like to remove the uniform noise that 
#                     Fluidigm software adds to each protein measurement for aesthetic
#                     and visualization purposes. Default = TRUE. 
#     - transform_fun = function to apply to each protein value for variance stabilization. 
#                      default is arcsinh (with a co-factor of 5). 
#
# Outputs: 
#     - processed_tof_tibble = a tibble that performs the specified preprocessing procedure on each 
#                              cell in the input tof_tibble. 
#
# Dependencies: 
#     - tidyverse library

tof_preprocess <- 
  function(
    tof_tibble = NULL, 
    metadata_vars = NULL, 
    channel_vars = where(is.numeric), 
    undo_noise = TRUE, 
    transform_fun = function(x) asinh(x/5)
  ) {
    #channel_vars <- enexprs(channel_vars)
    if (undo_noise) {
    tof_tibble <- 
      tof_tibble %>% 
      mutate(across({{channel_vars}}, ~ floor(.x) + 1))
    }
    tof_tibble <- 
      tof_tibble %>% 
      mutate(across({{channel_vars}}, transform_fun))
    
    return(tof_tibble)
  }




