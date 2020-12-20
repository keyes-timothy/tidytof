#### geom_tof_ridges.R ####
# Description: A function that encodes a geom specifically for generating univariate histograms
#              of CyTOF data. 
#             Currently VERY EXPERIMENTAL and NOT OPTIMIZED
#
# Inputs: 
#     - rel_min_height, size, and scale = arguments to pass to geom_density_ridges_gradient() 
#                                         (with nice defaults set)
#     - ... = additional arguments to pass to geom_density_ridges_gradient()
#
# Outputs: None
#
# Dependencies: 
#     - ggridges package and all of its dependencies

geom_tof_ridges <- function(rel_min_height = 0.00, size = 0.5, scale = 1.4, ...) {
  geom_density_ridges_gradient(
    rel_min_height = rel_min_height, 
    size = size, 
    scale = scale, 
    ...
  )
}
