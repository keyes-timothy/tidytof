#extract_featureMatrix.R
#
#function: takes a tidy matrix of clustered single-cell measurements and extracts cluster features for use 
#          in machine learning models. Returns a list of two feature matrices similar to that calculated by 
#          CITRUS (one is in long-data format, the other is in short-data format.)
#
#input: 1) expression.matrix: an [m x n] matrix of m cells and n variables (a combination of protein measurements
#          and tag variables that describe the patient, condition, stim, etc. that each cell belongs to. 
#       2) surface_markers: a character vector of all protein measurements that are NOT signaling markers 
#       3) signaling_markers: a character vector of all protein measurements that ARE signaling markers
#       4) grouping_vars: a character vector containing the names of all grouping variables over which to 
#          compute summary statistics. Passed to dplyr::group_by. 
#       5) cluster_var: a string giving the name of the column of expression.matrix in which cluster 
#          ids are located
#       6) stim_var: a string giving the name of the column of expression.matrix in which the stimulation
#          information for each cell is located.
#       7) central_tendency: a function indicating what summary statistic you'd like to use
#          when cluster features are calculated. The default is the mean. 
#
#output: 1) final_features, a list of two feature matrices. 
#             a) long_data: the first element of final_features. It is an [m x n] feature matrix of m samples and n 
#                           variables, where n contains both the calculated features from measured cyTOF parameters 
#                           and the grouping_vars as well as the cluster_var. 
#             b) wide_data: the second element of final_features. It is an [m x n] 
#        a feature matrix for all unique groups defined by the group_by variables. Features include the 
#        abundance of each population and the mean or median of each population for each marker in 
#        expression.matrix
#
#Dependencies: tidyverse, rlang

#helper function that will calculate the % of values in a vector over the inputted threshold value 
DDPR_thresholding <- function(x, threshold = 10){ 
  return(sum(x>threshold)/length(x))
}

#helper function that calculates the relative abundance of each sample in each cluster 
find_abundances <- function(expression_matrix = NULL, cluster_var = NULL){ 
  cells.clusters <- expression_matrix %>% group_by(!!sym(cluster_var), add = TRUE) %>% summarize(num.cells = n())
  cells.groups <- expression_matrix %>% summarize(total.cells = n())
  if (is.null(groups(expression_matrix))) { 
    cells.clusters <- cells.clusters %>% mutate(cell.percentage = num.cells/cells.groups %>% as.numeric())
    my.abundance <- cells.clusters %>% select(one_of(c(cluster_var, "cell.percentage")))
  } 
  else { 
    my.abundance <- left_join(x = cells.clusters, y = cells.groups) %>% 
      dplyr::mutate(cell.percentage = num.cells/total.cells) %>% 
      dplyr::select(cell.percentage, !!sym(cluster_var))
  }
}

#helper function
find_centroids <- 
  function(
    expression_matrix, 
    cluster_var, 
    surface_markers, 
    signaling_markers, 
    central_tendency, 
    ...
  ) { 
    my.centroids <- 
      expression_matrix %>% 
      dplyr::group_by(!!sym(cluster_var), add = TRUE) %>% 
      summarize_at(
        .funs = central_tendency, 
        .vars = vars(one_of(c(surface_markers, signaling_markers)))
      )
    return(my.centroids)
  } 

#helper function that will find the difference between each stimulation condition and the basal 
#condition for a given sample 
find_stim_features <- 
  function(
    expression_matrix = NULL, 
    cluster_var = NULL, 
    stim_var = NULL, 
    signaling_markers = NULL,
    central_tendency = NULL, 
    ...
  ) { 
    stim.summaries <- 
      expression_matrix %>% 
      dplyr::group_by(!!sym(stim_var), !!sym(cluster_var), add = TRUE) %>% 
      dplyr::summarize_at(signaling_markers,.funs = central_tendency)
    
    basal.summaries <- 
      stim.summaries %>% 
      dplyr::filter(!!sym(stim_var) == "Basal") %>% 
      ungroup() %>% dplyr::select(-stimulation)
    
    my.signals <- which(colnames(basal.summaries) %in% signaling_markers)
    names(basal.summaries)[my.signals] <- 
      paste0(colnames(basal.summaries)[my.signals], "_basal")
    
    stim.summaries <- stim.summaries %>% left_join(y = basal.summaries) %>% group_by(!!sym(cluster_var), add = TRUE) %>% 
      mutate_all(coalesce, 0) #replace NAs with 0 - the grouping beforehand is not optimal (hack-y way of not getting an error)
    for (signal in signaling_markers) {
      stim.summaries[signal] <- 
        stim.summaries[signal] - stim.summaries[paste0(signal, "_basal")]
    }
    stim.summaries <- stim.summaries %>% dplyr::select(-contains("_basal"))
    
    #find locations of signaling channels in colnames and change names 
    my.signals <- which(colnames(stim.summaries) %in% signaling_markers)
    names(stim.summaries)[my.signals] <- 
      paste0("delta_", colnames(stim.summaries)[my.signals])
    
    stim.summaries <- 
      stim.summaries %>% 
      ungroup() %>% 
      gather(key = signal, value = `change from basal`, contains("delta_")) %>% 
      mutate(stim.signal = paste(signal, !!sym(stim_var), sep = "@")) %>% 
      dplyr::select(-signal, -stimulation) %>% 
      tidyr::spread(key = stim.signal, value = `change from basal`, fill = 0) 
    return(stim.summaries)
  }

#actual function
extract_featureMatrix <- 
  function(expression.matrix = NULL, 
           surface_markers = NULL, 
           signaling_markers = NULL,
           grouping_vars = NULL, 
           cluster_var = NULL, 
           stim_var = NULL, 
           central_tendency = mean, 
           ...
  ) {
    
    expression.matrix <- 
      expression.matrix %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(!!!syms(grouping_vars))
    
    abundances <- find_abundances(expression.matrix, cluster_var)
    centroids <- 
      find_centroids(
        expression.matrix, 
        cluster_var, 
        surface_markers, 
        signaling_markers,
        central_tendency, 
        ...
      )
    
    if (!is.null(stim_var)) {
      stim_features <- 
        find_stim_features(
          expression.matrix, 
          cluster_var, 
          stim_var, 
          signaling_markers, 
          central_tendency, 
          ...
        )
      my_features <- left_join(abundances, centroids) %>% left_join(stim_features)
    }
    else{
      my_features <- left_join(abundances, centroids)
    }
    if (is.null(grouping_vars)) {
      final_features <- list(long_data = my_features, wide_data = NULL)
      return(final_features)
    }
    else{
      feature_matrix <- my_features %>% ungroup() %>%
        gather(
          key = "variable", 
          value = "value", 
          -one_of(c(grouping_vars, cluster_var))
        ) %>% 
        mutate(
          cluster.variable =  paste(!!sym(cluster_var), variable, sep = "_")
        ) %>% 
        dplyr::select(one_of(grouping_vars), cluster.variable, value) %>% 
        tidyr::spread(key = cluster.variable, value = value, fill = 0)
      
      final_features <- list(long_data = my_features, wide_data = feature_matrix)
      
      return(final_features)
    }
  }


#test
# my.test <- extract_featureMatrix(expression.matrix = expression.matrix, 
#                                  surface_markers = c(SURFACE.MARKERS, OTHER.MARKERS), 
#                                  signaling_markers = SIGNALING.MARKERS, 
#                                  grouping_vars = c("patient", "condition"), 
#                                  cluster_var = "phenograph.metacluster", 
#                                  stim_var = "stimulation", 
#                                  central_tendency = DDPR_thresholding, threshold = 0.1)





