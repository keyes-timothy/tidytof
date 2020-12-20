################   tof_cluster_flowSOM.R
### Description: 
### Clusters the rows of a tof_tibble using the flowSOM algorithm. 
#
# Inputs: 
#     - tof_tibble = a tibble, data.frame, or something that can be coerced into either 
#       of these
#     - cluster_markers = character vector of markers to use for clustering
#     - num_clusters = number of clusters you want in the final result
#
# Outputs: 
#     - flowSOM_metaclusters = a character vector of cluster membership for each row of `tof_tibble.` 
#
# Dependencies: 
#     - tidyverse
#     - flowSOM

tof_cluster_flowSOM <- 
  function(
    tof_tibble = NULL, 
    clustering_markers = NULL, 
    num_clusters = 20
  ) {
    my_FSO <- 
      list(
        data = 
          tof_tibble %>% 
          select(all_of(clustering_markers)) %>% 
          data.matrix(), 
        compensate = FALSE, 
        spillover = NULL, 
        transform = FALSE, 
        scale = NULL, 
        prettyColnames = clustering_markers
      )
    
    my_SOM <- 
      BuildSOM(
        fsom = my_FSO, 
        colsToUse = clustering_markers
      )
    
    #my_clusters <- my_SOM$map$mapping[,1]
    
    my_MST <- BuildMST(my_SOM, tSNE = FALSE)
    
    flowSOM_metacluster_object <- 
      MetaClustering(
        data = my_MST$map$codes, 
        method = "metaClustering_consensus", 
        nClus = num_clusters
      )
    
    flowSOM_metaclusters <- 
      flowSOM_metacluster_object[my_MST$map$mapping[,1]] %>% 
      as.character()
    
    return(flowSOM_metaclusters)
  }
