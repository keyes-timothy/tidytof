# flowsom_metacluster_helpers.R
#
# This file contains helper functions that adapt the usage of the build-in
# metaclustering functions of the FlowSOM package. Due to a few small bugs
# in the current release, we copy-and-paste (with only small changes) some of their
# non-exported functions so that things can work as intended.
#
# For details, see https://github.com/saeyslab/FlowSOM and
# https://github.com/saeyslab/FlowSOM/issues/11

# functions that restore FlowSOM's metaclustering functions (other than
# consensus clustering, everything in their package throws a bug right now)
metaClustering_consensus <- function(data, k = 7, seed = NULL) {

  rlang::check_installed(pkg = "ConsensusClusterPlus")

  if (!rlang::is_installed(pkg = "ConsensusClusterPlus")) {
    stop("Consensus clustering requires the ConsensusClusterPlus package")
  }

  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1,
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = "hc", # "hc","km","kmdist","pam"
    distance = "euclidean" ,
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))

  results[[k]]$consensusClass
}

tof_metaClustering_hclust <-
  function(data, k = 7, seed) {
    d <- stats::dist(data, method = "minkowski")
    fit <- stats::hclust(d, method = "ward.D2")
    stats::cutree(fit, k = k)
  }

tof_metaClustering_kmeans <-
  function(data, k = 7, seed) {
    stats::kmeans(data, centers = k)$cluster
  }

tof_metaClustering_som <-
  function(data, k = 7, seed) {
    s <- suppressMessages(FlowSOM::SOM(data, xdim = k, ydim = 1, rlen = 100))
    result <- s$mapping[,1]
    return(result)
  }

## copied-and-pasted internal functions from the FlowSOM package that
## aren't exported and thus need to be included here

MetaClustering <- function(data, method, max = 20, seed = NULL, ...){
  res <- DetermineNumberOfClusters(data, max, method, seed = seed, ...)
  method <- get(method)
  method(data, k = res, seed = seed)
}

DetermineNumberOfClusters <-
  function(
    data,
    max,
    method,
    plot = FALSE,
    smooth = 0.2,
    seed = NULL,
    ...
  ){
    # Try out a clustering algorithm for several numbers of clusters and
    # select optimal
    #
    # Args:
    #     data:     Matrix containing the data to cluster
    #     max:        Maximum number of clusters to try
    #     method: Clustering method to use
    #     plot:     Whether to plot the results for different k
    #     smooth: Smoothing option to find elbow:
    #             0: no smoothing, 1: maximal smoothing
    #     seed:   Seed to pass on to given method
    #
    # Returns:
    #     Optimal number of clusters
    if(method ==  "metaClustering_consensus"){
      results <- consensus(data,max, seed, ...)
      res <- rep(0,max)
      res[1] <- SSE(data,rep(1,nrow(data)))
      for(i in 2:max){
        c <- results[[i]]$consensusClass
        res[i] <- SSE(data, c)
      }
    } else {
      method <- get(method)
      res <- rep(0, max)
      for(i in 1:max){
        c <- method(data, k = i,...)
        res[i] <- SSE(data, c)
      }
    }

    for(i in 2:(max - 1)){
      res[i] <- (1 - smooth) * res[i] +
        (smooth / 2) * res[i - 1] +
        (smooth / 2) * res[i + 1]
    }

    if(plot) plot(1:max, res, type = "b", xlab = "Number of Clusters",
                  ylab = "Within groups sum of squares")
    findElbow(res)
  }

findElbow <- function(data){
  n <- length(data)
  data <- as.data.frame(cbind(1:n,data))
  colnames(data) <- c("X","Y")

  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  optimal
}

consensus <- function(data, max, seed = NULL, ...){
  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = max, reps = 100, pItem = 0.9, pFeature = 1,
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = "hc", # "hc","km","kmdist","pam"
    distance = "euclidean",
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))
}

#' @importFrom methods is
SSE <- function(data, clustering){
  if(!methods::is(clustering, "numeric"))
    clustering <- as.numeric(as.factor(clustering))
  c_wss <- 0
  for(j in seq_along(clustering)){
    if(sum(clustering == j) > 1){
      c_wss <- c_wss + (nrow(data[clustering == j, , drop = FALSE]) - 1) *
        sum(apply(data[clustering == j, , drop = FALSE], 2, stats::var))
    }
  }
  c_wss
}
