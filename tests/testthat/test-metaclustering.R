library(dplyr)
library(tidytof)
library(testthat)

data(phenograph_data)

# setup
clust_data <-
  phenograph_data %>%
  tof_preprocess() %>%
  tof_cluster(method = "kmeans")

# tof_metacluster_hierarchical -------------------------------------------------
test_that("metaclustering result is the right shape", {
  result <-
    clust_data %>%
    tof_metacluster_hierarchical(
      cluster_col = .kmeans_cluster
    )

  expect_equal(nrow(result), nrow(clust_data))
  expect_equal(ncol(result), 1L)
})

test_that("metaclustering result has the right name", {
  result <-
    clust_data %>%
    tof_metacluster_hierarchical(
      cluster_col = .kmeans_cluster
    )

  expect_equal(colnames(result), ".hierarchical_metacluster")
})



# tof_metacluster_kmeans -------------------------------------------------------
test_that("metaclustering result is the right shape", {
  result <-
    clust_data %>%
    tof_metacluster_kmeans(
      cluster_col = .kmeans_cluster
    )

  expect_equal(nrow(result), nrow(clust_data))
  expect_equal(ncol(result), 1L)
})

test_that("metaclustering result has the right name", {
  result <-
    clust_data %>%
    tof_metacluster_kmeans(
      cluster_col = .kmeans_cluster
    )

  expect_equal(colnames(result), ".kmeans_metacluster")
})


# tof_metacluster_phenograph ---------------------------------------------------
test_that("metaclustering result is the right shape", {
  result <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster
    )

  expect_equal(nrow(result), nrow(clust_data))
  expect_equal(ncol(result), 1L)
})

test_that("metaclustering result has the right name", {
  result <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster
    )

  expect_equal(colnames(result), ".phenograph_metacluster")
})

test_that("metaclustering result works with various kinds of tidyselection", {

  result_1 <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster,
      metacluster_cols = contains("cd")
    )

  result_2 <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster,
      metacluster_cols = c(contains("cd"), -cd45)
    )

  result_3 <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster,
      metacluster_cols = c(cd45)
    )

  result_4 <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster,
      metacluster_cols = c(cd45, -sample_name)
    )

  expect_equal(nrow(result_1), nrow(clust_data))
  expect_equal(ncol(result_1), 1L)
  expect_equal(nrow(result_2), nrow(clust_data))
  expect_equal(ncol(result_2), 1L)
  expect_equal(nrow(result_3), nrow(clust_data))
  expect_equal(ncol(result_3), 1L)
  expect_equal(nrow(result_4), nrow(clust_data))
  expect_equal(ncol(result_4), 1L)
})

# tof_metacluster_consensus ----------------------------------------------------
test_that("metaclustering result is the right shape", {
  result <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "kmeans"
    )

  result_2 <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "pam"
    )

  result_3 <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "hierarchical"
    )

  result_4 <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "hierarchical",
      distance_function = "pearson"
    )

  expect_equal(nrow(result), nrow(clust_data))
  expect_equal(ncol(result), 1L)

  expect_equal(nrow(result_2), nrow(clust_data))
  expect_equal(ncol(result_2), 1L)

  expect_equal(nrow(result_3), nrow(clust_data))
  expect_equal(ncol(result_3), 1L)

  expect_equal(nrow(result_4), nrow(clust_data))
  expect_equal(ncol(result_4), 1L)
})

test_that("metaclustering result has the right name", {
  result <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster
    )

  expect_equal(colnames(result), ".consensus_metacluster")
})


# tof_metacluster_flowsom ------------------------------------------------------
test_that("metaclustering result is the right shape", {

  result <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "consensus"
    )

  result_2 <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "hierarchical"
    )

  result_3 <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "som"
    )

  result_4 <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster,
      clustering_algorithm = "kmeans"
    )

  expect_equal(nrow(result), nrow(clust_data))
  expect_equal(ncol(result), 1L)

  expect_equal(nrow(result_2), nrow(clust_data))
  expect_equal(ncol(result_2), 1L)

  expect_equal(nrow(result_3), nrow(clust_data))
  expect_equal(ncol(result_3), 1L)

  expect_equal(nrow(result_4), nrow(clust_data))
  expect_equal(ncol(result_4), 1L)
})

test_that("metaclustering result has the right name", {
  result <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster
    )

  expect_equal(colnames(result), ".flowsom_metacluster")
})

# tof_metacluster --------------------------------------------------------------

test_that("hclust tof_metacluster results are identical to subroutine", {
  hclust_sub <-
    clust_data %>%
    tof_metacluster_hierarchical(
      cluster_col = .kmeans_cluster
    )

  hclust_wrap <-
    clust_data %>%
    tof_metacluster(
      cluster_col = .kmeans_cluster,
      augment = FALSE,
      method = "hierarchical"
    )

  expect_equal(hclust_sub, hclust_wrap)
})

test_that("kmeans tof_metacluster results are identical to subroutine", {
  set.seed(2020)
  kmeans_sub <-
    clust_data %>%
    tof_metacluster_kmeans(
      cluster_col = .kmeans_cluster
    )

  set.seed(2020)
  kmeans_wrap <-
    clust_data %>%
    tof_metacluster(
      cluster_col = .kmeans_cluster,
      augment = FALSE,
      method = "kmeans"
    )

  expect_equal(kmeans_sub, kmeans_wrap)
})

test_that("phenograph tof_metacluster results are identical to subroutine", {
  set.seed(2020)
  pheno_sub <-
    clust_data %>%
    tof_metacluster_phenograph(
      cluster_col = .kmeans_cluster
    )

  set.seed(2020)
  pheno_wrap <-
    clust_data %>%
    tof_metacluster(
      cluster_col = .kmeans_cluster,
      augment = FALSE,
      method = "phenograph"
    )

  expect_equal(pheno_sub, pheno_wrap)
})

test_that("consensus tof_metacluster results are identical to subroutine", {
  ccp_sub <-
    clust_data %>%
    tof_metacluster_consensus(
      cluster_col = .kmeans_cluster,
      seed = 2020L
    )

  ccp_wrap <-
    clust_data %>%
    tof_metacluster(
      cluster_col = .kmeans_cluster,
      seed = 2020L,
      augment = FALSE,
      method = "consensus"
    )

  expect_equal(ccp_sub, ccp_wrap)
})

test_that("flowsom tof_metacluster results are identical to subroutine", {
  set.seed(2020)
  flow_sub <-
    clust_data %>%
    tof_metacluster_flowsom(
      cluster_col = .kmeans_cluster,
    )

  set.seed(2020)
  flow_wrap <-
    clust_data %>%
    tof_metacluster(
      cluster_col = .kmeans_cluster,
      augment = FALSE,
      method = "flowsom"
    )

  expect_equal(flow_sub, flow_wrap)
})




