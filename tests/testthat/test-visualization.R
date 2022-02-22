library(tidytof)
library(testthat)

tt <-
  ddpr_data %>%
  tof_preprocess() %>%
  tof_cluster(
    cluster_cols = starts_with("CD"),
    perform_metaclustering = FALSE,
    method = "flowsom"
  ) %>%
  tof_metacluster(
    cluster_col = .flowsom_cluster,
    metacluster_cols = starts_with("CD"),
    num_metaclusters = 3,
    method = "kmeans"
  ) %>%
  tof_downsample(
    group_cols = .kmeans_metacluster,
    num_cells = 500,
    method = "constant"
  ) %>%
  tof_reduce_dimensions(tsne_cols = starts_with("CD"), method = "tsne") %>%
  tof_reduce_dimensions(pca_cols = starts_with("CD"), num_comp = 2L, method = "pca")

# tof_plot_cells_embedding -----------------------------------------------------

test_that("embedding visualizations run and return ggplot objects", {
  embed_plot_tsne <-
    tt %>%
    tof_plot_cells_embedding(
      embedding_cols = starts_with(".tsne"),
      color_col = .kmeans_metacluster
    )

  embed_plot_pca <-
    tt %>%
    tof_plot_cells_embedding(
      embedding_cols = starts_with(".pc"),
      color_col = .kmeans_metacluster
    )

  expect_s3_class(embed_plot_tsne, "ggplot")
  expect_s3_class(embed_plot_pca, "ggplot")
})

# tof_plot_cells_layout --------------------------------------------------------

# knn_graph <-
#   tof_make_knn_graph(
#     tof_tibble = tt,
#     knn_cols = starts_with("CD"),
#     num_neighbors = 5L,
#     graph_type = "unweighted"
#   )

test_that("layout visualizations run and return ggplot objects", {

  layout_plot <-
    tt %>%
    tof_plot_cells_layout(
      knn_cols = starts_with("CD"),
      color_col = .kmeans_metacluster,
      num_neighbors = 5L,
      graph_type = "unweighted"
    )

  expect_s3_class(layout_plot, "ggplot")

})



# tof_plot_cells_histogram -----------------------------------------------------




# tof_plot_cluster_mst ---------------------------------------------------------

test_that("layout visualizations run and return ggplot objects", {

expect_s3_class(
  tt %>%
    tof_plot_cluster_mst(
      cluster_col = .flowsom_cluster,
      color_col = sample_name
    ),
  "ggplot"
)

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = .kmeans_metacluster
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = .kmeans_metacluster,
        node_size = "cluster_size"
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        group_cols = .kmeans_metacluster,
        color_col = .kmeans_metacluster,
        graph_type = "unweighted",
        node_size = "cluster_size"
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = .kmeans_metacluster,
        graph_type = "unweighted"
      ),
    "ggplot"
  )
})

