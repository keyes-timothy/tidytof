library(tidytof)
library(testthat)

# setup

tt <-
  ddpr_data %>%
  tof_preprocess() %>%
  tof_cluster(
    cluster_cols = starts_with("CD"),
    group_cols = sample_name,
    perform_metaclustering = FALSE,
    method = "flowsom"
  ) %>%
  tof_metacluster(
    cluster_col = .flowsom_cluster,
    metacluster_cols = starts_with("CD"),
    method = "phenograph"
  ) %>%
  tof_metacluster(
    cluster_col = .flowsom_cluster,
    metacluster_cols = starts_with("CD"),
    num_metaclusters = 3,
    method = "kmeans"
  ) %>%
  tof_downsample(
    group_cols = .phenograph_metacluster,
    num_cells = 500,
    method = "constant"
  ) %>%
  tof_reduce_dimensions(tsne_cols = starts_with("CD"), method = "tsne") %>%
  tof_reduce_dimensions(
    pca_cols = starts_with("CD"),
    num_comp = 2L,
    method = "pca"
  ) %>%
  mutate(
    replicate = sample(c("a", "b", "c", "d"), size = n(), replace = TRUE),
    sample_name_2 =
      if_else(
        sample_name == "Healthy1_Basal",
        sample(c("x", "y"), size = n(), replace = TRUE),
        "z"
      )
  )

daa_ttest <-
  tt %>%
  tof_daa_ttest(
    cluster_col = .kmeans_metacluster,
    effect_col = sample_name,
    group_cols = replicate,
    min_samples = 2
  )

daa_diffcyt <-
  tt %>%
  mutate(sample = paste0(sample_name, replicate)) %>%
  tof_daa_diffcyt(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    fixed_effect_cols = sample_name,
    diffcyt_method = "voom",
    min_samples = 2
  )

daa_diffcyt_2 <-
  tt %>%
  mutate(sample = paste0(sample_name_2, replicate)) %>%
  tof_daa_diffcyt(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    fixed_effect_cols = sample_name_2,
    diffcyt_method = "voom",
    min_samples = 2
  )

daa_glmm <-
  tt %>%
  mutate(sample = paste0(sample_name, replicate)) %>%
  tof_daa_glmm(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    fixed_effect_cols = sample_name,
    min_samples = 2
  )

dea_ttest <-
 suppressWarnings(
   tt %>%
     tof_dea_ttest(
       cluster_col = .kmeans_metacluster,
       effect_col = sample_name,
       group_cols = replicate,
       min_samples = 2
     )
 )

dea_lmm <-
  tt %>%
  mutate(sample = paste0(sample_name, replicate)) %>%
  tof_dea_lmm(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    fixed_effect_cols = sample_name,
    min_samples = 2
  )

dea_diffcyt <-
  tt %>%
  mutate(sample = paste0(sample_name, replicate)) %>%
  tof_dea_diffcyt(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = sample_name,
    diffcyt_method = "lmm",
    min_samples = 2
  )

dea_diffcyt_2 <-
  tt %>%
  mutate(sample = paste0(sample_name_2, replicate)) %>%
  tof_dea_diffcyt(
    sample_col = sample,
    cluster_col = .kmeans_metacluster,
    marker_cols = c(cd45, cd19),
    fixed_effect_cols = sample_name_2,
    diffcyt_method = "limma",
    min_samples = 2
  )

# tof_plot_cells_embedding -----------------------------------------------------

test_that("embedding visualizations run and return ggplot objects when embeddings are precomputed", {
  embed_plot_tsne <-
    tt %>%
    tof_plot_cells_embedding(
      embedding_cols = starts_with(".tsne"),
      color_col = .phenograph_metacluster
    )

  embed_plot_pca <-
    tt %>%
    tof_plot_cells_embedding(
      embedding_cols = starts_with(".pc"),
      color_col = .phenograph_metacluster
    )

  expect_s3_class(embed_plot_tsne, "ggplot")
  expect_s3_class(embed_plot_pca, "ggplot")
})

test_that("embedding visualizations run and return ggplot objects when embeddings are computed de novo", {
  embed_plot_tsne <-
    tt %>%
    tof_plot_cells_embedding(
      color_col = .phenograph_metacluster,
      embedding_method = "tsne",
      tsne_cols = c(cd22, cd34, cd45, cd123)
    )

  embed_plot_pca <-
    embed_plot_tsne <-
    tt %>%
    tof_plot_cells_embedding(
      color_col = .phenograph_metacluster,
      embedding_method = "pca",
      pca_cols = c(cd22, cd34, cd45, cd123)
    )

  embed_plot_umap <-
    tt %>%
    tof_plot_cells_embedding(
      color_col = .phenograph_metacluster,
      embedding_method = "umap",
      umap_cols = c(cd22, cd34, cd45, cd123)
    )

  expect_s3_class(embed_plot_tsne, "ggplot")
  expect_s3_class(embed_plot_pca, "ggplot")
  expect_s3_class(embed_plot_umap, "ggplot")
})

# tof_plot_cells_layout --------------------------------------------------------


test_that("layout visualizations run and return ggplot objects", {

  layout_plot <-
    tt %>%
    tof_plot_cells_layout(
      knn_cols = starts_with("CD"),
      color_col = .phenograph_metacluster,
      num_neighbors = 5L,
      graph_type = "unweighted"
    )

  expect_s3_class(layout_plot, "ggplot")

})

test_that("layout visualizations can use old plots as a template", {

  layout_plot <-
    tt %>%
    tof_downsample_constant(group_cols = .phenograph_metacluster, num_cells = 20) %>%
    tof_plot_cells_layout(
      knn_cols = starts_with("CD"),
      color_col = .phenograph_metacluster,
      num_neighbors = 5L,
      graph_type = "unweighted"
    )

  expect_s3_class(layout_plot, "ggplot")

})





# tof_plot_cells_density -------------------------------------------------------

test_that("density plots run and result in ggplot objects", {

  density_1 <-
    tt %>%
    tof_plot_cells_density(
      marker_col = cd43,
      group_col = .kmeans_metacluster,
      use_ggridges = TRUE
    )

  density_2 <-
    tt %>%
    tof_plot_cells_density(
      marker_col = ps6,
      group_col = .kmeans_metacluster,
      use_ggridges = FALSE
    )

  density_3 <-
    tt %>%
    tof_plot_cells_density(
      marker_col = cd45,
      use_ggridges = FALSE
    )

  expect_s3_class(density_1, "ggplot")
  expect_s3_class(density_2, "ggplot")
  expect_s3_class(density_3, "ggplot")
})


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
        color_col = .phenograph_metacluster
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = .phenograph_metacluster,
        node_size = "cluster_size"
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        group_cols = .phenograph_metacluster,
        color_col = .phenograph_metacluster,
        graph_type = "unweighted",
        node_size = "cluster_size"
      ),
    "ggplot"
  )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = .phenograph_metacluster,
        graph_type = "unweighted"
      ),
    "ggplot"
  )

  # use an existing plot as a template layout
  layout_plot <-
    tt %>%
    tof_plot_cluster_mst(
      cluster_col = .flowsom_cluster,
      color_col = .phenograph_metacluster,
      graph_type = "unweighted"
    )

  expect_s3_class(
    tt %>%
      tof_plot_cluster_mst(
        cluster_col = .flowsom_cluster,
        color_col = cd34,
        graph_type = "unweighted",
        graph_layout = layout_plot
      ),
    "ggplot"
  )

})

# tof_plot_cluster_heatmap ---------------------------------------------------------

test_that("heatmap function runs and returns ggplot object", {
  heatmap_1 <-
    tt %>%
    tof_plot_cluster_heatmap(
      marker_cols = starts_with("cd"),
      cluster_col = .flowsom_cluster
    )
  expect_s3_class(heatmap_1, "ggplot")

  heatmap_2 <-
    tt %>%
    tof_plot_cluster_heatmap(
      marker_cols = starts_with("cd"),
      cluster_col = .flowsom_cluster,
      scale_markerwise = TRUE,
      scale_clusterwise = TRUE
    )
  expect_s3_class(heatmap_2, "ggplot")


  heatmap_3 <-
    tt %>%
    tof_plot_cluster_heatmap(
      marker_cols = starts_with("cd"),
      cluster_col = .flowsom_cluster,
      scale_markerwise = TRUE,
      scale_clusterwise = TRUE
    )
  expect_s3_class(heatmap_3, "ggplot")


  heatmap_4 <-
    tt %>%
    tof_plot_cluster_heatmap(
      marker_cols = starts_with("cd"),
      cluster_col = .phenograph_metacluster,
      scale_markerwise = FALSE,
      scale_clusterwise = FALSE
    )
  expect_s3_class(heatmap_4, "ggplot")


  heatmap_5 <-
    tt %>%
    tof_plot_heatmap(
      y_col = sample_name
    )
  expect_s3_class(heatmap_5, "ggplot")


  heatmap_6 <-
    tt %>%
    tof_plot_sample_heatmap(
      sample_col = sample_name_2
    )
  expect_s3_class(heatmap_6, "ggplot")

})

# tof_plot_cluster_volcano -----------------------------------------------------

test_that("tof_plot_cluster_volcano throws an error if diffcyt_lmm was the dea_method", {
  expect_error(
    dea_diffcyt %>%
      tof_plot_cluster_volcano(
      )
  )
})

test_that("tof_plot_cluster_volcano runs and creates ggplot objects", {
  volcano_diffcyt_2 <-
    dea_diffcyt_2 %>%
    tof_plot_cluster_volcano(point_size = 4, use_ggrepel = FALSE)
  expect_s3_class(volcano_diffcyt_2, "ggplot")

  volcano_lmm <-
    suppressWarnings(
      dea_lmm %>%
        tof_plot_cluster_volcano(use_ggrepel = TRUE)
    )
  expect_s3_class(volcano_lmm, "ggplot")

  dea_ttest <-
    suppressWarnings(
      tt %>%
        select(-starts_with(".tsne"), -starts_with(".pc")) %>%
        tof_dea_ttest(
          cluster_col = .kmeans_metacluster,
          marker_cols = where(tof_is_numeric), #starts_with("CD"),
          effect_col = sample_name,
          group_cols = replicate,
          min_samples = 2
        )
    )

  volcano_ttest <-
    suppressWarnings(
      dea_ttest %>%
        tof_plot_cluster_volcano(use_ggrepel = TRUE)
    )
  expect_s3_class(volcano_ttest, "ggplot")

})


# tof_plot_cluster_abundance ---------------------------------------------------

test_that("tof_plot_cluster_abundance runs and returns ggplot objects", {
  expect_s3_class(
    daa_ttest %>%
      tof_plot_cluster_abundance(
        tof_tibble = tt,
        cluster_col = .kmeans_metacluster,
        group_cols = c(replicate, sample_name),
        color_col = sample_name
      ),
    "ggplot"
  )

  expect_s3_class(
    daa_glmm %>%
      tof_plot_cluster_abundance(
        tof_tibble = tt,
        cluster_col = .kmeans_metacluster,
        group_cols = c(replicate, sample_name),
        color_col = sample_name
      ),
    "ggplot"
  )
})

# tof_plot_sample_model --------------------------------------------------------

feature_tibble <-
  dplyr::tibble(
    sample = as.character(1:100),
    cd45 = runif(n = 100),
    pstat5 = runif(n = 100),
    cd34 = runif(n = 100),
    outcome = (3 * cd45) + (4 * pstat5) + rnorm(100),
    class =
      if_else(outcome > median(outcome), "class1", "class2") %>%
      as.factor(),
    multiclass =
      if_else(class == "class1", "class1", sample(c("class2", "class3"), size = 100, replace = TRUE)) %>%
      as.factor(),
       # as.factor(c(rep("class1", 30), rep("class2", 30), rep("class3", 40))),
    event = c(rep(0, times = 50), rep(1, times = 50)),
    time_to_event =
      if_else(
        outcome > median(outcome),
        rnorm(n = 100, mean = 8, sd = 2),
        rnorm(n = 100, mean = 14, sd = 2)
      )
  )

# permute all columns to get "new" data
new_data <-
  feature_tibble %>%
  mutate(across(everything(), ~ sample(x = .x, size = length(.x), replace = FALSE)))

linear_tof_model <-
  feature_tibble %>%
  tof_split_data(split_method = "bootstrap", num_bootstraps = 1L) %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = outcome,
    model_type = "linear"
  )

logistic_tof_model <-
  feature_tibble %>%
  tof_split_data(split_method = "bootstrap", num_bootstraps = 1L) %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = class,
    model_type = "two-class"
  )

multinomial_tof_model <-
  feature_tibble %>%
  tof_split_data(split_method = "bootstrap", num_bootstraps = 1L) %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    response_col = multiclass,
    model_type = "multiclass"
  )

survival_tof_model <-
  feature_tibble %>%
  tof_split_data(split_method = "bootstrap", num_bootstraps = 1L) %>%
  tof_train_model(
    predictor_cols = c(cd45, pstat5, cd34),
    event_col = event,
    time_col = time_to_event,
    model_type = "survival"
  )

linear_plot <-
  linear_tof_model %>%
  tof_plot_model()

linear_plot_new <-
  linear_tof_model %>%
  tof_plot_model(new_data = new_data)

linear_plot_tuning <-
  linear_tof_model %>%
  tof_plot_model(new_data = "tuning")

logistic_plot <-
  logistic_tof_model %>%
  tof_plot_model()

logistic_plot_new <-
  logistic_tof_model %>%
  tof_plot_model(new_data = new_data)

logistic_plot_tuning <-
  logistic_tof_model %>%
  tof_plot_model(new_data = "tuning")

multinomial_plot <-
  multinomial_tof_model %>%
  tof_plot_model()

multinomial_plot_new <-
  multinomial_tof_model %>%
  tof_plot_model(new_data = new_data)

multinomial_plot_tuning <-
  multinomial_tof_model %>%
  tof_plot_model(new_data = "tuning")

survival_plot <-
  survival_tof_model %>%
  tof_plot_model()

survival_plot_new <-
  survival_tof_model %>%
  tof_plot_model(new_data = new_data)

survival_plot_tuning <-
  survival_tof_model %>%
  tof_plot_model(new_data = "tuning")

test_that("all results are ggplot objects", {
  expect_s3_class(linear_plot, "ggplot")
  expect_s3_class(logistic_plot, "ggplot")
  expect_s3_class(multinomial_plot, "ggplot")
  expect_s3_class(survival_plot, "ggplot")

  expect_s3_class(linear_plot_new, "ggplot")
  expect_s3_class(logistic_plot_new, "ggplot")
  expect_s3_class(multinomial_plot_new, "ggplot")
  expect_s3_class(survival_plot_new, "ggplot")

  expect_s3_class(linear_plot_tuning, "ggplot")
  expect_s3_class(logistic_plot_tuning, "ggplot")
  expect_s3_class(multinomial_plot_tuning, "ggplot")
  expect_s3_class(survival_plot_tuning, "ggplot")
})



