library(tidytof)

# for interoperability with Bioconductor ---------------------------------------

data(ddpr_data)



## SingleCellExperiment -----------

ddpr_sce <-
  ddpr_data |>
  as_SingleCellExperiment()

ddpr_sce_dr <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  as_SingleCellExperiment(
    reduced_dimensions_cols = starts_with(".pc")
  )

ddpr_sce_dr_auto <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  as_SingleCellExperiment()

ddpr_sce_dr_auto_multiple <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  tof_reduce_dimensions(method = "umap") |>
  as_SingleCellExperiment()

ddpr_sce_dr_split <-
  ddpr_sce_dr_auto_multiple |>
  tof_split_tidytof_reduced_dimensions()

test_that("SCEs have the right number of rows", {
  expect_equal(nrow(ddpr_sce), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_sce_dr), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_sce_dr_auto), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_sce_dr_auto_multiple), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_sce_dr_split), ncol(ddpr_data) - 1)
})

test_that("SCEs have the right number of columns", {
  expect_equal(ncol(ddpr_sce), nrow(ddpr_data))
  expect_equal(ncol(ddpr_sce_dr), nrow(ddpr_data))
  expect_equal(ncol(ddpr_sce_dr_auto), nrow(ddpr_data))
  expect_equal(ncol(ddpr_sce_dr_auto_multiple), nrow(ddpr_data))
  expect_equal(ncol(ddpr_sce_dr_split), nrow(ddpr_data))
})


test_that("SCEs have the right marker structure", {
  expect_setequal(
    names(ddpr_sce@rowRanges),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    names(ddpr_sce_dr@rowRanges),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    names(ddpr_sce_dr_auto@rowRanges),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    names(ddpr_sce_dr_auto_multiple@rowRanges),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    names(ddpr_sce_dr_split@rowRanges),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
})


test_that("SCEs have the right dimensionality reduction structure", {
  # rows
  expect_equal(
    nrow(ddpr_sce@int_colData$reducedDims),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_sce_dr@int_colData$reducedDims),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_sce_dr_auto@int_colData$reducedDims),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_sce_dr_auto_multiple@int_colData$reducedDims),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_sce_dr_split@int_colData$reducedDims),
    nrow(ddpr_data)
  )

  # columns
  expect_equal(
    ncol(ddpr_sce@int_colData$reducedDims),
    0L
  )

  expect_equal(
    ncol(ddpr_sce_dr@int_colData$reducedDims),
    1L
  )
  expect_equal(
    ncol(ddpr_sce_dr@int_colData$reducedDims$tidytof_reduced_dimensions),
    5L
  )

  expect_equal(
    ncol(ddpr_sce_dr_auto@int_colData$reducedDims),
    1L
  )
  expect_equal(
    ncol(ddpr_sce_dr_auto@int_colData$reducedDims$tidytof_reduced_dimensions),
    5L
  )

  expect_equal(
    ncol(ddpr_sce_dr_auto_multiple@int_colData$reducedDims),
    1L
  )
  expect_equal(
    ncol(ddpr_sce_dr_auto_multiple@int_colData$reducedDims$tidytof_reduced_dimensions),
    7L
  )

  expect_equal(
    ncol(ddpr_sce_dr_split@int_colData$reducedDims),
    2L
  )

})

test_that("SCEs with more than one dimensionality reduction embedding can be split.", {
  expect_setequal(
    colnames(ddpr_sce_dr_split@int_colData$reducedDims),
    c("tidytof_pca", "tidytof_umap")
  )
})


## SeuratObjects ----------


ddpr_seurat <-
  ddpr_data |>
  as_seurat()

ddpr_seurat_dr <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  as_seurat(
    reduced_dimensions_cols = starts_with(".pc")
  )

ddpr_seurat_dr_auto <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  as_seurat()

ddpr_seurat_dr_auto_multiple <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  tof_reduce_dimensions(method = "umap") |>
  as_seurat()

ddpr_seurat_dr_split <-
  ddpr_data |>
  tof_reduce_dimensions(num_comp = 5L, method = "pca") |>
  tof_reduce_dimensions(method = "umap") |>
  as_seurat(split_reduced_dimensions = TRUE)



test_that("SeuratObjects have the right number of rows", {
  expect_equal(nrow(ddpr_seurat), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_seurat_dr), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_seurat_dr_auto), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_seurat_dr_auto_multiple), ncol(ddpr_data) - 1)
  expect_equal(nrow(ddpr_seurat_dr_split), ncol(ddpr_data) - 1)
})


test_that("SeuratObjects have the right number of columns", {
  expect_equal(ncol(ddpr_seurat), nrow(ddpr_data))
  expect_equal(ncol(ddpr_seurat_dr), nrow(ddpr_data))
  expect_equal(ncol(ddpr_seurat_dr_auto), nrow(ddpr_data))
  expect_equal(ncol(ddpr_seurat_dr_auto_multiple), nrow(ddpr_data))
  expect_equal(ncol(ddpr_seurat_dr_split), nrow(ddpr_data))
})


test_that("SeuratObjects have the right marker structure", {
  expect_setequal(
    rownames(ddpr_seurat@assays$cytometry),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    rownames(ddpr_seurat_dr@assays$cytometry),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    rownames(ddpr_seurat_dr_auto@assays$cytometry),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    rownames(ddpr_seurat_dr_auto_multiple@assays$cytometry),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )
  expect_setequal(
    rownames(ddpr_seurat_dr_split@assays$cytometry),
    colnames(ddpr_data[-which(colnames(ddpr_data) == "sample_name")])
  )

})


test_that("SeuratObjects have the right dimensionality reduction structure", {
  # rows
  expect_equal(ddpr_seurat@reductions, list())
  expect_equal(
    nrow(ddpr_seurat_dr@reductions$tidytof_reduced_dimensions),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_seurat_dr_auto@reductions$tidytof_reduced_dimensions),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_seurat_dr_auto_multiple@reductions$tidytof_reduced_dimensions),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_seurat_dr_split@reductions$tidytof_pca),
    nrow(ddpr_data)
  )
  expect_equal(
    nrow(ddpr_seurat_dr_split@reductions$tidytof_umap),
    nrow(ddpr_data)
  )

  # columns
  expect_equal(
    ncol(ddpr_seurat_dr@reductions$tidytof_reduced_dimensions),
    5L
  )
  expect_equal(
    ncol(ddpr_seurat_dr_auto@reductions$tidytof_reduced_dimensions),
    5L
  )
  expect_equal(
    ncol(ddpr_seurat_dr_auto_multiple@reductions$tidytof_reduced_dimensions),
    7L
  )
  expect_equal(
    ncol(ddpr_seurat_dr_split@reductions$tidytof_pca),
    5L
  )
  expect_equal(
    ncol(ddpr_seurat_dr_split@reductions$tidytof_umap),
    2L
  )

})

test_that("SeuratObjects with more than one dimensionality reduction embedding can be split.", {
  expect_equal(length(ddpr_seurat_dr_split@reductions), 2L)
})










