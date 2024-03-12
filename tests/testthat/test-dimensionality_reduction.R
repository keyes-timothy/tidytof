library(dplyr)
library(purrr)
library(tidytof)
library(testthat)

# setup

set.seed(2020)

data(ddpr_data)

ddpr <-
    ddpr_data |>
    group_by(sample_name) |>
    slice_sample(n = 100) |>
    ungroup()

# tof_reduce_pca ---------------------------------------------------------------

test_that("Shape of PCA result is correct", {
    dr_ddpr <-
        ddpr |>
        tof_reduce_pca()

    dr_ddpr_2 <-
        ddpr |>
        tof_reduce_pca(num_comp = 4L)

    expect_equal(nrow(ddpr), nrow(dr_ddpr))
    expect_equal(ncol(dr_ddpr), 5L)
    expect_equal(ncol(dr_ddpr_2), 4L)
})

test_that("colnames of PCA result are correct", {
    pca_colnames <-
        ddpr |>
        tof_reduce_pca() |>
        colnames()

    expect_true(all(grepl("^\\.pc", pca_colnames)))
})

test_that("return_recipe argument works", {
    pca_recipe <-
        ddpr |>
        tof_reduce_pca(return_recipe = TRUE)

    expect_s3_class(object = pca_recipe, "recipe")
})

# tof_reduce_tsne --------------------------------------------------------------

test_that("Shape of tsne result is correct", {
    dr_ddpr <-
        ddpr |>
        tof_reduce_tsne()

    dr_ddpr_2 <-
        ddpr |>
        tof_reduce_tsne(num_comp = 3L)

    expect_equal(nrow(ddpr), nrow(dr_ddpr))
    expect_equal(ncol(dr_ddpr), 2L)
    expect_equal(ncol(dr_ddpr_2), 3L)
})

test_that("colnames of tSNE result are correct", {
    tsne_colnames <-
        ddpr |>
        tof_reduce_tsne() |>
        colnames()

    expect_true(all(grepl("^\\.tsne", tsne_colnames)))
})


# tof_reduce_umap --------------------------------------------------------------
test_that("Shape of umap result is correct", {
    dr_ddpr <-
        ddpr |>
        tof_reduce_umap()

    dr_ddpr_2 <-
        ddpr |>
        tof_reduce_umap(num_comp = 3L)

    expect_equal(nrow(ddpr), nrow(dr_ddpr))
    expect_equal(ncol(dr_ddpr), 2L)
    expect_equal(ncol(dr_ddpr_2), 3L)
})

test_that("colnames of UMAP result are correct", {
    umap_colnames <-
        ddpr |>
        tof_reduce_umap() |>
        colnames()

    expect_true(all(grepl("^\\.umap", umap_colnames)))
})

test_that("return_recipe argument works", {
    umap_recipe <-
        ddpr |>
        tof_reduce_umap(return_recipe = TRUE)

    expect_s3_class(object = umap_recipe, "recipe")
})

# tof_reduce_dimensions --------------------------------------------------------

test_that("pca results are identical to tof_reduce_pca", {
    expect_equal(
        tof_reduce_pca(ddpr),
        tof_reduce_dimensions(ddpr, method = "pca", augment = FALSE)
    )
})

test_that("tsne results are identical to tof_reduce_tsne", {
    set.seed(2020)
    tsne_1 <- tof_reduce_tsne(ddpr)
    set.seed(2020)
    tsne_2 <- tof_reduce_dimensions(ddpr, method = "tsne", augment = FALSE)
    expect_equal(
        tsne_1,
        tsne_2
    )
})

test_that("umap results are identical to tof_reduce_umap", {
    set.seed(2020)
    umap_1 <- tof_reduce_umap(ddpr)
    set.seed(2020)
    umap_2 <- tof_reduce_dimensions(ddpr, method = "umap", augment = FALSE)
    expect_equal(
        umap_1,
        umap_2
    )
})

test_that("augment argument produces output with the correct shape", {
    ddpr_pca <- tof_reduce_dimensions(ddpr, method = "pca", num_comp = 3)

    expect_equal(ncol(ddpr_pca), ncol(ddpr) + 3L)
    expect_equal(nrow(ddpr_pca), nrow(ddpr))
})

test_that("invalid method throws an error", {
    expect_error(tof_reduce_dimensions(ddpr, method = "foo"))
})
