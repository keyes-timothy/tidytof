library(dplyr)
library(stringr)
library(tidytof)
library(testthat)

data(ddpr_data)
data(phenograph_data)

# setup
clust_data <-
    phenograph_data |>
    group_by(sample_name) |>
    slice_sample(n = 100) |>
    ungroup()

synth_data <-
    dplyr::tibble(
        x = c(
            runif(n = 1000, min = 0, max = 1),
            runif(n = 1000, min = 1.5, max = 3)
        ),
        y = c(
            runif(n = 1000, min = -4, max = -3),
            runif(n = 1000, min = -5.5, max = -5)
        ),
        z = c(
            rnorm(n = 1000, mean = 1, sd = 0.5),
            rnorm(n = 1000, mean = 3, sd = 1)
        )
    )


# tof_cluster_flowsom ----------------------------------------------------------


test_that("flowsom result is a tibble with a single character vector column of correct length", {
    flowsom <-
        clust_data |>
        tof_cluster_flowsom(cluster_cols = c(cd34, cd45, cd123, cd11b))

    flowsom_2 <-
        clust_data |>
        tof_cluster_flowsom(
            cluster_cols = c(cd34, cd45, cd123, cd11b),
            perform_metaclustering = FALSE
        )

    expect_equal(nrow(clust_data), nrow(flowsom))
    expect_true(is.character(flowsom$.flowsom_metacluster))
    expect_equal(ncol(flowsom), 1L)

    expect_equal(nrow(clust_data), nrow(flowsom_2))
    expect_true(is.character(flowsom_2$.flowsom_cluster))
    expect_equal(ncol(flowsom_2), 1L)
})

test_that("som_xdim and som_ydim arguments work", {
    flowsom_10 <-
        clust_data |>
        tof_cluster_flowsom(
            cluster_cols = c(cd34, cd45, cd123, cd11b),
            som_xdim = 3,
            som_ydim = 3,
            perform_metaclustering = FALSE
        )

    flowsom_5 <-
        clust_data |>
        tof_cluster_flowsom(
            cluster_cols = c(cd34, cd45, cd123, cd11b),
            som_xdim = 2,
            som_ydim = 2,
            perform_metaclustering = FALSE
        )

    expect_equal(nrow(distinct(flowsom_10)), 9)
    expect_equal(nrow(distinct(flowsom_5)), 4)
})

# tof_cluster_phenograph -------------------------------------------------------

test_that("phenograph result is a tibble with a single character vector column of correct length", {
    phenograph <-
        clust_data |>
        tof_cluster_phenograph(cluster_cols = c(cd34, cd45, cd123, cd11b))

    expect_equal(nrow(clust_data), nrow(phenograph))
    expect_true(is.character(phenograph$.phenograph_cluster))
    expect_equal(ncol(phenograph), 1L)
})

test_that("phenograph result columns are named correctly", {
    phenograph <-
        clust_data |>
        tof_cluster_phenograph(cluster_cols = c(cd34, cd45, cd123, cd11b))

    expect_equal(colnames(phenograph), ".phenograph_cluster")
})

# tof_cluster_kmeans -----------------------------------------------------------

test_that("k-means result is a tibble with a single character vector column of correct length", {
    k <-
        clust_data |>
        tof_cluster_kmeans(cluster_cols = c(cd34, cd45, cd123, cd11b))

    expect_equal(nrow(clust_data), nrow(k))
    expect_true(is.character(k$.kmeans_cluster))
    expect_equal(ncol(k), 1L)
})

test_that("num_clusters argument works", {
    k_10 <-
        clust_data |>
        tof_cluster_kmeans(
            cluster_cols = c(cd34, cd45, cd123, cd11b),
            num_clusters = 10
        )

    k_5 <-
        clust_data |>
        tof_cluster_kmeans(
            cluster_cols = c(cd34, cd45, cd123, cd11b),
            num_clusters = 5
        )

    expect_equal(nrow(distinct(k_10)), 10)
    expect_equal(nrow(distinct(k_5)), 5)
})




# tof_cluster_ddpr -------------------------------------------------------------

test_that("ddpr result is a tibble with a single character vector column of correct length", {
    healthy <-
        ddpr_data |>
        dplyr::filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        dplyr::mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    # in series
    ddpr <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19)
        )

    # in parallel
    ddpr_2 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            parallel_cols = sample_name
        )

    expect_equal(nrow(ddpr_data), nrow(ddpr))
    expect_true(is.character(ddpr$.mahalanobis_cluster))
    expect_equal(ncol(ddpr), 1L)

    expect_equal(nrow(ddpr_data), nrow(ddpr_2))
    expect_true(is.character(ddpr_2$.mahalanobis_cluster))
    expect_equal(ncol(ddpr_2), 1L)
})

test_that("return_distances argument works", {
    healthy <-
        ddpr_data |>
        filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    # in series
    ddpr <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            return_distances = TRUE
        )

    # in parallel
    ddpr_2 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            parallel_cols = sample_name,
            return_distances = TRUE
        )

    expect_equal(ncol(ddpr), 3L)
    expect_equal(ncol(ddpr_2), 3L)
})

test_that("output columns are named correctly in series", {
    healthy <-
        ddpr_data |>
        filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    ddpr_1 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19)
        )

    ddpr_1b <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            return_distances = TRUE
        )

    ddpr_2 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            distance_function = "cosine"
        )

    ddpr_3 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            distance_function = "pearson"
        )

    expect_true(stringr::str_detect(colnames(ddpr_1), "\\.mahalanobis"))
    expect_true(all(stringr::str_detect(colnames(ddpr_1b), "\\.mahalanobis")))
    expect_true(stringr::str_detect(colnames(ddpr_2), "\\.cosine"))
    expect_true(stringr::str_detect(colnames(ddpr_3), "\\.pearson"))
})

test_that("output columns are named correctly in parallel", {
    healthy <-
        ddpr_data |>
        filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    ddpr_1 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            parallel_cols = sample_name
        )

    ddpr_1b <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            parallel_cols = sample_name,
            return_distances = TRUE
        )

    ddpr_2 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            distance_function = "cosine",
            parallel_cols = sample_name
        )

    ddpr_3 <-
        ddpr_data |>
        tof_cluster_ddpr(
            healthy_tibble = healthy,
            healthy_label_col = cluster,
            cluster_cols = c(cd34, cd45, cd19),
            distance_function = "pearson",
            parallel_cols = sample_name
        )

    expect_true(stringr::str_detect(colnames(ddpr_1), "\\.mahalanobis"))
    expect_true(all(stringr::str_detect(colnames(ddpr_1b), "\\.mahalanobis")))
    expect_true(stringr::str_detect(colnames(ddpr_2), "\\.cosine"))
    expect_true(stringr::str_detect(colnames(ddpr_3), "\\.pearson"))
})


# tof_cluster_xshift -----------------------------------------------------------


# tof_cluster_tibble -----------------------------------------------------------

test_that("clustering output is identical for all methods", {
    healthy <-
        ddpr_data |>
        filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_flowsom(clust_data, perform_metaclustering = FALSE)
        },
        {
            set.seed(2020)
            tof_cluster_tibble(clust_data, method = "flowsom", augment = FALSE, perform_metacluster = FALSE)
        }
    )

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_kmeans(clust_data)
        },
        {
            set.seed(2020)
            tof_cluster_tibble(clust_data, method = "kmeans", augment = FALSE)
        }
    )

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_phenograph(clust_data)
        },
        {
            set.seed(2020)
            tof_cluster_tibble(clust_data, augment = FALSE, method = "phenograph")
        }
    )

    expect_equal(
        tof_cluster_ddpr(
            ddpr_data,
            healthy_tibble = healthy,
            healthy_label_col = cluster
        ),
        tof_cluster_tibble(
            ddpr_data,
            method = "ddpr",
            augment = FALSE,
            healthy_tibble = healthy,
            healthy_label_col = cluster
        )
    )
})


# tof_cluster_grouped ----------------------------------------------------------

test_that("result has the correct number of columns", {
    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            method = "phenograph"
        ) |>
            ncol(),
        ncol(clust_data) + 1
    )

    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            augment = FALSE,
            method = "phenograph"
        ) |>
            ncol(),
        1
    )
})

test_that("result has the correct number of rows", {
    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            num_clusters = 3,
            method = "kmeans"
        ) |>
            nrow(),
        nrow(clust_data)
    )

    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            augment = FALSE,
            method = "phenograph"
        ) |>
            nrow(),
        nrow(clust_data)
    )
})

test_that("result has the correct number unique cluster ids", {
    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            num_clusters = 3,
            method = "kmeans"
        ) |>
            distinct(.kmeans_cluster) |>
            nrow(),
        3 * nrow(distinct(clust_data, sample_name))
    )

    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            num_clusters = 5,
            method = "kmeans"
        ) |>
            distinct(.kmeans_cluster) |>
            nrow(),
        5 * nrow(distinct(clust_data, sample_name))
    )

    expect_equal(
        tof_cluster_grouped(
            tof_tibble = clust_data,
            group_cols = sample_name,
            num_clusters = 5,
            method = "kmeans"
        ) |>
            distinct(.kmeans_cluster) |>
            nrow(),
        5 * nrow(distinct(clust_data, sample_name))
    )
})




# tof_cluster ------------------------------------------------------------------
test_that("clustering output is identical for all methods", {
    healthy <-
        ddpr_data |>
        filter(stringr::str_detect(sample_name, "Healthy"))
    healthy <-
        healthy |>
        mutate(cluster = sample(c("1", "2"), size = nrow(healthy), replace = TRUE))

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_flowsom(clust_data, perform_metaclustering = FALSE)
        },
        {
            set.seed(2020)
            tof_cluster(clust_data, method = "flowsom", augment = FALSE, perform_metacluster = FALSE)
        }
    )

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_kmeans(clust_data)
        },
        {
            set.seed(2020)
            tof_cluster(clust_data, method = "kmeans", augment = FALSE)
        }
    )

    expect_equal(
        {
            set.seed(2020)
            tof_cluster_phenograph(synth_data)
        },
        {
            set.seed(2020)
            tof_cluster(synth_data, augment = FALSE, method = "phenograph")
        }
    )

    expect_equal(
        tof_cluster_ddpr(
            ddpr_data,
            healthy_tibble = healthy,
            healthy_label_col = cluster
        ),
        tof_cluster(
            ddpr_data,
            method = "ddpr",
            augment = FALSE,
            healthy_tibble = healthy,
            healthy_label_col = cluster
        )
    )
})

# tof_annotate_clusters --------------------------------------------------------

sim_data <-
    dplyr::tibble(
        cd45 = rnorm(n = 1000),
        cd38 = c(rnorm(n = 500), rnorm(n = 500, mean = 2)),
        cd34 = c(rnorm(n = 500), rnorm(n = 500, mean = 4)),
        cd19 = rnorm(n = 1000),
        cluster_id = c(rep("a", 500), rep("b", 500))
    )


test_that("Annotation result is the correct shape", {
    # using named character vector
    annotation_result <-
        sim_data |>
        tof_annotate_clusters(
            cluster_col = cluster_id,
            annotations = c("macrophage" = "a", "dendritic cell" = "b")
        )

    vector_shape <-
        annotation_result |>
        dplyr::count(cluster_id, cluster_id_annotation) |>
        nrow()

    expect_equal(vector_shape, expected = 2L)

    # using tibble
    annotation_data_frame <-
        data.frame(
            cluster_id = c("a", "b"),
            cluster_annotation = c("macrophage", "dendritic cell")
        )

    annotation_result <-
        sim_data |>
        tof_annotate_clusters(
            cluster_col = cluster_id,
            annotations = annotation_data_frame
        )

    tibble_shape <-
        annotation_result |>
        dplyr::count(cluster_id, cluster_annotation) |>
        nrow()

    expect_equal(tibble_shape, expected = 2L)
})
