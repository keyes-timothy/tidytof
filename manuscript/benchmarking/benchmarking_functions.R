# reading data -----------------------------------------------------------------

# tidytof
read_data_tidytof <-
  function(file_names) {
    return(tof_read_data(file_names))
  }

# base (read FCS or CSV file into a tidy data.frame)
read_data_base <-
  function(file_names) {
    extension_type <- tools::file_ext(file_names[[1]])
    if (extension_type == "csv") {
      data_list <- lapply(X = file_names, read.csv)
      output_data <- rbind(data_list)
    } else {
      flowset <- read.flowSet(file_names)
      raw_exprs <- fsApply(flowset, exprs)
      output_data <- as.data.frame(raw_exprs)
      channel_names <- as.character(flowset[[1]]@parameters[["desc"]])
      metal_names <- as.character(flowset[[1]]@parameters[["name"]])
      metals_to_change <- grepl(pattern = "Di", x = metal_names)
      metal_names[metals_to_change] <-
        gsub(
          pattern = "\\(|\\)|(Di)",
          replacement = "",
          x = metal_names[metals_to_change]
        )
      new_names <- paste(channel_names, metal_names, sep = "_")
      new_names[!metals_to_change] <- channel_names[!metals_to_change]
      colnames(output_data) <- new_names
    }
    return(output_data)
  }

# read a flowset with names cleaned up in the same way that tidytof does it
read_data_flowcore <- function(file_names) {
  flowset <- read.flowSet(file_names)
  channel_names <- as.character(flowset[[1]]@parameters[["desc"]])
  metal_names <- as.character(flowset[[1]]@parameters[["name"]])
  metals_to_change <- grepl(pattern = "Di", x = metal_names)
  metal_names[metals_to_change] <-
    gsub(
      pattern = "\\(|\\)|(Di)",
      replacement = "",
      x = metal_names[metals_to_change]
    )
  new_names <- paste(channel_names, metal_names, sep = "_")
  new_names[!metals_to_change] <- channel_names[!metals_to_change]
  num_files <- length(flowset)
  flowframe_list <- list()
  for (i in 1:num_files) {
    flowframe <- flowset[[i]]
    flowframe@parameters[["desc"]] <- new_names
    flowframe@parameters[["name"]] <- new_names
    flowframe_list[[i]] <- flowframe
  }
  result <- flowSet(flowframe_list)
  return(result)
}


# downsampling -----------------------------------------------------------------

# tidytof
downsample_tidytof <-
  function(data_frame) {
    result <-
      data_frame %>%
      tof_downsample_constant(group_cols = file_name, num_cells = 200)
    return(result)
  }

# base R
downsample_base <-
  function(data_frame) {
    file_names <- unique(data_frame$file_name)
    num_file_names <- length(file_names)
    final_subset <- list()
    for (i in 1:num_file_names) {
      file_name <- file_names[[i]]
      file_indices <- which(data_frame$file_name == file_name)
      subset_indices <- sample(x = file_indices, size = 200L)
      final_subset[[file_name]] <- subset_indices
    }
    final_subset <- sort(as.numeric(c(final_subset, recursive = TRUE)))
    result <- data_frame[final_subset, ]
    return(result)
  }

# flowcore
downsample_flowcore <-
  function(flowset) {
    subset_flowframe <- function(flowframe) {
      num_cells <- nrow(flowframe)
      sample_indices <- sort(sample(x = 1:num_cells, size = 200L))
      subsampled_flowframe <- flowframe[sample_indices]
      return(subsampled_flowframe)
    }
    subsampled_flowset <- flowCore::fsApply(x = flowset, FUN = subset_flowframe)
  }

# preprocessing ----------------------------------------------------------------

# tidytof
preprocess_tidytof <-
  function(data_frame) {
    return(tof_preprocess(data_frame, undo_noise = FALSE))
  }

# base R
preprocess_base <-
  function(data_frame) {
    num_cols <- ncol(data_frame)
    for (i in 1:num_cols) {
      col_values <- data_frame[[i]]
      if (is.numeric(col_values)) {
        new_values <- asinh(col_values / 5)
        data_frame[[i]] <- new_values
      }
    }
    return(data_frame)
  }

# flowCore
preprocess_flowcore <-
  function(flowset) {
    numeric_columns <-
      as.logical(apply(X = exprs(flowset[[1]]), MARGIN = 2, FUN = is.numeric))
    asinh_transform <- arcsinhTransform(a = 0, b = 0.2, c = 0)
    translist <-
      transformList(colnames(flowset)[numeric_columns], asinh_transform)
    preprocessed_flowset <- transform(flowset, translist)
    return(preprocessed_flowset)
  }

# dimensionality reduction -----------------------------------------------------

### tsne

# tidytof
tsne_tidytof <-
  function(data_frame) {
    return(tof_reduce_tsne(data_frame, tsne_cols = starts_with("CD")))
  }

# base R
tsne_base <-
  function(data_frame) {
    dr_colnames <- grepl(pattern = "^CD", x = colnames(data_frame))
    tsne_df <- data_frame[, dr_colnames]
    tsne_result <- Rtsne::Rtsne(X = as.matrix(tsne_df))
    tsne_embeddings <- tsne_result$Y
    final_result <-
      data.frame(tsne_1 = tsne_embeddings[, 1], tsne_2 = tsne_embeddings[, 2])
    return(final_result)
  }

# flowcore
tsne_flowcore <-
  function(flowset) {
    channel_names <- flowset[[1]]@parameters@data$desc
    dr_columns <- grepl(pattern = "^CD", x = channel_names)
    my_exprs <- fsApply(x = flowset, FUN = exprs)
    dr_exprs <- my_exprs[, dr_columns]
    tsne_result <- Rtsne::Rtsne(X = as.matrix(dr_exprs))
    tsne_embeddings <- tsne_result$Y
    colnames(tsne_embeddings) <- paste0(".tsne", 1:2)
    flowframe_num_cells <- as.numeric(fsApply(x = flowset, FUN = nrow))
    starting_index <- 1L
    flowframe_list <- list()
    for (i in 1:length(flowframe_num_cells)) {
      flowframe <- flowset[[i]]
      num_cells <- flowframe_num_cells[[i]]
      ending_index <- starting_index + num_cells - 1
      flowframe_tsne <- tsne_embeddings[starting_index:ending_index, ]
      new_flowframe <- fr_append_cols(flowframe, flowframe_tsne)
      flowframe_list[[i]] <- new_flowframe
      starting_index <- ending_index + 1
    }
    result <- flowSet(flowframe_list)
    return(result)
  }

### pca

# tidytof
pca_tidytof <-
  function(data_frame) {
    return(tof_reduce_pca(data_frame, pca_cols = starts_with("CD")))
  }

# base R
pca_base <-
  function(data_frame) {
    dr_colnames <- grepl(pattern = "^CD", x = colnames(data_frame))
    pca_df <- data_frame[, dr_colnames]
    column_variances <-
      apply(X = pca_df, MARGIN = 2, FUN = function(x) length(unique(x)))
    zv_columns <- as.logical(round(column_variances) == 1)
    pca_df <- pca_df[, !zv_columns]
    pca_result <- prcomp(x = as.matrix(pca_df), center = TRUE, scale. = TRUE)
    pca_embeddings <- pca_result$x
    final_result <-
      data.frame(PC1 = pca_embeddings[, 1], PC2 = pca_embeddings[, 2])
    return(final_result)
  }

# flowcore
pca_flowcore <- function(flowset) {
  channel_names <- flowset[[1]]@parameters@data$desc
  dr_columns <- grepl(pattern = "^CD", x = channel_names)
  my_exprs <- fsApply(x = flowset, FUN = exprs)
  dr_exprs <- my_exprs[, dr_columns]
  column_variances <-
    apply(X = dr_exprs, MARGIN = 2, FUN = function(x) length(unique(x)))
  zv_columns <- as.logical(round(column_variances) == 1)
  dr_exprs <- dr_exprs[, !zv_columns]
  pca_result <- prcomp(x = as.matrix(dr_exprs), center = TRUE, scale. = TRUE)
  pca_embeddings <- pca_result$x[, 1:2]
  colnames(pca_embeddings) <- paste0("PC", 1:2)
  flowframe_num_cells <- as.numeric(fsApply(x = flowset, FUN = nrow))
  starting_index <- 1L
  flowframe_list <- list()
  for (i in 1:length(flowframe_num_cells)) {
    flowframe <- flowset[[i]]
    num_cells <- flowframe_num_cells[[i]]
    ending_index <- starting_index + num_cells - 1
    flowframe_pca <- pca_embeddings[starting_index:ending_index, ]
    new_flowframe <- fr_append_cols(flowframe, flowframe_pca)
    flowframe_list[[i]] <- new_flowframe
    starting_index <- ending_index + 1
  }
  result <- flowSet(flowframe_list)
  return(result)
}

### umap

# tidytof
umap_tidytof <-
  function(data_frame) {
    return(tof_reduce_umap(data_frame, umap_cols = starts_with("CD")))
  }

# base R
umap_base <-
  function(data_frame) {
    dr_colnames <- grepl(pattern = "^CD", x = colnames(data_frame))
    umap_df <- data_frame[, dr_colnames]
    column_variances <-
      apply(X = umap_df, MARGIN = 2, FUN = function(x) length(unique(x)))
    zv_columns <- as.logical(round(column_variances) == 1)
    umap_df <- umap_df[, !zv_columns]
    umap_result <-
      uwot::umap(X = umap_df, n_neighbors = 5, scale = TRUE)
    final_result <-
      data.frame(UMAP1 = umap_result[, 1], UMAP2 = umap_result[, 2])
    return(final_result)
  }

# flowcore
umap_flowcore <- function(flowset) {
  channel_names <- flowset[[1]]@parameters@data$desc
  dr_columns <- grepl(pattern = "^CD", x = channel_names)
  my_exprs <- fsApply(x = flowset, FUN = exprs)
  dr_exprs <- my_exprs[, dr_columns]
  column_variances <-
    apply(X = dr_exprs, MARGIN = 2, FUN = function(x) length(unique(x)))
  zv_columns <- as.logical(round(column_variances) == 1)
  dr_exprs <- dr_exprs[, !zv_columns]
  umap_result <-
    uwot::umap(X = dr_exprs, n_neighbors = 5, scale = TRUE)
  colnames(umap_result) <- paste0("UMAP", 1:2)
  flowframe_num_cells <- as.numeric(fsApply(x = flowset, FUN = nrow))
  starting_index <- 1L
  flowframe_list <- list()
  for (i in 1:length(flowframe_num_cells)) {
    flowframe <- flowset[[i]]
    num_cells <- flowframe_num_cells[[i]]
    ending_index <- starting_index + num_cells - 1
    flowframe_umap <- umap_result[starting_index:ending_index, ]
    new_flowframe <- fr_append_cols(flowframe, flowframe_umap)
    flowframe_list[[i]] <- new_flowframe
    starting_index <- ending_index + 1
  }
  result <- flowSet(flowframe_list)
  return(result)
}


# clustering -------------------------------------------------------------------

# tidytof
flowsom_tidytof <-
  function(file_names) {
    clusters <-
      file_names %>%
      tof_read_data() %>%
      tof_preprocess() %>%
      tof_cluster_flowsom(
        cluster_cols = starts_with("CD", ignore.case = FALSE),
        perform_metaclustering = TRUE
      )
    return(clusters)
  }

# base R (+ flowcore)
flowsom_base <-
  function(file_names) {
    flowset <-
      flowCore::read.flowSet(
        files = file_names,
        transformation = FALSE,
        truncate_max_range = FALSE
      )
    raw_exprs <- flowCore::fsApply(flowset, flowCore::exprs, simplify = FALSE)
    asinh_exprs <- lapply(X = raw_exprs, FUN = function(x) asinh(x / 5))
    for (i in 1:length(asinh_exprs)) {
      flowCore::exprs(flowset[[i]]) <- asinh_exprs[[i]]
    }
    channel_names <- as.character(flowset[[1]]@parameters[["desc"]])
    cluster_colnames <-
      grepl(pattern = "^CD", x = channel_names, ignore.case = FALSE)
    flowsom_object <- FlowSOM::ReadInput(input = flowset)
    flowsom_som <-
      FlowSOM::BuildSOM(fsom = flowsom_object, colsToUse = cluster_colnames)
    cluster_labels <- flowsom_som$map$mapping[, 1]
    metaclusters <-
      FlowSOM::MetaClustering(
        data = flowsom_som$map$codes,
        method = "metaClustering_consensus",
        max = 20
      )
    metacluster_labels <- metaclusters[flowsom_som$map$mapping[, 1]]
    clusters <- data.frame(.flowsom_metacluster = metacluster_labels)
  }

# feature extraction -----------------------------------------------------------

# tidytof
extract_tidytof <-
  function(data_frame) {
    data_frame %>%
      tof_extract_features(
        cluster_col = cluster,
        group_cols = file_name,
        lineage_cols = starts_with("CD"),
        signaling_cols = starts_with("p")
      )
  }

# base R
extract_base <-
  function(data_frame) {
    lineage_cols <- grepl(pattern = "^CD", x = colnames(data_frame))
    signaling_cols <- grepl(pattern = "^p", x = colnames(data_frame))
    lineage_df <- data_frame[, lineage_cols]
    signaling_df <- data_frame[, signaling_cols]
    grouped_lineage_df <-
      split(
        x = lineage_df,
        f = list(data_frame$file_name, data_frame$cluster),
        sep = "@"
      )
    grouped_signaling_df <-
      split(
        x = signaling_df,
        f = list(data_frame$file_name, data_frame$cluster),
        sep = "@"
      )
    lineage_medians <-
      sapply(
        X = grouped_lineage_df,
        FUN = function(x) apply(X = x, MARGIN = 2, FUN = median)
      )
    signaling_thresh <-
      sapply(
        X = grouped_signaling_df,
        FUN =
          function(x) {
            apply(X = x, MARGIN = 2, FUN = function(y) mean(y > asinh(10 / 5)))
          }
      )
    lineage_result <- as.data.frame(t(lineage_medians))
    signaling_result <- as.data.frame(t(signaling_thresh))
    lineage_result$cluster <-
      substr(
        x = row.names(lineage_result),
        start = nchar(row.names(lineage_result)),
        stop = nchar(row.names(lineage_result))
      )
    lineage_result$file_name <-
      sub(
        pattern = "@.$",
        replacement = "",
        x = row.names(lineage_result)
      )
    signaling_result$cluster <-
      substr(
        x = row.names(signaling_result),
        start = nchar(row.names(signaling_result)),
        stop = nchar(row.names(signaling_result))
      )
    signaling_result$file_name <-
      sub(
        pattern = "@.$",
        replacement = "",
        x = row.names(signaling_result)
      )
    num_files <- length(unique(signaling_result$file_name))
    num_clusters <- length(unique(signaling_result$cluster))
    num_signaling_channels <-
      length(setdiff(colnames(signaling_result), c("file_name", "cluster")))
    num_lineage_channels <-
      length(setdiff(colnames(lineage_result), c("file_name", "cluster")))

    lineage_final <-
      matrix(
        data = 0,
        nrow = num_files,
        ncol = num_clusters * num_lineage_channels
      )
    row.names(lineage_final) <- unique(lineage_result$file_name)
    colname_grid <-
      expand.grid(
        setdiff(colnames(lineage_result), c("file_name", "cluster")),
        unique(lineage_result$cluster)
      )
    colnames(lineage_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(lineage_final)) {
      for (j in 1:ncol(lineage_final)) {
        file <- row.names(lineage_final)[[i]]
        colname <- colnames(lineage_final)[[j]]
        cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          lineage_result[
            lineage_result$cluster == cluster & lineage_result$file_nam == file,
            channel
          ]
        lineage_final[i, j] <- current_value
      }
    }
    signaling_final <-
      matrix(
        data = 0,
        nrow = num_files,
        ncol = num_clusters * num_signaling_channels
      )
    row.names(signaling_final) <- unique(lineage_result$file_name)
    colname_grid <-
      expand.grid(
        setdiff(colnames(signaling_result), c("file_name", "cluster")),
        unique(signaling_result$cluster)
      )
    colnames(signaling_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(signaling_final)) {
      for (j in 1:ncol(signaling_final)) {
        file <- row.names(signaling_final)[[i]]
        colname <- colnames(signaling_final)[[j]]
        cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          signaling_result[
            signaling_result$cluster == cluster & signaling_result$file_nam == file,
            channel
          ]
        signaling_final[i, j] <- current_value
      }
    }
    result <- as.data.frame(cbind(lineage_final, signaling_final))
    return(result)
  }
