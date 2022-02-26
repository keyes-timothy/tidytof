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

# cytofkit
read_data_cytofkit <-
  function(file_names) {
    extension_type <- tools::file_ext(file_names[[1]])
    if (extension_type == "csv") {
      data_list <- lapply(X = file_names, read.csv)
      output_data <- rbind(data_list)
    } else {
      output_data <- cytof_exprsMerge(file_names, mergeMethod="all")
      l <- colnames(output_data)
      channel_names  <- sub('.*\\<(.*)\\>.*', '\\1', l)
      metal_names  <- sub('.*\\((.*)\\).*', '\\1', l)
      new_names <- paste(channel_names, metal_names, sep = "_")
      colnames(output_data) <- new_names
    }
    return(output_data)
  }

# immunoCluster
read_data_immunocluster <-
  function(file_names, group) {
    column_names <-
      c(
        "time", "cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "CD235_61",
        "CD45", "empty_1", "cPARP", "empty_2", "pPLCg1_2", "CD19", "CD22",
        "p4EBP1", "tIkaros", "CD79b", "CD20", "CD34", "CD179a", "pSTAT5",
        "CD123", "Ki67", "IgMi", "Kappa_lambda", "pIkaros", "CD10", "empty_3",
        "CD179b", "pAkt", "CD24", "TSLPr", "CD127", "RAG1", "TdT", "Pax5",
        "pSyk", "CD43", "CD38", "CD58", "CD3", "FITC_myeloid", "pS6", "pErk",
        "HLADR", "IgMs", "pCreb", "DNA1", "DNA2", "viability"
      )
    metadata <-
      data.frame(
        file = file_names,
        group = group,
        row.names = file_names,
        stringsAsFactors = FALSE
      )
    output_data <-
      immunoCluster::processFCS(
        files = file_names,
        metadata = metadata,
        transformation = FALSE,
        downsample = NULL,
        downsample_grouping = NULL,
        colsDiscard = NULL,
        newColnames = column_names
      )
    return(output_data)
  }

# spectre
read_data_spectre <-
  function(file_names, file_type) {
    # file_directory <- base::dirname(file_names[[1]])
    new_dir <-
      file.path(
        tempdir(),
        paste0(length(file_names), gsub("\\.", "_", file_type))
      )
    dir.create(path = new_dir, recursive = TRUE)
    for (i in 1:length(file_names)) {
      file_name <- file_names[[i]]
      base_name <- base::basename(file_names[[i]])
      file.copy(from = file_name, to = file.path(new_dir, base_name))
    }
    input_data <-
      Spectre::read.files(
        file.loc = new_dir,
        file.type = file_type,
        do.embed.file.names = TRUE
      )
    data_table <- Spectre::do.merge.files(dat = input_data)
    channel_names <-
      gsub(
        pattern = "^.+Di_|\\(|\\)|Di",
        replacement = "",
        x = colnames(data_table)
      )
    metal_names <-
      gsub(
        pattern = "_.+$|\\(|\\)\\Di",
        replacement = "",
        x = colnames(data_table)
      )
    new_names <- paste(channel_names, metal_names, sep = "_")
    repeat_names <- metal_names == channel_names
    new_names[repeat_names] <- metal_names[repeat_names]
    colnames(data_table) <- new_names
    return(data_table)
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
    return(subsampled_flowset)
  }

# cytofkit (not modular, so must be done by reading in files again)
downsample_cytofkit <-
  function(file_names) {
    output_data <- cytof_exprsMerge(file_names, mergeMethod="ceil", fixedNum = 200)
    l <- colnames(output_data)
    channel_names  <- sub('.*\\<(.*)\\>.*', '\\1', l)
    metal_names  <- sub('.*\\((.*)\\).*', '\\1', l)
    new_names <- paste(channel_names, metal_names, sep = "_")
    colnames(output_data) <- new_names
    return(output_data)
  }

# immunoCluster (not modular, so must be done by reading in files again)
downsample_immunocluster <-
  function(file_names, group, num_cells) {
    column_names <-
      c(
        "time", "cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "CD235_61",
        "CD45", "empty_1", "cPARP", "empty_2", "pPLCg1_2", "CD19", "CD22",
        "p4EBP1", "tIkaros", "CD79b", "CD20", "CD34", "CD179a", "pSTAT5",
        "CD123", "Ki67", "IgMi", "Kappa_lambda", "pIkaros", "CD10", "empty_3",
        "CD179b", "pAkt", "CD24", "TSLPr", "CD127", "RAG1", "TdT", "Pax5",
        "pSyk", "CD43", "CD38", "CD58", "CD3", "FITC_myeloid", "pS6", "pErk",
        "HLADR", "IgMs", "pCreb", "DNA1", "DNA2", "viability"
      )
    metadata <-
      data.frame(
        file = file_names,
        group = group,
        row.names = file_names,
        stringsAsFactors = FALSE
      )
    output_data <-
      immunoCluster::processFCS(
        files = file_names,
        metadata = metadata,
        transformation = FALSE,
        downsample = num_cells,
        downsample_grouping = "group",
        colsDiscard = NULL,
        newColnames = column_names
      )
    return(output_data)
  }

# spectre
downsample_spectre <-
  function(data_table, num_cells) {
    num_groups <- length(unique(data_table[["FileName"]]))
    targets <- rep(x = num_cells, times = num_groups)
    output_data <-
      Spectre::do.subsample(
        dat = data_table,
        targets = targets,
        divide.by = "FileName"
      )
    return(output_data)
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

# cytofkit
preprocess_cytofkit <-
  function(file_names) {
    output_data <-
      cytof_exprsMerge(
        file_names,
        mergeMethod = "all",
        transformMethod = "cytofAsinh"
      )
    l <- colnames(output_data)
    channel_names  <- sub('.*\\<(.*)\\>.*', '\\1', l)
    metal_names  <- sub('.*\\((.*)\\).*', '\\1', l)
    new_names <- paste(channel_names, metal_names, sep = "_")
    colnames(output_data) <- new_names
    return(output_data)
  }

# immunocluster
preprocess_immunocluster <-
  function(file_names, group) {
    column_names <-
      c(
        "time", "cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "CD235_61",
        "CD45", "empty_1", "cPARP", "empty_2", "pPLCg1_2", "CD19", "CD22",
        "p4EBP1", "tIkaros", "CD79b", "CD20", "CD34", "CD179a", "pSTAT5",
        "CD123", "Ki67", "IgMi", "Kappa_lambda", "pIkaros", "CD10", "empty_3",
        "CD179b", "pAkt", "CD24", "TSLPr", "CD127", "RAG1", "TdT", "Pax5",
        "pSyk", "CD43", "CD38", "CD58", "CD3", "FITC_myeloid", "pS6", "pErk",
        "HLADR", "IgMs", "pCreb", "DNA1", "DNA2", "viability"
      )
    metadata <-
      data.frame(
        file = file_names,
        group = group,
        row.names = file_names,
        stringsAsFactors = FALSE
      )
    output_data <-
      immunoCluster::processFCS(
        files = file_names,
        metadata = metadata,
        transformation = TRUE,
        transFun = function(x) asinh(x),
        asinhFactor = 5,
        downsample = NULL,
        downsample_grouping = NULL,
        colsDiscard = NULL,
        newColnames = column_names
      )
    return(output_data)
  }

# spectre
preprocess_spectre <-
  function(data_table) {
    column_indices <- grepl(pattern = "^CD", x = colnames(data_table))
    column_names <- colnames(data_table)[column_indices]
    output_data <-
      Spectre::do.asinh(
        dat = data_table,
        cofactor = 5,
        use.cols = column_names
      )
    return(output_data)
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

#cytofkit
tsne_cytofkit <-
  function(data_frame) {
    dr_colnames <- grepl(pattern = "^CD", x = colnames(data_frame))
    tsne_df <- data_frame[, dr_colnames]
    tsne_embeddings <- cytof_dimReduction(data = tsne_df, method = "tsne")
    final_result <-
      data.frame(TSNE1 = tsne_embeddings[, 1], TSNE2 = tsne_embeddings[, 2])
    return(final_result)
  }

# immunocluster
tsne_immunocluster <-
  function(SCE) {
    dr_indices <- grepl(pattern = "^CD", x = names(SCE))
    dr_names <- names(SCE)[dr_indices]
    embeddings <-
      immunoCluster::performTSNE(
        indata = SCE,
        reducedDim = NULL,
        newDimName = ".tsne",
        useMarkers = dr_names,
        perplexity = 30
      )
    return(embeddings)
  }

# spectre
tsne_spectre <-
  function(data_table) {
    dr_indices <- grepl(pattern = "^CD", x = colnames(data_table))
    dr_names <- colnames(data_table)[dr_indices]
    embeddings <-
      Spectre::run.tsne(
        dat = data_table,
        use.cols = dr_names,
        tsne.x.name = ".tsne1",
        tsne.y.name = ".tsne2",
        verbose = FALSE
      )
    result <- embeddings[, c(".tsne1", ".tsne2")]
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
pca_flowcore <-
  function(flowset) {
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

# cytofkit
pca_cytofkit <-
  function(data_frame) {
    dr_colnames <- grepl(pattern = "^CD", x = colnames(data_frame))
    pca_df <- data_frame[, dr_colnames]
    column_variances <-
      apply(X = pca_df, MARGIN = 2, FUN = function(x) length(unique(x)))
    zv_columns <- as.logical(round(column_variances) == 1)
    pca_df <- pca_df[, !zv_columns]
    pca_embeddings <- cytof_dimReduction(data = pca_df, method = "pca")
    final_result <-
      data.frame(PC1 = pca_embeddings[, 1], PC2 = pca_embeddings[, 2])
    return(final_result)
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
umap_flowcore <-
  function(flowset) {
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

#immunoCluster
umap_immunocluster <-
  function(SCE) {
    dr_indices <- grepl(pattern = "^CD", x = names(SCE))
    dr_names <- names(SCE)[dr_indices]
    embeddings <-
      immunoCluster::performUMAP(
        indata = SCE,
        reducedDim = NULL,
        newDimName = ".umap",
        useMarkers = dr_names
      )
    return(embeddings)
  }

#spectre
umap_spectre <-
  function(data_table) {
    dr_indices <- grepl(pattern = "^CD", x = colnames(data_table))
    dr_names <- colnames(data_table)[dr_indices]
    embeddings <-
      Spectre::run.umap(
        dat = data_table,
        use.cols = dr_names,
        umap.x.name = ".umap1",
        umap.y.name = ".umap2",
        neighbours = 5,
        min_dist = 0.01,
        verbose = FALSE
      )
    result <- embeddings[, c(".umap1", ".umap2")]
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

# cytofkit
flowsom_cytofkit <-
  function(file_names) {
    data_frame <-
      cytof_exprsMerge(
        fcsFiles = file_names,
        transformMethod = "cytofAsinh",
        mergeMethod = "all"
      )
    cluster_colnames <- grepl(pattern = "CD", x = colnames(data_frame))
    data_frame <- data_frame[, cluster_colnames]
    cluster_FlowSOM <-
      cytof_cluster(
        xdata = data_frame,
        method = "FlowSOM",
        FlowSOM_k = 20
      )
    result <- data.frame(.flowsom_metacluster = as.character(cluster_FlowSOM))
    return(result)
  }

# immunocluster
flowsom_immunocluster <-
  function(file_names, group) {
    column_names <-
      c(
        "time", "cell_length", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "CD235_61",
        "CD45", "empty_1", "cPARP", "empty_2", "pPLCg1_2", "CD19", "CD22",
        "p4EBP1", "tIkaros", "CD79b", "CD20", "CD34", "CD179a", "pSTAT5",
        "CD123", "Ki67", "IgMi", "Kappa_lambda", "pIkaros", "CD10", "empty_3",
        "CD179b", "pAkt", "CD24", "TSLPr", "CD127", "RAG1", "TdT", "Pax5",
        "pSyk", "CD43", "CD38", "CD58", "CD3", "FITC_myeloid", "pS6", "pErk",
        "HLADR", "IgMs", "pCreb", "DNA1", "DNA2", "viability"
      )
    metadata <-
      data.frame(
        file = file_names,
        group = group,
        row.names = file_names,
        stringsAsFactors = FALSE
      )
    sce <-
      immunoCluster::processFCS(
        files = file_names,
        metadata = metadata,
        transformation = FALSE,
        downsample = NULL,
        downsample_grouping = NULL,
        colsDiscard = NULL,
        newColnames = column_names
      )
    dr_indices <- grepl(pattern = "^CD", x = names(sce))
    dr_names <- names(sce)[dr_indices]
    result <-
      immunoCluster::runFlowSOM(
        sct = sce,
        k = 20,
        markers = dr_names
      )
    return(result)
  }

flowsom_spectre <-
  function(file_names) {
    new_dir <-
      file.path(
        tempdir(),
        paste0(length(file_names), "_flowsom")
      )
    dir.create(path = new_dir, recursive = TRUE)
    for (i in 1:length(file_names)) {
      file_name <- file_names[[i]]
      base_name <- base::basename(file_names[[i]])
      file.copy(from = file_name, to = file.path(new_dir, base_name))
    }
    input_data <-
      Spectre::read.files(
        file.loc = new_dir,
        file.type = ".fcs",
        do.embed.file.names = TRUE
      )
    data_table <- Spectre::do.merge.files(dat = input_data)
    channel_names <-
      gsub(
        pattern = "^.+Di_|\\(|\\)|Di",
        replacement = "",
        x = colnames(data_table)
      )
    metal_names <-
      gsub(
        pattern = "_.+$|\\(|\\)\\Di",
        replacement = "",
        x = colnames(data_table)
      )
    new_names <- paste(channel_names, metal_names, sep = "_")
    repeat_names <- metal_names == channel_names
    new_names[repeat_names] <- metal_names[repeat_names]
    colnames(data_table) <- new_names
    col_indices <- grepl(pattern = "^CD", x = colnames(data_table))
    col_names <- colnames(data_table)[col_indices]
    result <-
      Spectre::run.flowsom(
        dat = data_table,
        use.cols = col_names,
        xdim = 10,
        ydim = 10,
        max.meta = 20,
        clust.name = ".flowsom_cluster",
        meta.clust.name = ".flowsom_metacluster"
      )
    return(result[, c(".flowsom_cluster", ".flowsom_metacluster")])
  }

# feature extraction -----------------------------------------------------------

# tidytof
extract_tidytof <-
  function(data_frame) {
    data_frame %>%
      tof_extract_features(
        cluster_col = my_cluster,
        group_cols = file_name,
        lineage_cols = starts_with("CD"),
        signaling_cols = starts_with("p", ignore.case = FALSE)
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
        f = list(data_frame$file_name, data_frame$my_cluster),
        sep = "@"
      )
    grouped_signaling_df <-
      split(
        x = signaling_df,
        f = list(data_frame$file_name, data_frame$my_cluster),
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
    lineage_result$my_cluster <-
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
    signaling_result$my_cluster <-
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
    num_clusters <- length(unique(signaling_result$my_cluster))
    num_signaling_channels <-
      length(setdiff(colnames(signaling_result), c("file_name", "my_cluster")))
    num_lineage_channels <-
      length(setdiff(colnames(lineage_result), c("file_name", "my_cluster")))

    lineage_final <-
      matrix(
        data = 0,
        nrow = num_files,
        ncol = num_clusters * num_lineage_channels
      )
    row.names(lineage_final) <- unique(lineage_result$file_name)
    colname_grid <-
      expand.grid(
        setdiff(colnames(lineage_result), c("file_name", "my_cluster")),
        unique(lineage_result$my_cluster)
      )
    colnames(lineage_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(lineage_final)) {
      for (j in 1:ncol(lineage_final)) {
        file <- row.names(lineage_final)[[i]]
        colname <- colnames(lineage_final)[[j]]
        my_cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          lineage_result[
            lineage_result$my_cluster == my_cluster & lineage_result$file_nam == file,
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
        setdiff(colnames(signaling_result), c("file_name", "my_cluster")),
        unique(signaling_result$my_cluster)
      )
    colnames(signaling_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(signaling_final)) {
      for (j in 1:ncol(signaling_final)) {
        file <- row.names(signaling_final)[[i]]
        colname <- colnames(signaling_final)[[j]]
        my_cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          signaling_result[
            signaling_result$my_cluster == my_cluster & signaling_result$file_nam == file,
            channel
          ]
        signaling_final[i, j] <- current_value
      }
    }
    result <- as.data.frame(cbind(lineage_final, signaling_final))
    return(result)
  }

# immunoCluster
# no great way to compute the mean for each cluster within each patient
extract_immunocluster <-
  function(sce) {
    data_matrix <- t(sce@assays@data$scaled)
    data_frame <- as.data.frame(data_matrix)
    data_frame$my_cluster <- sce@metadata$my_cluster
    data_frame$file_name <- sce@metadata$file
    lineage_cols <- grepl(pattern = "^CD", x = colnames(data_frame))
    signaling_cols <- grepl(pattern = "^p", x = colnames(data_frame))
    lineage_df <- data_frame[, lineage_cols]
    signaling_df <- data_frame[, signaling_cols]
    grouped_lineage_df <-
      split(
        x = lineage_df,
        f = list(data_frame$file_name, data_frame$my_cluster),
        sep = "@"
      )
    grouped_signaling_df <-
      split(
        x = signaling_df,
        f = list(data_frame$file_name, data_frame$my_cluster),
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
    lineage_result$my_cluster <-
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
    signaling_result$my_cluster <-
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
    num_clusters <- length(unique(signaling_result$my_cluster))
    num_signaling_channels <-
      length(setdiff(colnames(signaling_result), c("file_name", "my_cluster")))
    num_lineage_channels <-
      length(setdiff(colnames(lineage_result), c("file_name", "my_cluster")))
    lineage_final <-
      matrix(
        data = 0,
        nrow = num_files,
        ncol = num_clusters * num_lineage_channels
      )
    row.names(lineage_final) <- unique(lineage_result$file_name)
    colname_grid <-
      expand.grid(
        setdiff(colnames(lineage_result), c("file_name", "my_cluster")),
        unique(lineage_result$my_cluster)
      )
    colnames(lineage_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(lineage_final)) {
      for (j in 1:ncol(lineage_final)) {
        file <- row.names(lineage_final)[[i]]
        colname <- colnames(lineage_final)[[j]]
        my_cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          lineage_result[
            lineage_result$my_cluster == my_cluster & lineage_result$file_nam == file,
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
        setdiff(colnames(signaling_result), c("file_name", "my_cluster")),
        unique(signaling_result$my_cluster)
      )
    colnames(signaling_final) <-
      paste(colname_grid$Var1, colname_grid$Var2, sep = "@")
    for (i in 1:nrow(signaling_final)) {
      for (j in 1:ncol(signaling_final)) {
        file <- row.names(signaling_final)[[i]]
        colname <- colnames(signaling_final)[[j]]
        my_cluster <-
          substr(x = colname, start = nchar(colname), stop = nchar(colname))
        channel <- sub(pattern = "@.$", replacement = "", x = colname)
        current_value <-
          signaling_result[
            signaling_result$my_cluster == my_cluster & signaling_result$file_nam == file,
            channel
          ]
        signaling_final[i, j] <- current_value
      }
    }
    result <- as.data.frame(cbind(lineage_final, signaling_final))
    return(result)
  }

# spectre
extract_spectre <-
  function(data_table) {
    lineage_cols <- grepl(pattern = "^CD", x = colnames(data_table))
    signaling_cols <- grepl(pattern = "^p", x = colnames(data_table))
    perc_positive <-
      data.table(
        marker_names = colnames(data_table)[signaling_cols],
        threshold = asinh(10 / 5)
      )
    result <-
      Spectre::create.sumtable(
        dat = data_table,
        sample.col = "FileName",
        pop.col = "my_cluster",
        use.cols = colnames(data_table)[lineage_cols],
        perc.pos = perc_positive,
        func = "median"
      )
    return(result)
  }


# memory benchmarking ----------------------------------------------------------
# return the size of a cytof dataset in megabytes as either a tof_tibble


# or as a flowSet
get_memory <- function(file_paths, mode) {
  if (mode == "tidytof") {
    result <-
      file_paths %>%
      tof_read_data() %>%
      lobstr::obj_size() %>%
      as.numeric()

  } else if (mode == "flowcore") {
    result <-
      file_paths %>%
      flowCore::read.flowSet() %>%
      lobstr::obj_size() %>%
      as.numeric()

  } else if (mode == "immunocluster") {
    result <-
      file_paths %>%
      read_data_immunocluster(group = file_paths) %>%
      lobstr::obj_size() %>%
      as.numeric()

  } else if (mode == "spectre") {
    result <-
      file_paths %>%
      read_data_spectre(file_type = ".fcs") %>%
      lobstr::obj_size() %>%
      as.numeric()

  } else if (mode == "cytofkit") {
    result <-
      file_paths %>%
      read_data_cytofkit() %>%
      lobstr::obj_size() %>%
      as.numeric()

  } else {
    stop("not a valid mode")
  }

  return(result / 1e6)
}

# code style benchmarking ------------------------------------------------------


# count the lines of code for a given workflow
get_lines <- function(function_object) {
  if (is.na(function_object)) {
    return(NA)
  } else {
    lines <- capture.output(print(function_object))
    return(length(lines) - 2L)
  }
}

# count the number of assigned variables for a given workflow
get_assignments <- function(function_object) {
  if (is.na(function_object)) {
    return(NA)
  } else {
    lines <- capture.output(print(function_object))
    assignments <- max(sum(str_count(lines, pattern = "<-")), 1L)
    return(assignments)
  }
}
