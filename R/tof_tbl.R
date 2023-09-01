
# new_tof_tibble ---------------------------------------------------------------

#' Constructor for a tof_tibble.
#'
#' @param x A data.frame or tibble containing single-cell mass cytometry data
#' such that rows are cells and columns are CyTOF measurements.
#'
#' @param panel A data.frame or tibble containing information about the panel
#' for the mass cytometry data in x.
#'
#' @return A `tof_tbl`, an tibble extension that tracks a few other attributes
#' that are useful for CyTOF data analysis.
#'
#' @family tof_tbl utilities
#'
new_tof_tibble <- function(x = dplyr::tibble(), panel = dplyr::tibble()) {

  stopifnot(inherits(x, "tbl_df"))
  stopifnot(inherits(panel, "tbl_df"))

  if("grouped_df" %in% class(x)) {
    subclasses <- c("grouped_tof_tbl", "grouped_df", "tof_tbl")
  } else {
    subclasses <- "tof_tbl"
  }

  tibble::new_tibble(
    x,
    panel = panel,
    nrow = nrow(x),
    class = subclasses
  )
}

# tof_get_panel ----------------------------------------------------------------

#' Get panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @return A tibble containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in `tof_tibble`.
#'
#' @export
#'
#' @family tof_tbl utilities
#'
#' @examples
#' input_file <- dir(tidytof_example_data("aml"), full.names = TRUE)[[1]]
#' tof_tibble <- tof_read_data(input_file)
#' tof_get_panel(tof_tibble)
#'
#'
tof_get_panel <- function(tof_tibble) {
  panel <-
    tof_tibble %>%
    attr(which = "panel")

  return(panel)
}


# tof_set_panel ----------------------------------------------------------------

#' Set panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @param panel A tibble containing two columns (`metals` and `antigens`) representing
#' the information about a panel
#'
#' @return A `tof_tibble` containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in the input `tof_tibble`. Two columns are required: "metals" and "antigens".
#'
#' @family tof_tbl utilities
#'
#' @export
#'
#' @examples
#' # get current panel from an .fcs file
#' input_file <- dir(tidytof_example_data("aml"), full.names = TRUE)[[1]]
#' tof_tibble <- tof_read_data(input_file)
#' current_panel <- tof_get_panel(tof_tibble)
#'
#' # create a new panel (remove empty channels)
#' new_panel <- dplyr::filter(current_panel, antigens != "empty")
#' tof_set_panel(tof_tibble = tof_tibble, panel = new_panel)
#'
#'
tof_set_panel <- function(tof_tibble, panel) {
  attr(tof_tibble, which = "panel") <- panel
  return(tof_tibble)
}


# tof_tbl methods --------------------------------------------------------------

## tidyr methods

#' @export
nest.tof_tbl <- function(.data, ..., .names_sep = NULL) {#, .key = deprecated()) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
unnest.tof_tbl <- function(data, ...) {
  start_panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = start_panel))
}

#' @export
pivot_longer.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
pivot_wider.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
#'
#' @importFrom dplyr group_by_drop_default
group_by.tof_tbl <-
  function(.data, ..., .add = FALSE, .drop = dplyr::group_by_drop_default(.data)) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }


## dplyr methods

#' @export
mutate.tof_tbl <- function(.data, ...) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
slice_sample.tof_tbl <-
  function(.data, ..., n, prop, weight_by = NULL, replace = FALSE) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }

# grouped_tof_tbl methods ------------------------------------------------------

## tidyr methods

#' @export
nest.grouped_tof_tbl <- nest.tof_tbl

#' @export
ungroup.grouped_tof_tbl <- function(x, ...) {
  panel <- tof_get_panel(x)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

# for interoperability with flowCore -------------------------------------------

#' Convert an object into a tof_tbl
#'
#' @param flow_data A FlowSet
#'
#' @param sep A string to use to separate the antigen name and its associated
#' metal in the column names of the output tibble. Defaults to "|".
#'
#' @export
#'
#' @importFrom flowCore fsApply
#' @importFrom flowCore exprs
#'
#' @importFrom dplyr as_tibble
#'
#' @return a `tof_tbl`
#'
#'
as_tof_tbl.flowSet <- function(flow_data, sep = "|") {
  # check if flowset is empty
  if (length(flow_data) < 1) {
    stop("This flowSet is empty.")
  }
  panel_info <-
    flow_data[[1]] %>%
    tof_find_panel_info()

  flowset_exprs <-
    flow_data %>%
    flowCore::fsApply(FUN = flowCore::exprs) %>%
    dplyr::as_tibble()

  col_names <-
    base::paste(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  colnames(flowset_exprs) <- col_names

  result <- new_tof_tibble(x = flowset_exprs, panel = panel_info)

  return(result)
}

#' @export
#'
#' @importFrom dplyr as_tibble
#' @importFrom flowCore exprs
#'
as_tof_tbl.flowFrame <- function(flow_data, sep = "|") {
  panel_info <-
    flow_data %>%
    tof_find_panel_info()

  col_names <-
    #stringr::str_c(panel_info$antigens, panel_info$metals, sep = sep)
    base::paste(panel_info$antigens, panel_info$metals, sep = sep)

  # prevent repeating names twice when antigen and metal are identical
  repeat_indices <-
    which(panel_info$metals == panel_info$antigens)
  col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

  flowframe_exprs <-
    setNames(
      object = dplyr::as_tibble(flowCore::exprs(flow_data)),
      nm = col_names
    )

  result <-
    new_tof_tibble(
      x = flowframe_exprs,
      panel = panel_info
    )

  return(result)

}

#' Coerce flowFrames or flowSets into tof_tbl's.
#'
#' @param flow_data A flowFrame or flowSet
#'
#' @param sep A string indicating which symbol should be used to separate
#' antigen names and metal names in the columns of the output tof_tbl.
#'
#' @export
#'
#' @return A tof_tbl.
#'
#' @examples
#' input_file <- dir(tidytof_example_data("aml"), full.names = TRUE)[[1]]
#'
#' input_flowframe <- flowCore::read.FCS(input_file)
#'
#' tof_tibble <- as_tof_tbl(input_flowframe)
#'
as_tof_tbl <- function(flow_data, sep = "|") {
  UseMethod("as_tof_tbl")
}

# for interoperability with Bioconductor ---------------------------------------

#' Coerce an object into a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @param x An object.
#'
#' @param ... Method-specific arguments
#'
#' @export
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @examples
#' NULL
#'
as_SingleCellExperiment <- function(x, ...) {
  UseMethod("as_SingleCellExperiment")
}



#' Coerce a tof_tbl into a \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'
#' @param x A tof_tbl
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#' If nothing is specified, the default is all numeric columns.
#'
#' @param reduced_dimensions_cols Unquoted column names representing columns that contain
#' dimensionality reduction embeddings, such as tSNE or UMAP embeddings.
#' Supports tidyselect helpers.
#'
#' @param metadata_cols Unquoted column names representing columns that contain
#' metadata about the samples from which each cell was collected. If nothing
#' is specified, the default is all non-numeric columns.
#'
#' @param split_reduced_dimensions A boolean value indicating whether the
#' dimensionality results in x should be split into separate slots in the resulting
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}. If FALSE (the default),
#' the split will not be performed and the
#' \code{\link[SingleCellExperiment]{reducedDims}} slot in the result will have
#' a single entry ("tidytof_reduced_dimensions"). If TRUE, the split will be
#' performed and the \code{\link[SingleCellExperiment]{reducedDims}} slot in
#' the result will have 1-4 entries depending on which dimensionality reduction
#' results are present in x ("tidytof_pca", "tidytof_tsne", "tidytof_umap",
#' and "tidytof_reduced_dimensions"). Note that "tidytof_reduced_dimensions" will
#' include all dimensionality reduction results that are not named according to
#' tidytof's pca, umap, and tsne conventions.
#'
#' @param ... Unused.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#'
#'
#' @rdname as_SingleCellExperiment
#'
#' @export
#'
#'
#' @examples
#' NULL
#'
as_SingleCellExperiment.tof_tbl <-
  function(
    x,
    channel_cols = where(tof_is_numeric),
    reduced_dimensions_cols,
    metadata_cols = where(\(.x) !tof_is_numeric(.x)),
    split_reduced_dimensions = FALSE,
    ...
  ) {

    # check to see if SingleCellExperiment is installed
    rlang::check_installed(pkg = "SingleCellExperiment")

    if (!requireNamespace(package = "SingleCellExperiment")) {
      stop("as_SingleCellExperiment requires the SingleCellExperiment package to be installed from Bioconductor.")
    }

    # extract column names
    channel_colnames <-
      x |>
      dplyr::select({{ channel_cols }}) |>
      colnames()

    reduced_dimensions_colnames <-
      x |>
      dplyr::select({{ reduced_dimensions_cols }}) |>
      colnames()

    metadata_colnames <-
      x |>
      dplyr::select( {{ metadata_cols }}) |>
      colnames()

    # extract embedding (reduced dimension) data
    if (length(reduced_dimensions_colnames) > 1) {
      cytometry_reduced_dimensions <-
        x |>
        dplyr::select({{ reduced_dimensions_cols }}) |>
        as.matrix()

    } else {
      reduced_dimensions_colnames <-
        x |>
        dplyr::select(dplyr::matches("^.pc\\d+|^.tsne\\d+|^.umap\\d+")) |>
        colnames()

      if (length(reduced_dimensions_colnames) == 0) {
        cytometry_reduced_dimensions <- NULL
      } else {
        cytometry_reduced_dimensions <-
          x |>
          dplyr::select(dplyr::any_of(reduced_dimensions_colnames)) |>
          as.matrix()
      }
    }

    # extract marker data
    ## remove any dimensionality reduction columns from the cytometry assay
    ## in the event that they were accidentally included.
    channel_colnames <- setdiff(channel_colnames, reduced_dimensions_colnames)
    cytometry_data <-
      x |>
      dplyr::select({{ channel_cols }}) |>
      dplyr::select(dplyr::any_of(channel_colnames)) |>
      as.matrix() |>
      t()
    row.names(cytometry_data) <- channel_colnames

    # extract metadata
    cytometry_metadata <-
      x |>
      dplyr::select({{ metadata_cols }}) |>
      as.data.frame()

    # assemble SCE object
    tof_assays <- list(cytometry = cytometry_data)

    if (is.null(cytometry_reduced_dimensions)) {
      reduced_dims <- list()
    } else {

      reduced_dims <-
        list(tidytof_reduced_dimensions = cytometry_reduced_dimensions)
    }

    row_data <- data.frame(marker_name = channel_colnames)

    col_data <- cytometry_metadata
    row.names(col_data) <- paste0("cell_", 1:nrow(col_data))

    result <-
      SingleCellExperiment::SingleCellExperiment(
       assays = tof_assays,
       rowData = row_data,
       colData = col_data,
       reducedDims = reduced_dims
      )

    if (split_reduced_dimensions) {
      result <- tof_split_tidytof_reduced_dimensions(result)
    }

    return(result)
  }


#' Coerce an object into a \code{\link[SeuratObject]{SeuratObject}}
#'
#' @param x An object
#'
#' @param ... Method-specific arguments
#'
#' @export
#'
#' @return A \code{\link[SeuratObject]{SeuratObject}}
#'
#' @examples
#' NULL
#'
as_seurat <- function(x, ...) {
  UseMethod("as_seurat")
}



#' Coerce a tof_tbl into a \code{\link[SeuratObject]{SeuratObject}}
#'
#' @param x A tof_tbl
#'
#' @param channel_cols Unquoted column names representing columns that contain
#' single-cell protein measurements. Supports tidyselect helpers.
#' If nothing is specified, the default is all numeric columns.
#'
#' @param reduced_dimensions_cols Unquoted column names representing columns that contain
#' dimensionality reduction embeddings, such as tSNE or UMAP embeddings.
#' Supports tidyselect helpers.
#'
#' @param metadata_cols Unquoted column names representing columns that contain
#' metadata about the samples from which each cell was collected. If nothing
#' is specified, the default is all non-numeric columns.
#'
#' @param split_reduced_dimensions A boolean value indicating whether the
#' dimensionality results in x should be split into separate slots in the resulting
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}. If FALSE (the default),
#' the split will not be performed and the
#' \code{\link[SingleCellExperiment]{reducedDims}} slot in the result will have
#' a single entry ("tidytof_reduced_dimensions"). If TRUE, the split will be
#' performed and the \code{\link[SingleCellExperiment]{reducedDims}} slot in
#' the result will have 1-4 entries depending on which dimensionality reduction
#' results are present in x ("tidytof_pca", "tidytof_tsne", "tidytof_umap",
#' and "tidytof_reduced_dimensions"). Note that "tidytof_reduced_dimensions" will
#' include all dimensionality reduction results that are not named according to
#' tidytof's pca, umap, and tsne conventions.
#'
#' @param ... Unused.
#'
#' @return A \code{\link[SeuratObject]{SeuratObject}}.
#'
#' @rdname as_seurat
#'
#' @export
#'
#' @examples
#' NULL
#'
as_seurat.tof_tbl <-
  function(
    x,
    channel_cols = where(tof_is_numeric),
    reduced_dimensions_cols,
    metadata_cols = where(\(.x) !tof_is_numeric(.x)),
    split_reduced_dimensions = FALSE,
    ...
  ) {

    # check to see if SingleCellExperiment is installed
    rlang::check_installed(pkg = "SingleCellExperiment")

    if (!requireNamespace(package = "SingleCellExperiment")) {
      stop("as_seurat requires the SingleCellExperiment package to be installed from Bioconductor.")
    }

    # check to see if Seurat is installed
    rlang::check_installed(pkg = "Seurat")

    if (!requireNamespace(package = "Seurat")) {
      stop("as_seurat requires the Seurat package to be installed from CRAN.")
    }

    # check to see if SeuratObject is installed
    rlang::check_installed(pkg = "SeuratObject")

    if (!requireNamespace(package = "SeuratObject")) {
      stop("as_seurat requires the SeuratObject package to be installed from CRAN.")
    }


    if (missing(reduced_dimensions_cols)) {
      reduced_dimensions_cols <-
        x |>
        dplyr::select(dplyr::matches("^.pc\\d+|^.tsne\\d+|^.umap\\d+")) |>
        colnames()
    }

      sce <-
        x |>
        as_SingleCellExperiment(
          channel_cols = {{ channel_cols }},
          reduced_dimensions_cols = {{ reduced_dimensions_cols }},
          metadata_cols = {{ metadata_cols }}
        )

      if (split_reduced_dimensions) {
        sce <-
          sce |>
          tof_split_tidytof_reduced_dimensions()
      }


    #}

    suppressWarnings(suppressMessages(
      result <-
        sce |>
        Seurat::as.Seurat(
          counts = NULL,
          data = "cytometry",
          project = "cytometry"
        )
    ))
    # refine the default coercion a bit to make the Seurat object more intuitive
    # for the user

    suppressWarnings(suppressMessages(
      result <- SeuratObject::RenameAssays(result, originalexp = "cytometry")
    ))

    return(result)
  }


#' Split the dimensionality reduction data that tidytof combines during \code{\link[SingleCellExperiment]{SingleCellExperiment}} conversion
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with an
#' entry named "tidytof_reduced_dimensions" in its \code{\link[SingleCellExperiment]{reducedDims}} slot.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} with separate entries
#' named "tidytof_pca", "tidytof_umap", and "tidytof_tsne" in its
#' \code{\link[SingleCellExperiment]{reducedDims}} slots (one for each of the
#' dimensionality reduction methods for which tidytof has native support).
#'
#'
#' @importFrom rlang check_installed
#' @importFrom rlang is_installed
#'
#' @importFrom purrr discard
#'
#' @examples
#' NULL
tof_split_tidytof_reduced_dimensions <- function(sce) {

  # check to see if SingleCellExperiment is installed
  rlang::check_installed(pkg = "SingleCellExperiment")

  if (!requireNamespace(package = "SingleCellExperiment")) {
    stop("as_seurat requires the SingleCellExperiment package to be installed from Bioconductor.")
  }

  sce_reduced_dimensions <-
    sce |>
    SingleCellExperiment::reducedDim("tidytof_reduced_dimensions")

  sce_pca <-
    sce_reduced_dimensions[,
                           grepl(
                             pattern = "^\\.pc\\d+",
                             x = colnames(sce_reduced_dimensions)
                           )
    ]

  sce_umap <-
    sce_reduced_dimensions[,
                           grepl(
                             pattern = "^\\.umap\\d+",
                             x = colnames(sce_reduced_dimensions)
                           )
    ]

  sce_tsne <-
    sce_reduced_dimensions[,
                           grepl(
                             pattern = "^\\.tsne\\d+",
                             x = colnames(sce_reduced_dimensions)
                           )
    ]

  sce_leftover <-
    sce_reduced_dimensions[,
                           !grepl(
                             pattern = "^\\.pc\\d+|^\\.umap\\d+|^\\.tsne\\d+",
                             x = colnames(sce_reduced_dimensions)
                           )
    ]

  result_list <-
    list(
      tidytof_pca = sce_pca,
      tidytof_tsne = sce_tsne,
      tidytof_umap = sce_umap,
      tidytof_reduced_dimensions = sce_leftover
    ) |>
    purrr::discard(.p = \(.x) ncol(.x) == 0)

  SingleCellExperiment::reducedDim(sce, "tidytof_reduced_dimensions") <- NULL

  for (i in 1:length(result_list)) {
    name <- names(result_list)[[i]]
    result <- result_list[[i]]

    SingleCellExperiment::reducedDim(sce, name) <- result
  }
  return(sce)

}




#' Coerce an object into a \code{\link[flowCore]{flowFrame}}
#'
#' @param x An object.
#'
#' @param ... Method-specific arguments
#'
#' @export
#'
#' @return A \code{\link[flowCore]{flowFrame}}
#'
#' @examples
#' NULL
#'
as_flowFrame <- function(x, ...) {
  UseMethod("as_flowFrame")
}


#' Coerce a tof_tbl into a \code{\link[flowCore]{flowFrame}}
#'
#' @param x A tof_tbl.
#'
#' @param ... Unused.
#'
#' @return A \code{\link[flowCore]{flowFrame}}. Note that all non-numeric
#' columns in `x` will be removed.
#'
#' @rdname as_flowFrame
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#'
#' @importFrom flowCore flowFrame
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' NULL
as_flowFrame.tof_tbl <- function(x, ...) {
  tof_tibble <-
    x |>
    dplyr::select(where(tof_is_numeric))

  maxes_and_mins <-
    tof_tibble |>
    dplyr::summarize(
      dplyr::across(
        dplyr::everything(),
        .fns =
          list(max = ~ max(.x, na.rm = TRUE), min = ~ min(.x, na.rm = TRUE)),
        # use the many underscores because it's unlikely this will come up
        # in column names on their own
        .names = "{.col}_____{.fn}"
      )
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = c("antigen", "value_type"),
      values_to = "value",
      names_sep = "_____"
    )  |>
    tidyr::pivot_wider(
      names_from = "value_type",
      values_from = "value"
    )

  # extract the names of all columns to used in the flowFrame
  data_cols <- maxes_and_mins$antigen

  # make the AnnotatedDataFrame flowCore needs
  parameters <-
    make_flowcore_annotated_data_frame(maxes_and_mins = maxes_and_mins)

  # assemble flowFrame
  result <-
    tof_tibble |>
    dplyr::rename_with(
      .fn = stringr::str_replace,
      pattern = "\\|",
      replacement = "_"
    ) |>
    as.matrix() |>
    flowCore::flowFrame(
      parameters = parameters
    )

  return(result)
}





#' Coerce an object into a \code{\link[flowCore]{flowSet}}
#'
#' @param x An object.
#'
#' @param ... Method-specific arguments
#'
#' @export
#'
#' @return A \code{\link[flowCore]{flowSet}}
#'
#' @examples
#' NULL
#'
as_flowSet <- function(x, ...) {
  UseMethod("as_flowSet")
}


#' Coerce a tof_tbl into a \code{\link[flowCore]{flowSet}}
#'
#' @param x A tof_tbl.
#'
#' @param group_cols Unquoted names of the columns in `x` that should
#' be used to group cells into separate \code{\link[flowCore]{flowFrame}}s.
#' Supports tidyselect helpers. Defaults to
#' NULL (all cells are written into a single \code{\link[flowCore]{flowFrame}}).
#'
#' @param ... Unused.
#'
#' @return A \code{\link[flowCore]{flowSet}}. Note that all non-numeric
#' columns in `x` will be removed.
#'
#' @rdname as_flowSet
#'
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowCore flowSet
#' @importFrom flowCore phenoData
#' @importFrom flowCore sampleNames
#'
#' @importFrom purrr map
#'
#' @importFrom stringr str_replace
#'
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @examples
#' NULL
as_flowSet.tof_tbl <- function(x, group_cols, ...) {

  if (missing(group_cols)) {
    result <- as_flowFrame(x)

  } else {
    tof_tibble <-
      x |>
      dplyr::select({{group_cols}}, where(tof_is_numeric))

    maxes_and_mins <-
      tof_tibble |>
      dplyr::summarize(
        dplyr::across(
          -{{group_cols}},
          .fns =
            list(max = ~ max(.x, na.rm = TRUE), min = ~ min(.x, na.rm = TRUE)),
          # use the many underscores because it's unlikely this will come up
          # in column names on their own
          .names = "{.col}_____{.fn}"
        )
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = c("antigen", "value_type"),
        values_to = "value",
        names_sep = "_____"
      )  |>
      tidyr::pivot_wider(
        names_from = "value_type",
        values_from = "value"
      )

    # extract the names of all non-grouping columns to be saved to the .fcs file
    data_cols <- maxes_and_mins$antigen

    tof_tibble <-
      suppressWarnings(
        tof_tibble |>
          tidyr::nest(.by = {{group_cols}})
       )

    # make the AnnotatedDataFrame flowCore needs
    parameters <-
      make_flowcore_annotated_data_frame(maxes_and_mins = maxes_and_mins)

    tof_tibble <-
      tof_tibble |>
      dplyr::mutate(
        flowFrames =
          purrr::map(
            .x = data,
            ~ flowCore::flowFrame(
              exprs =
                # have to change any instances of "|" in column names to another
                # separator, as "|" has special meaning as an .fcs file delimiter
                as.matrix(
                  dplyr::rename_with(
                    .x,
                    stringr::str_replace,
                    pattern = "\\|",
                    replacement = "_"
                  )
                ),
              parameters = parameters
            )
          )
      ) |>
      dplyr::select(
        {{group_cols}},
        "flowFrames"
      )

    metadata_frame <-
      tof_tibble |>
      dplyr::select({{group_cols}}) |>
      as.data.frame()

    # store group_cols metadata in an annotated data frame for the flowSet
    row.names(metadata_frame) <- paste0('Sample_', 1:nrow(metadata_frame))
    annotated_metadata_frame <- as(metadata_frame, "AnnotatedDataFrame")

    result <- flowCore::flowSet(tof_tibble$flowFrames)
    flowCore::sampleNames(result) <- paste0('Sample_', 1:length(result))

    flowCore::phenoData(result) <- annotated_metadata_frame

  }

  return(result)
}




#' Make the AnnotatedDataFrame needed for the flowFrame class
#'
#' @param maxes_and_mins a data.frame containing information about the max
#' and min values of each channel to be saved in the flowFrame.
#'
#' @return An AnnotatedDataFrame.
#'
#' @importFrom dplyr transmute
#'
#' @importFrom methods new
#'
#' @importFrom stringr str_replace
#' @importFrom stringr str_c
#'
#' @examples
#' NULL
make_flowcore_annotated_data_frame <- function(maxes_and_mins) {
  fcs_varMetadata <-
    data.frame(
      labelDescription =
        c(
          "Name of Parameter",
          "Description of Parameter",
          "Range of Parameter",
          "Minimum Parameter Value after Transformation",
          "Maximum Parameter Value after Transformation"
        )
    )

  fcs_data <-
    maxes_and_mins |>
    dplyr::transmute(
      # have to change any instances of "|" in column names to another
      # separator, as "|" has special meaning as an .fcs file delimiter
      name = stringr::str_replace(.data$antigen, "\\|", "_"),
      desc = stringr::str_replace(.data$antigen, "\\|", "_"),
      range = max - min,
      minRange = min,
      maxRange = max
    ) |>
    as.data.frame()

  row.names(fcs_data) <- stringr::str_c("$", "P", 1:nrow(fcs_data))

  # make the AnnotatedDataFrame
  parameters <-
    methods::new(
      "AnnotatedDataFrame",
      data = fcs_data,
      varMetadata = fcs_varMetadata
    )

  return(parameters)
}

