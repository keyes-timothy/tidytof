# benchmarking_script.R
#
# This script creates a function that will perform {tidytof} speed, memory, and
# code style benchmarking using user-specified parameters (i.e. which components of the
# benchmarking to run, how many antigen channels and compute repetitions should be used
# for each workflow.

benchmark_master <-
  function(
    set_up_csvs = FALSE,
    run_io = FALSE,
    run_downsample = FALSE,
    run_preprocess = FALSE,
    run_tsne = FALSE,
    run_pca = FALSE,
    run_umap = FALSE,
    run_cluster = FALSE,
    run_extract = FALSE,
    run_memory = FALSE,
    run_style = FALSE,
    num_repeats = 10,
    num_channels = 10,
    base_path,
    verbose = FALSE
  ) {

    # handle file paths
    ddpr_path <-
      file.path(base_path, "DDPR_Data")
    csv_path <-
      file.path(base_path, "DDPR_Data_csv")
    mini_path <-
      file.path(base_path, "DDPR_Data_mini")
    ddpr_files <-
      dir(ddpr_path, full.names = TRUE)

    # read in information about the dataset being used
    ddpr_data_info <-
      read_rds(
        file = here::here("manuscript", "benchmarking", "ddpr_data_info.rds")
      )

    # Find the largest cohort of DDPR files with the same panel
    largest_cohort_names <-
      ddpr_data_info %>%
      filter(num_files == 288L) %>%
      pull(file_names) %>%
      pluck(1)

    largest_cohort_paths <-
      file.path(ddpr_path, largest_cohort_names)

    if (
      any(
        run_io, run_cluster, run_downsample, run_extract,
        run_pca, run_tsne, run_umap, run_preprocess, run_memory
      )
    ) {
      if (verbose) {message("Reading in data")}

      # put together a small benchmarking dataset that includes a set of 1 to
      # 20 FCS files from the largest ddpr cohort
      ddpr_datasets_mini <-
        tibble(
          num_files = 1:20,
          file_paths = map(.x = num_files, .f = ~ largest_cohort_paths[1:.x])
        ) %>%
        mutate(
          tidytof_memory = map_dbl(.x = file_paths, .f = get_memory, mode = "tidytof"),
          flowcore_memory = map_dbl(.x = file_paths, .f = get_memory, mode = "flowcore"),
          immunocluster_memory =
            map_dbl(.x = file_paths, .f = get_memory, mode = "immunocluster"),
          cytofkit_memory =
            map_dbl(.x = file_paths, .f = get_memory, mode = "cytofkit"),
          spectre_memory =
            map_dbl(.x = file_paths, .f = get_memory, mode = "spectre"),
          num_cells = map_int(.x = file_paths, .f = ~nrow(tof_read_data(.x)))
        )

      ddpr_spectre <-
        read_data_spectre(file_names = largest_cohort_paths[1:20], file_type = ".fcs")

      ddpr_spectre_mini <-
        tibble(
          num_files = ddpr_datasets_mini$num_files,
          file_names =
            map(
              .x = num_files,
              .f = ~ str_remove(largest_cohort_names, pattern = ".fcs")[1:.x]
            ),
          data_table =
            map(
              .x = file_names,
              .f = ~ ddpr_spectre[FileName %in% .x, ]
            )
        )

      # if needed, create CSV versions of the files in ddpr_datasets_mini
      if (set_up_csvs) {
        ddpr_large <-
          largest_cohort_paths %>%
          tof_read_data()

        ddpr_large %>%
          mutate(file_name = str_remove(string = file_name, pattern = "\\.fcs")) %>%
          tof_write_data(
            group_cols = file_name,
            out_path = csv_path,
            format = "csv"
          )
      }

      # set up pre-read data sets for downsampling, dimensionality reduction,
      # feature extraction, etc.

      ## every file in the 20-FCS dataset is a row, with single-cell data stored
      ## in the "data" column
      ddpr_mini_tibble <-
        tof_read_data(ddpr_datasets_mini$file_paths[[nrow(ddpr_datasets_mini)]]) %>%
        group_by(file_name) %>%
        nest() %>%
        ungroup()

      ## same as above, but in flowSet form (flowSet of 20 experiments)
      ddpr_mini_flowset <-
        flowCore::read.flowSet(
          files = ddpr_datasets_mini$file_paths[[nrow(ddpr_datasets_mini)]]
        )

      # same as above, but in SCE form (SCE of 20 experiments)
      ddpr_mini_sce <-
        read_data_immunocluster(
          file_names = ddpr_datasets_mini$file_paths[[nrow(ddpr_datasets_mini)]],
          group = basename(ddpr_datasets_mini$file_paths[[nrow(ddpr_datasets_mini)]])
        )

      # store each benchmark iteration as a single row in a tibble.
      # Each row is an iteration. The "tibbles" column contains a single tibble
      # with all the cells for that iteration. The "flowSets" column contains a single
      # flowSet with all the cells for that iteration.
      ddpr_data_mini <-
        ddpr_datasets_mini %>%
        mutate(
          tibbles =
            map(
              .x = 1:nrow(ddpr_datasets_mini),
              .f = ~ ddpr_mini_tibble %>%
                slice_head(n = .x) %>%
                unnest(cols = data)
            ),
          flowSets =
            map(
              .x = 1:nrow(ddpr_datasets_mini),
              .f = ~ddpr_mini_flowset[1:.x]
            ),
          data_tables = ddpr_spectre_mini$data_table,
          SCEs =
            map(
              .x = 1:nrow(ddpr_datasets_mini),
              .f = ~
                ddpr_mini_sce[
                  ,
                  ddpr_mini_sce@metadata$group %in% unique(ddpr_mini_sce@metadata$group)[1:.x]
              ]
            )
        )
    }

    # handle the file paths for saving/loading benchmarking results
    in_out_paths <-
      tibble(
        process =
          c("fcs", "csv", "downsampling", "preprocessing", "tsne", "pca", "umap",  "flowsom", "extraction"),
        pattern =
          paste0(
            c(
              "io_benchmarking_speed_",
              "io_benchmarking_speed_csv_",
              "downsample_benchmarking_speed_",
              "preprocess_benchmarking_speed_",
              "tsne_benchmarking_speed_",
              "pca_benchmarking_speed_",
              "umap_benchmarking_speed_",
              "flowsom_benchmarking_speed_",
              "extract_benchmarking_speed_"
            ),
            "[:digit:]+"
          ),
        largest_size =
          map_chr(
            .x = pattern,
            .f = ~
              here::here("manuscript", "benchmarking") %>%
              dir() %>%
              str_extract(.x) %>%
              na.omit() %>%
              str_extract("[:digit:]+$") %>%
              unique() %>%
              sort() %>%
              pluck(1),
          ),
        file_name =
          paste0(
            str_replace(
              string = pattern, pattern = "\\[.+\\]\\+",
              as.character(num_repeats)
            ),
            ".rds"
          ),
        alternate_name =
          paste0(
            str_replace(string = pattern, pattern = "\\[.+\\]\\+", largest_size),
            ".rds"
          )
      )

    # benchmarking calculations ------------------------------------------------

    # reading files
    if (run_io) {
      if (verbose) {message("benchmark io")}
      ### Benchmarking function
      benchmark_io <-
        function(file_paths) {
          microbenchmark(
            read_data_tidytof(file_names = file_paths),
            read_data_base(file_names = file_paths),
            read_data_flowcore(file_names = file_paths),
            read_data_cytofkit(file_names = file_paths),
            read_data_immunocluster(file_names = file_paths, group = file_paths),
            read_data_spectre(file_names = file_paths, file_type = ".fcs"),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      benchmark_io_csv <-
        function(file_paths) {
          microbenchmark(
            read_data_tidytof(file_names = file_paths),
            read_data_base(file_names = file_paths),
            read_data_spectre(file_names = file_paths, file_type = ".csv"),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      # file paths
      file_name <-
        in_out_paths %>%
        filter(process == "fcs") %>%
        pull(file_name)

      file_name_csv <-
        in_out_paths %>%
        filter(process == "csv") %>%
        pull(file_name)

      ### Perform benchmarking
      io_benchmarking <-
        ddpr_datasets_mini %>%
        #slice_head(n = 2L) %>%
        mutate(io_benchmarks = map(.x = file_paths, .f = benchmark_io))

      io_benchmarking_csv <-
        ddpr_datasets_mini %>%
        mutate(
          file_paths =
            map(.x = file_paths, .f = str_replace, pattern = "\\.fcs", replacement = ".csv") %>%
            map(.f = str_replace, pattern = "DDPR_Data", replacement = "DDPR_Data_csv")
        ) %>%
        mutate(io_benchmarks = map(.x = file_paths, .f = benchmark_io_csv))

      # save benchmarking output as .rds files
      io_benchmarking %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name))

      io_benchmarking_csv %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name_csv))
      }

    # downsampling
    if (run_downsample) {
      if (verbose) {message("benchmark downsampling")}
      ### Benchmarking function
      benchmark_downsample <-
        function(data_frame, flowset, file_names, data_table) {
          microbenchmark(
            downsample_tidytof(data_frame),
            downsample_base(data_frame),
            downsample_flowcore(flowset),
            downsample_cytofkit(file_names),
            downsample_immunocluster(file_names = file_names, group = file_names, num_cells = 200L),
            downsample_spectre(data_table = data_table, num_cells = 200L),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name <-
        in_out_paths %>%
        filter(process == "downsampling") %>%
        pull(file_name)

      downsample_benchmarking <-
        ddpr_data_mini %>%
        mutate(
          downsample_benchmarks =
            pmap(
              .l = list(tibbles, flowSets, file_paths, data_tables),
              .f = benchmark_downsample
            )
        ) %>%
        select(-tibbles, -flowSets, -data_tables, -SCEs)

      # save results as .rds
      downsample_benchmarking %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name))
    }

    # preprocessing
    if (run_preprocess) {
      if (verbose) {message("benchmark preprocessing")}

      benchmark_preprocess <-
        function(data_frame, flowset, file_names, data_table) {
          microbenchmark(
            preprocess_tidytof(data_frame),
            preprocess_base(data_frame),
            preprocess_flowcore(flowset),
            preprocess_cytofkit(file_names),
            preprocess_immunocluster(file_names, group = file_names),
            preprocess_spectre(data_table),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name <-
        in_out_paths %>%
        filter(process == "preprocessing") %>%
        pull(file_name)

      ### Perform benchmarking
      preprocess_benchmarking <-
        ddpr_data_mini %>%
        mutate(
          preprocess_benchmarks =
            pmap(
              .l = list(tibbles, flowSets, file_paths, data_tables),
              .f = benchmark_preprocess
            )
          ) %>%
        select(-tibbles, -flowSets, -data_tables, -SCEs)

      # save benchmarking output as an .rds
      preprocess_benchmarking %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name))
    }

    # set up dimensionality reduction data
    if (any(run_pca, run_tsne, run_umap)) {

      ## pca data
      pca_tibbles <-
        tibble(
          sample_cells = (seq(1000, 10^6, length.out = 10)),
          data_frames =
            map(
              .x = sample_cells,
              .f = ~ unnest(ddpr_mini_tibble, cols = data) %>%
                tof_downsample_constant(num_cells = .x) %>%
                tof_preprocess(undo_noise = FALSE)
            )
        )

      flowset_path <- file.path(base_path, "tidytof_pca_fcs_files")

      pca_tibbles %>%
        select(sample_cells, data_frames) %>%
        unnest(cols = data_frames) %>%
        mutate(file_name = str_remove(file_name, "\\.fcs")) %>%
        tof_write_data(
          group_cols = c(sample_cells, file_name),
          out_path = flowset_path,
          format = "fcs"
        )

      pca_flowsets <-
        tibble(
          num_cells = pca_tibbles$sample_cells
        ) %>%
        mutate(
          file_paths =
            map(.x = as.character(num_cells), .f = ~ dir(flowset_path, pattern = .x, full.names = TRUE))
        )

      pca_flowsets$file_paths[nrow(pca_flowsets)] <-
        list(dir(flowset_path, pattern = "1e+", full.names = TRUE))

      pca_flowsets <-
        pca_flowsets %>%
        mutate(
          flowSets = map(.x = file_paths, .f = read.flowSet)
        )

      pca_tibbles$pca_flowSets <- pca_flowsets$flowSets

      # tsne
      tsne_tibbles <-
        tibble(
          sample_cells = (seq(1000, 10^4, length.out = 10)),
          data_frames =
            map(
              .x = sample_cells,
              .f = ~ unnest(ddpr_mini_tibble, cols = data) %>%
                tof_downsample_constant(num_cells = .x) %>%
                tof_preprocess(undo_noise = FALSE)
            )
        )

      flowset_path <- file.path(base_path, "tidytof_tsne_fcs_files")

      tsne_tibbles %>%
        select(sample_cells, data_frames) %>%
        unnest(cols = data_frames) %>%
        mutate(file_name = str_remove(file_name, "\\.fcs")) %>%
        tof_write_data(
          group_cols = c(sample_cells, file_name),
          out_path = flowset_path,
          format = "fcs"
        )

      tsne_flowsets <-
        tibble(
          num_cells = tsne_tibbles$sample_cells
        ) %>%
        mutate(
          file_paths =
            map(.x = paste0(as.character(num_cells), "_"), .f = ~ dir(flowset_path, pattern = .x, full.names = TRUE))
        )

      tsne_flowsets <-
        tsne_flowsets %>%
        mutate(
          flowSets = map(.x = file_paths, .f = read.flowSet)
        )

      tsne_sces <-
        map(
          .x = tsne_flowsets$file_paths,
          .f = ~ read_data_immunocluster(file_names = .x, group = .x)
        )

      new_dir <-
        file.path(
          tempdir(),
          "tsne_files"
        )
      dir.create(path = new_dir, recursive = TRUE)

      tsne_data_tables <-
        map2(
          .x = tsne_flowsets$file_paths,
          .y = tsne_flowsets$num_cells,
          .f = function(.x, .y) {
            current_dir <- file.path(new_dir, .y)
            dir.create(current_dir)
            for (i in 1:length(.x)) {
              file_name <- .x[[i]]
              base_name <- base::basename(.x[[i]])
              file.copy(from = file_name, to = file.path(current_dir, base_name))
            }
            input_data <-
              Spectre::read.files(
                file.loc = current_dir,
                file.type = ".fcs",
                do.embed.file.names = TRUE
              )
            result <- Spectre::do.merge.files(input_data)
            return(result)
          }
        )

      tsne_tibbles$tsne_flowSets <- tsne_flowsets$flowSets
      tsne_tibbles$data_tables <- tsne_data_tables
      tsne_tibbles$SCEs <- tsne_sces
    }

    # tsne
    if (run_tsne) {
      if (verbose) {message("benchmark tsne")}

      ### Benchmarking function
      benchmark_tsne <-
        function(data_frame, flowset, data_table, SCE) {
          microbenchmark(
            tsne_tidytof(data_frame),
            tsne_base(data_frame),
            tsne_flowcore(flowset),
            tsne_cytofkit(data_frame),
            tsne_spectre(data_table),
            tsne_immunocluster(SCE),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name_tsne <-
        in_out_paths %>%
        filter(process == "tsne") %>%
        pull(file_name)

      ### Perform benchmarking
      tsne_benchmark <-
        tsne_tibbles %>%
        mutate(
          tsne_benchmarks =
            pmap(
              .l = list(data_frames, tsne_flowSets, data_tables, SCEs),
              .f = benchmark_tsne
            )
        ) %>%
        select(-data_frames, -tsne_flowSets, -data_tables, -SCEs) %>%
        dplyr::rename(num_cells = sample_cells)

      # save benchmarking output as an .rds
      tsne_benchmark %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name_tsne))

    }

    # pca
    if (run_pca) {
      if (verbose) {message("benchmark pca")}

      benchmark_pca <-
        function(data_frame, flowset) {
          microbenchmark(
            pca_tidytof(data_frame),
            pca_base(data_frame),
            pca_flowcore(flowset),
            pca_cytofkit(data_frame),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name_pca <-
        in_out_paths %>%
        filter(process == "pca") %>%
        pull(file_name)

      # pca benchmarking
      pca_benchmark <-
        pca_tibbles %>%
        mutate(
          pca_benchmarks = map2(.x = data_frames, .y = pca_flowSets, .f = benchmark_pca)
        ) %>%
        select(-data_frames, -pca_flowSets) %>%
        rename(num_cells = sample_cells)

      # save result as .rds file
      pca_benchmark %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name_pca))
    }

    # umap
    if (run_umap) {
      if (verbose) {message("benchmark umap")}

      benchmark_umap <-
        function(data_frame, flowset, data_table, SCE) {
          microbenchmark(
            umap_tidytof(data_frame),
            #umap_base(data_frame),
            #umap_flowcore(flowset),
            umap_spectre(data_table),
            umap_immunocluster(SCE),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name_umap <-
        in_out_paths %>%
        filter(process == "umap") %>%
        pull(file_name)

      # umap benchmarking
      umap_benchmark <-
        tsne_tibbles %>%
        mutate(
          umap_benchmarks =
            pmap(
              .l = list(data_frames, tsne_flowSets, data_tables, SCEs),
              .f = benchmark_umap
            )
        ) %>%
        select(-data_frames, -tsne_flowSets, -data_tables, -SCEs) %>%
        dplyr::rename(num_cells = sample_cells)

      # save benchmarking output as an .rds
      umap_benchmark %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name_umap))

    }

    # clustering
    if (run_cluster) {
      if (verbose) {message("benchmark clustering")}

      benchmark_clustering <-
        function(file_paths) {
          microbenchmark(
            flowsom_tidytof(file_names = file_paths),
            flowsom_base(file_names = file_paths),
            flowsom_cytofkit(file_names = file_paths),
            flowsom_immunocluster(file_paths, group = file_paths),
            flowsom_spectre(file_names = file_paths),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name <-
        in_out_paths %>%
        filter(process == "flowsom") %>%
        pull(file_name)

      ### perform benchmarking
      cluster_benchmarking <-
        ddpr_datasets_mini %>%
        mutate(
          cluster_benchmarks = map(.x = file_paths, .f = benchmark_clustering)
        )

      # save result as an .rds file
      cluster_benchmarking %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name))
    }

    if (run_extract) {
      if (verbose) {message("benchmark feature extraction")}

      benchmark_extract <-
        function(data_frame, data_table, sce) {
          microbenchmark(
            extract_tidytof(data_frame = data_frame),
            #extract_base(data_frame = data_frame),
            extract_immunocluster(sce),
            extract_spectre(data_table),
            times = num_repeats,
            unit = "s" # seconds
          ) %>%
            as_tibble() %>%
            transmute(
              engine =
                case_when(
                  str_detect(expr, "tidytof")  ~ "tidytof",
                  str_detect(expr, "flowcore") ~ "flowcore",
                  str_detect(expr, "cytofkit") ~ "cytofkit",
                  str_detect(expr, "immunocluster") ~ "immunocluster",
                  str_detect(expr, "spectre") ~ "spectre",
                  TRUE                         ~ "base"
                ),
              time = time / (10^9)
            )
        }

      file_name <-
        in_out_paths %>%
        filter(process == "extraction") %>%
        pull(file_name)

      # perform benchmarking
      extract_benchmarking <-
        ddpr_data_mini %>%
        mutate(
          tibbles =
            map(
              .x = tibbles,
              .f = ~
                mutate(
                  .x,
                  my_cluster =
                    sample(c("a", "b", "c", "d"), size = nrow(.x), replace = TRUE)
                )
            ),
          data_tables = map(
            .x = data_tables,
            .f = ~
              .x[,
                 my_cluster :=
                   sample(c("a", "b", "c", "d"), size = nrow(.x), replace = TRUE)
              ]
          ),
          SCEs =
            map2(
              .x = SCEs,
              .y = file_paths,
              .f = function(sce, file_names) {
                metadata <-
                  sce@metadata %>%
                  dplyr::filter(file %in% file_names) %>%
                  mutate(
                    my_cluster =
                      sample(c("a", "b", "c", "d"), size = ncol(sce), replace = TRUE)
                  )
                sce@metadata <- metadata
                return(sce)
              }
            )
        ) %>%
        mutate(
          extract_benchmarks =
            pmap(
              .l = list(tibbles, data_tables, SCEs),
              .f = benchmark_extract
            )
        ) %>%
        dplyr::select(-tibbles, -flowSets, -data_tables, -SCEs)

      # save result as an .rds file
      extract_benchmarking %>%
        write_rds(file = here::here("manuscript", "benchmarking", file_name))
      }

    # memory
    if (run_memory) {
      if (verbose) {message("benchmark memory")}

      memory_tibble <-
        ddpr_datasets_mini %>%
        select(num_files, num_cells, ends_with("_memory"))

      memory_tibble %>%
        write_rds(file = here::here("manuscript", "benchmarking", "memory.rds"))
    }

    # style
    if (run_style) {
      if (verbose) {message("benchmark style")}

      prefixes <-
        c("read_data", "preprocess", "downsample", "pca", "tsne", "umap", "flowsom", "extract")

      analyses <-
        c("reading files", "preprocessing", "downsampling", "pca", "tsne", "umap", "clustering", "feature extraction")

      base_functions <-
        mget(x = paste0(prefixes, "_base"), envir = globalenv())

      tidytof_functions <-
        mget(x = paste0(prefixes, "_tidytof"), envir = globalenv())

      flowcore_functions <-
        mget(x = paste0(prefixes, "_flowcore"), ifnotfound = NA, envir = globalenv())

      cytofkit_functions <-
        mget(x = paste0(prefixes, "_cytofkit"), ifnotfound = NA, envir = globalenv())

      immunocluster_functions <-
        mget(x = paste0(prefixes, "_immunocluster"), ifnotfound = NA, envir = globalenv())

      spectre_functions <-
        mget(x = paste0(prefixes, "_spectre"), ifnotfound = NA, envir = globalenv())

      # make tibble
      code_tibble <-
        tibble(
          analysis = factor(analyses, levels = analyses),
          base_lines = map_int(.x = base_functions, .f = get_lines),
          tidytof_lines = map_int(.x = tidytof_functions, .f = get_lines),
          flowcore_lines = map_int(.x = flowcore_functions, .f = get_lines),
          immunocluster_lines = map_int(.x = immunocluster_functions, .f = get_lines),
          cytofkit_lines = map_int(.x = cytofkit_functions, .f = get_lines),
          spectre_lines = map_int(.x = spectre_functions, .f = get_lines),
          base_variables = map_int(.x = base_functions, .f = get_assignments),
          tidytof_variables = map_int(.x = tidytof_functions, .f = get_assignments),
          flowcore_variables = map_int(.x = flowcore_functions, .f = get_assignments),
          cytofkit_variables = map_int(.x = cytofkit_functions, .f = get_assignments),
          immunocluster_variables = map_int(.x = immunocluster_functions, .f = get_assignments),
          spectre_variables = map_int(.x = spectre_functions, .f = get_assignments)
        ) %>%
        pivot_longer(
          cols = -analysis,
          names_to = c("engine", "measure"),
          values_to = "number",
          names_sep = "_"
        )

      code_tibble %>%
        write_rds(file = here::here("manuscript", "benchmarking", "code_style.rds"))
    }

  }
