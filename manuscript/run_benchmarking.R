# run_benchmarking.R
#
# This script is used to render the rmarkdown document named
# "tidytof_performance_benchmarking.Rmd".
#
rmarkdown::render(
  input = "tidytof_performance_benchmarking.Rmd",
  output_format = "pdf_document",
  output_file = "tidytof_benchmarking",
  output_dir = here::here("manuscript", "benchmarking")
)
