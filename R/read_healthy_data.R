# Description

# Author: Timothy Keyes
# Version: 2020-02-17

# Libraries
library(flowCore)
library(tidyverse)

# Parameters

input_path <- file.path("~", "Desktop", "new_healthies")
out_data <- here::here("c01-own", "data")
marker_path <- here::here("c01-own", "docs", "ALL_panel.csv")

#===============================================================================

#read in names of markers in the dataset 
marker_names <- 
  marker_path %>% 
  read_csv() %>% 
  mutate(Metal = str_replace_all(Metal, "[()]", ""))

str_extract(string = temp, pattern = "\\d{2}_([:alpha:]+_?)+(\\d?)+")

my_data <- 
  input_path %>%
  list.files(path = ., full.names = TRUE) %>% 
  str_extract(string = ., pattern = "\\d{2}_([:alpha:]+_?)+(\\d?)+") %>% 
  str_sub(start = 4L) %>% 
  tibble(population = .) %>% 
  transmute(
    file_name = list.files(path = input_path, full.names = TRUE),
    population, 
    data = 
      map(
        file_name, 
        ~ 
          read.FCS(
            filename = ., 
            transformation = FALSE, 
            truncate_max_range = FALSE
          ) %>% 
          flowCore::exprs() %>% 
          as_tibble()
      )
  )

col_names <- 
  map(my_data$data, colnames) %>% 
  unlist() %>% 
  str_replace_all(pattern = "[()]", replacement = "") %>% 
  table() %>% 
  enframe() %>% 
  arrange(value)

#may have to change this data structure to use regular expressions
lookup_table <-
  setNames(object = marker_names$Metal, nm = marker_names$protein)

#lookup_table["CD34"] <- "Sm148Di"

tof_rename <- function(data, lookup_table) { 
  colnames(data) <- str_replace_all(colnames(data), "[()]", "")
  my_lookups <- (lookup_table[which(lookup_table %in% colnames(data))])
  
  
  data %>% 
    select(which(colnames(.) %in% lookup_table)) %>% 
    rename(!!! my_lookups)
}

my_data <- 
  my_data %>% 
  group_by(file_name, population) %>% 
  unnest(cols = data) %>% 
  ungroup() %>% 
  rename(
    Tdt = TdT
  )



write_rds(
  x = my_data, path = file.path(out_data, "population_data.rds"), compress =  "gz"
)


