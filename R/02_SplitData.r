library(tidyverse)
library(tidymodels)
library(stringr)

suppressMessages({
  full_data <- read_csv("./data/full_data.csv") |> 
    mutate(
      Testing = str_starts(path, "1024455") | str_starts(path, "658")
    )
})

full_data_train <- full_data |> 
  filter(!Testing) |> 
  select(-Testing)

full_data_test <- full_data |> 
  filter(Testing) |> 
  select(-Testing)

training_min_class_num <- full_data_train |> 
  pull(GOF) |> 
  table() |> 
  min()

full_data_train_balanced <- full_data_train |> 
  group_by(GOF) |> 
  slice_sample(n = training_min_class_num) |> 
  ungroup()

testing_min_class_num <- full_data_test |> 
  pull(GOF) |> 
  table() |> 
  min()

full_data_test_balanced <- full_data_test |> 
  group_by(GOF) |> 
  slice_sample(n = testing_min_class_num) |> 
  ungroup() |> 
  mutate(
    splitter = 1:n(),
    splitter = splitter %% 2,
    splitter = splitter |> as.logical()
  )

full_data_train_balanced |> 
  write_csv("./data/data_train_full.csv")

full_data_test_balanced |> 
  filter(splitter) |> 
  select(-splitter) |> 
  write_csv("./data/data_test_full")

full_data_test_balanced |> 
  filter(!splitter) |> 
  select(-splitter) |> 
  write_csv("./data/data_calib_full")
