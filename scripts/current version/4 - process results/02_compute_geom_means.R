# load libraries ####
library(tidyverse)
library(ggplot2)
library(gtable)
library(coda)
library(readr)
library(here)
library(purrr)
library(beepr)

# load data ####
highout <- read_csv(here("results", "highout.csv"))
medout <- read_csv(here("results", "medout.csv"))
lowout <- read_csv(here("results", "lowout.csv"))

# function to compute quantiles of every variable in dataframe ####
p <- c(0.025, 0.5, 0.975)
p_names <- map_chr(p, ~paste0(.x*100, "%"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

# reformat data #####
row.high <- highout %>%
  select(contains("lambda") | contains("sims") | contains("scenario")) %>% 
  group_by(sims, scenario, simscenarios) %>% 
  summarise(across(everything(), funs(!!!p_funs))) %>% # get quantiles
  pivot_longer( # begin reshaping
    cols = contains("%"), 
    names_to = "quantile",
    values_to = "lambda",
  ) %>% 
  mutate(Year = str_extract(quantile, "\\[\\d+\\]")) %>% 
  mutate(Year = str_extract(quantile, "\\d+")) %>% 
  mutate(Quantile = str_extract(quantile, "\\d+\\.?\\d?%")) %>% 
  select(-quantile) %>% 
  select(c(1:3, 6, 5, 4)) %>% 
  arrange(sims, scenario, simscenarios, Quantile, Year) %>% 
  pivot_wider(names_from = Year, values_from = lambda, names_prefix = "Year_") %>% 
  select(c(1:4, "Year_1", "Year_2", "Year_3", "Year_4", "Year_5", "Year_6", 
                "Year_7", "Year_8", "Year_9", "Year_10", "Year_11", "Year_12", 
                "Year_13", "Year_14")) # end reshaping

row.med <- medout %>%
  select(contains("lambda") | contains("sims") | contains("scenario")) %>% 
  group_by(sims, scenario, simscenarios) %>% 
  summarise(across(everything(), funs(!!!p_funs))) %>% # get quantiles
  pivot_longer( # begin reshaping
    cols = contains("%"), 
    names_to = "quantile",
    values_to = "lambda",
  ) %>% 
  mutate(Year = str_extract(quantile, "\\[\\d+\\]")) %>% 
  mutate(Year = str_extract(quantile, "\\d+")) %>% 
  mutate(Quantile = str_extract(quantile, "\\d+\\.?\\d?%")) %>% 
  select(-quantile) %>% 
  select(c(1:3, 6, 5, 4)) %>% 
  arrange(sims, scenario, simscenarios, Quantile, Year) %>% 
  pivot_wider(names_from = Year, values_from = lambda, names_prefix = "Year_") %>% 
  select(c(1:4, "Year_1", "Year_2", "Year_3", "Year_4", "Year_5", "Year_6", 
           "Year_7", "Year_8", "Year_9", "Year_10", "Year_11", "Year_12", 
           "Year_13", "Year_14")) # end reshaping

row.low <- lowout %>%
  select(contains("lambda") | contains("sims") | contains("scenario")) %>% 
  group_by(sims, scenario, simscenarios) %>% 
  summarise(across(everything(), funs(!!!p_funs))) %>% # get quantiles
  pivot_longer( # begin reshaping
    cols = contains("%"), 
    names_to = "quantile",
    values_to = "lambda",
  ) %>% 
  mutate(Year = str_extract(quantile, "\\[\\d+\\]")) %>% 
  mutate(Year = str_extract(quantile, "\\d+")) %>% 
  mutate(Quantile = str_extract(quantile, "\\d+\\.?\\d?%")) %>% 
  select(-quantile) %>% 
  select(c(1:3, 6, 5, 4)) %>% 
  arrange(sims, scenario, simscenarios, Quantile, Year) %>% 
  pivot_wider(names_from = Year, values_from = lambda, names_prefix = "Year_") %>% 
  select(c(1:4, "Year_1", "Year_2", "Year_3", "Year_4", "Year_5", "Year_6", 
           "Year_7", "Year_8", "Year_9", "Year_10", "Year_11", "Year_12", 
           "Year_13", "Year_14")) # end reshaping

# create new variables for geometric mean by year
for (i in 1:14) {
  row.low <- row.low %>%
    mutate("geomean.{i}" :=  NA_real_)
  row.med <- row.med %>%
    mutate("geomean.{i}" :=  NA_real_)
  row.high <- row.high %>%
    mutate("geomean.{i}" :=  NA_real_)
}

beep(sound = 8)
rm("highout", "medout", "lowout")

# compute geometric means #######

# slow
for(i in 1:dim(row.low)[1]) {
  print(paste("row", i))
  for(j in 1:((ncol(row.low) - 4 )/2)) {
    row.low[i, ((ncol(row.low) + 4)/2 + j)] <- exp(mean(unlist(log(row.low[i, 4 + 1:j]))))
  }
}

# slow
for(i in 1:dim(row.med)[1]) {
  print(paste("row", i))
  for(j in 1:((ncol(row.med) - 4 )/2)) {
    row.med[i, ((ncol(row.med) + 4)/2 + j)] <- exp(mean(unlist(log(row.med[i, 4 + 1:j]))))
  }
}

# slow
for(i in 1:dim(row.high)[1]) {
  print(paste("row", i))
  for(j in 1:((ncol(row.high) - 4 )/2)) {
    row.high[i, ((ncol(row.high) + 4)/2 + j)] <- exp(mean(unlist(log(row.high[i, 4 + 1:j]))))
  }
}

# save objects
write_csv(row.low, here("results", "row_low.csv"))
write_csv(row.med, here("results", "row_med.csv"))
write_csv(row.high, here("results", "row_high.csv"))
